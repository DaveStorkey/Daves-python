#! /usr/bin/env python

'''
Routine to plot vertical heat or salt budget.

NB. Assumes that the input file doesn't have multiple times.

To do: 
  1. Option to remove seasonal variability from the "eddy" trends.

@author: Dave Storkey
'''

import matplotlib
# Hack for running in batch mode on SPICE. 
# Note this disables plt.show()
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import numpy.ma as ma

def vert_tracer_budget(filename,coordsfile=None,fieldtype=None,fields=None,outfile=None,xmin=None,xmax=None,depthmax=None,
                       totaltrend=None,integralcheck=None,diagnostics=None,nomaskcheck=None):

    # field labels

    field_labels = { 'trd_adv_mean'           : 'advection by mean',
                     'trd_adv_eddy'           : 'advection by eddies',
                     'trd_adv_GM'             : 'advection by Gent-McWilliams bolus velocity',
                     'trd_adv_mean_plus_eddy' : 'sum of advection by mean and eddies',
                     'trd_adv_mean_plus_eddy_plus_GM' : 'sum of advection by mean and eddies and Gent-McWilliams',
                     'trd_totad'              : 'total advection (model diagnostic)',
                     'trd_xad'                : 'x-component of advection',
                     'trd_yad'                : 'y-component of advection',
                     'trd_zad'                : 'z-component of advection',
                     'trd_ldf'                : 'lateral diffusion',
                     'trd_iso_x'              : 'x-component of isopycnal diffusion',
                     'trd_iso_y'              : 'y-component of isopycnal diffusion',
                     'trd_iso_z'              : 'z-component of isopycnal diffusion',
                     'trd_iso_z1'             : 'partial z-component of isopycnal diffusion',
                     'trd_zdf'                : 'total vertical diffusion',
                     'trd_zdfp'               : 'pure vertical diffusion with convection',
                     'trd_zdfp_minus_evd'     : 'pure vertical diffusion without convection',
                     'trd_evd'                : 'EVD convection',
                     'trd_qsr'                : 'penetrating solar heating',
                     'trd_sum'                : 'sum of model trends',
                     'trd_tot'                : 'total trend (model diagnostic)',
                     'model_trend'             : 'total trend (end field minus start field)'
                   }

    # reference density and heat capacity: NEMO 3.6 values!
    rau0 = 1026.0            
    rcp  = 3991.86795711963  

    if fieldtype == 'T' or fieldtype == 't':
        prefix='t'
        names_in_file=['votemper','tn',prefix+fields[0]]
        GMtransport='heattr'
        xlabel = 'temperature trend (degC/year)'
    elif fieldtype == 'S' or fieldtype == 's':
        prefix='s'
        names_in_file=['vosaline','sn',prefix+fields[0]]
        GMtransport='salttr'
        xlabel = 'salinity trend (psu/year)'
    else: 
        raise Exception('Unrecognised (or unspecified) fieldtype. Should be T or S.')

    if fields is None:
        fields=['','','','']

    coords_in= nc.Dataset(coordsfile,'r')
    infile_T = nc.Dataset(filename+'_grid_T.nc','r')
    depth=infile_T.variables['deptht']

    e1t = coords_in.variables['e1t'][:]
    e2t = coords_in.variables['e2t'][:]
    # unmask the e3 arrays
    e3t = infile_T.variables['e3t'][:]; e3t[:].mask = ma.nomask
    # prototype T-grid field for shape and mask
    for fieldname in names_in_file:
        try:
            TSfield = infile_T.variables[fieldname][:]
        except KeyError:
            pass
        else:
            break
    else:
        raise Exception('Could not find prototype T-grid field in file.')

    if not nomaskcheck:
        # check for masking
        if ma.count_masked(TSfield) == 0:
            raise Exception("Error: prototype T-grid field is not mask, so spatial averages will be incorrect.")

    ######## Calculate advective trends from uT type diagnostics if required ########

    if 'trd_adv_mean' in fields or 'trd_adv_eddy' in fields or 'trd_adv_mean_plus_eddy' in fields or 'trd_adv_mean_plus_eddy_plus_GM' in fields: 

        infile_U = nc.Dataset(filename+'_grid_U.nc','r')
        infile_V = nc.Dataset(filename+'_grid_V.nc','r')
        infile_W = nc.Dataset(filename+'_grid_W.nc','r')
        # unmask the e3 arrays
        e3u = infile_U.variables['e3u'][:]; e3u[:].mask = ma.nomask
        e3v = infile_V.variables['e3v'][:]; e3v[:].mask = ma.nomask

        if diagnostics:
            # write 3D diagnostics fields if requested
            indims_T=infile_T.dimensions
            outdata_T = nc.Dataset(filename+"_out_grid_T.nc","w")
            for dim in infile_T.variables[names_in_file[0]].dimensions:
                outdata_T.createDimension(dim, len(indims_T[dim])) 

        # Create <U><T> term and <U'T'> fluxes. Need to horizontally and vertically interpolate T.
        # New arrays inherit mask from temp and are set to zero at unmasked points.
        TSfield_interp   = ma.zeros(TSfield.shape)
        meanU_meanT   = np.zeros(TSfield.shape) # Don't want these fields
        Uprime_Tprime = np.zeros(TSfield.shape) # to be masked.
        trd_adv_mean = ma.zeros(TSfield.shape)
        trd_adv_eddy = ma.zeros(TSfield.shape)
    
        print '>>> calculating x-fluxes and adding divergence to trend'
        Uvel = infile_U.variables['vozocrtx'][:]
        mean_UT = infile_U.variables['ut'][:]
        GM=True
        try: 
            ueiv_TStr = infile_U.variables['ueiv_'+GMtransport+'3d'][:]
            # unmask the array, putting zeroes at masked points
            ueiv_TStr[ueiv_TStr.mask] = 0.0
            ueiv_TStr.mask = ma.nomask
        except KeyError:
            GM=False
        else:
            print '>>>> found GM trends'
            trd_adv_GM = ma.zeros(TSfield.shape)
      
        # Calculate fluxes on internal interfaces between T-cells. 
        # Fluxes on external interfaces are zero.
        for ii in range(TSfield.shape[2]-3):
            TSfield_interp[:,:,ii+1]   = 0.5 * ( TSfield[:,:,ii+1] + TSfield[:,:,ii+2] )
            meanU_meanT[:,:,ii+1]   = Uvel[:,:,ii+1] * TSfield_interp[:,:,ii+1]
            Uprime_Tprime[:,:,ii+1] = mean_UT[:,:,ii+1] - meanU_meanT[:,:,ii+1]
        for ii in range(TSfield.shape[2]-2):
            trd_adv_mean[:,:,ii+1] =  trd_adv_mean[:,:,ii+1] + \
              ( e3u[:,:,ii]*meanU_meanT[:,:,ii] - e3u[:,:,ii+1]*meanU_meanT[:,:,ii+1] ) \
              / ( e1t[:,ii+1]*e3t[:,:,ii+1] )
            trd_adv_eddy[:,:,ii+1] =  trd_adv_eddy[:,:,ii+1] + \
              ( e3u[:,:,ii]*Uprime_Tprime[:,:,ii] - e3u[:,:,ii+1]*Uprime_Tprime[:,:,ii+1] ) \
              / ( e1t[:,ii+1]*e3t[:,:,ii+1] )
        if GM:
            for ii in range(TSfield.shape[2]-2):
                trd_adv_GM[:,:,ii+1] =  trd_adv_GM[:,:,ii+1] + \
                  ( ueiv_TStr[:,:,ii] - ueiv_TStr[:,:,ii+1] ) \
                  / ( rau0*rcp*e1t[:,ii+1]*e2t[:,ii+1]*e3t[:,:,ii+1] )

        print '>>> calculating y-fluxes and adding divergence to trend'
        # re-use arrays to save memory
        Uvel[:] = 0.0
        mean_UT[:] = 0.0
        Uvel = infile_V.variables['vomecrty'][:]
        mean_UT = infile_V.variables['vt'][:]
        if GM:
            ueiv_TStr[:] = 0.0
            ueiv_TStr = infile_V.variables['veiv_'+GMtransport+'3d'][:]
            # unmask the array, putting zeroes at masked points
            ueiv_TStr[ueiv_TStr.mask] = 0.0
            ueiv_TStr.mask = ma.nomask
        TSfield_interp[:] = 0.0
        meanU_meanT[:] = 0.0
        Uprime_Tprime[:] = 0.0
        # Calculate fluxes on internal interfaces between T-cells. 
        # Fluxes on external interfaces are zero.
        for ij in range(TSfield.shape[1]-3):
            TSfield_interp[:,ij+1,:] = 0.5 * ( TSfield[:,ij+1,:] + TSfield[:,ij+2,:] )
            meanU_meanT[:,ij+1,:] = Uvel[:,ij+1,:] * TSfield_interp[:,ij+1,:]
            Uprime_Tprime[:,ij+1,:] = mean_UT[:,ij+1,:] - meanU_meanT[:,ij+1,:]
        for ij in range(TSfield.shape[1]-2):
            trd_adv_mean[:,ij+1,:] =  trd_adv_mean[:,ij+1,:] + \
              ( e3v[:,ij,:]*meanU_meanT[:,ij,:] - e3v[:,ij+1,:]*meanU_meanT[:,ij+1,:] ) \
             / ( e2t[ij+1,:]*e3t[:,ij+1,:] )
            trd_adv_eddy[:,ij+1,:] =  trd_adv_eddy[:,ij+1,:] + \
              ( e3v[:,ij,:]*Uprime_Tprime[:,ij,:] - e3v[:,ij+1,:]*Uprime_Tprime[:,ij+1,:] ) \
              / ( e2t[ij+1,:]*e3t[:,ij+1,:] )
        if GM:
            for ij in range(TSfield.shape[1]-2):
                trd_adv_GM[:,ij+1,:] =  trd_adv_GM[:,ij+1,:] + \
                  ( ueiv_TStr[:,ij,:] - ueiv_TStr[:,ij+1,:] ) \
                  / ( rau0*rcp*e1t[ij+1,:]*e2t[ij+1,:]*e3t[:,ij+1,:] )

        print '>>> calculating z-fluxes'
        # re-use arrays to save memory
        Uvel[:] = 0.0
        mean_UT[:] = 0.0
        Uvel = infile_W.variables['vovecrtz'][:]
        mean_UT = infile_W.variables['wt'][:]
        if GM:
            ueiv_TStr[:] = 0.0
            ueiv_TStr = infile_W.variables['weiv_'+GMtransport+'3d'][:]
            # unmask the array, putting zeroes at masked points
            ueiv_TStr[ueiv_TStr.mask] = 0.0
            ueiv_TStr.mask = ma.nomask
        TSfield_interp[:] = 0.0
        meanU_meanT[:] = 0.0
        Uprime_Tprime[:] = 0.0
        # Calculate fluxes on internal interfaces between T-cells. 
        # Fluxes on external interfaces are zero.
        for ik in range(TSfield.shape[0]-3):
            TSfield_interp[ik+1,:,:] = 0.5 * ( TSfield[ik,:,:] + TSfield[ik+1,:,:] )
            meanU_meanT[ik+1,:,:] = Uvel[ik+1,:,:] * TSfield_interp[ik+1,:,:]
            Uprime_Tprime[ik+1,:,:] = mean_UT[ik+1,:,:] - meanU_meanT[ik+1,:,:]
        for ik in range(TSfield.shape[0]-1):
            trd_adv_mean[ik,:,:] =  trd_adv_mean[ik,:,:] + \
                        ( meanU_meanT[ik+1,:,:] - meanU_meanT[ik,:,:] ) / e3t[ik,:,:]
            trd_adv_eddy[ik,:,:] =  trd_adv_eddy[ik,:,:] + \
                        ( Uprime_Tprime[ik+1,:,:] - Uprime_Tprime[ik,:,:] ) / e3t[ik,:,:]
        if GM:
            for ik in range(TSfield.shape[0]-2):
                trd_adv_GM[ik,:,:] =  trd_adv_GM[ik,:,:] + \
                  ( ueiv_TStr[ik+1,:,:] - ueiv_TStr[ik,:,:] ) \
              / ( rau0*rcp*e1t[:,:]*e2t[:,:]*e3t[ik,:,:] )

        # mask the final 3D trends fields to ensure the area integrals are calculated correctly
        # (shouldn't include the zeroes in the outermost rows/columns in the integral).
        trd_adv_mean.mask = TSfield.mask
        trd_adv_eddy.mask = TSfield.mask
        if GM:
            trd_adv_GM.mask = TSfield.mask

        if diagnostics:
            outdata_T.createVariable(prefix+'trd_adv_mean',datatype='f',dimensions=outdata_T.dimensions.keys(),fill_value=-1.0e+20)
            outdata_T.variables[prefix+'trd_adv_mean'][:] = trd_adv_mean[:]
            outdata_T.createVariable(prefix+'trd_adv_eddy',datatype='f',dimensions=outdata_T.dimensions.keys(),fill_value=-1.0e+20)
            outdata_T.variables[prefix+'trd_adv_eddy'][:] = trd_adv_eddy[:]
            if GM:
                 outdata_T.createVariable(prefix+'trd_adv_GM',datatype='f',dimensions=outdata_T.dimensions.keys(),fill_value=-1.0e+20)
                 outdata_T.variables[prefix+'trd_adv_GM'][:] = trd_adv_GM[:]

    ######## Read in (or initialise) required fields ########
    # not very memory efficient to read in all 3D fields at once...

    trd = {} # will be dictionary of masked arrays
    for fieldname in fields:
        if fieldname == 'trd_adv_mean':
            trd[fieldname] = trd_adv_mean
        elif fieldname == 'trd_adv_eddy':
            trd[fieldname] = trd_adv_eddy
        elif fieldname == 'trd_adv_GM':
            trd[fieldname] = trd_adv_GM
        elif fieldname in ['trd_zdfp_minus_evd','trd_adv_mean_plus_eddy','trd_adv_mean_plus_eddy_plus_GM','trd_sum']:
            pass # calculate these later when all other fields read in
        else:   
            trd[fieldname] = infile_T.variables[prefix+fieldname][:]

    if 'trd_zdfp_minus_evd' in fields:
        for fieldname in ['trd_zdfp','trd_evd']:
            if fieldname not in fields:
                trd[fieldname] = infile_T.variables[prefix+fieldname][:]
        trd['trd_zdfp_minus_evd'] = trd['trd_zdfp'] - trd['trd_evd']

    if 'trd_adv_mean_plus_eddy' in fields:
        trd['trd_adv_mean_plus_eddy'] = trd['trd_adv_mean'] + trd['trd_adv_eddy']
        if diagnostics:
            outdata_T.createVariable(prefix+'trd_adv_mean_plus_eddy',datatype='f',dimensions=outdata_T.dimensions.keys(),fill_value=-1.0e+20)
            outdata_T.variables[prefix+'trd_adv_mean_plus_eddy'][:] = trd['trd_adv_mean_plus_eddy'][:]

    if 'trd_adv_mean_plus_eddy_plus_GM' in fields:
        if GM:
            trd['trd_adv_mean_plus_eddy_plus_GM'] = trd['trd_adv_mean'] + trd['trd_adv_eddy'] + trd['trd_adv_GM']
            if diagnostics:
                outdata_T.createVariable(prefix+'trd_adv_mean_plus_eddy_plus_GM',datatype='f',dimensions=outdata_T.dimensions.keys(),fill_value=-1.0e+20)
                outdata_T.variables[prefix+'trd_adv_mean_plus_eddy_plus_GM'][:] = trd['trd_adv_mean_plus_eddy_plus_GM'][:]
        else:
            print "Can't plot trd_adv_mean_plus_eddy_plus_GM because there aren't any GM trends in the file."

    # need to do this last of course
    if 'trd_sum' in fields:
        fields_to_sum = list(fields) # force it to copy the list not just point to the original
        for field_to_remove in ['trd_sum','trd_tot']:
            if field_to_remove in fields_to_sum:
                fields_to_sum.remove(field_to_remove)
        trd['trd_sum'] = ma.copy(trd[fields_to_sum[0]])
        for fieldname in fields_to_sum[1:]:
            trd['trd_sum'][:] = trd['trd_sum'][:] + trd[fieldname][:]

    if totaltrend is not None:
        totalfile = nc.Dataset(totaltrend,'r')
        for varname in names_in_file:
            try:
                trd['model_trend'] = totalfile.variables[varname][:]
            except KeyError:
                pass
            else:
                break
        else:
            raise Exception('Could not find temperature field in model trend file.')
        fields=fields+['model_trend']

    ######## Horizontal average trends and plot ########

    # Re-mask the e3 fields otherwise uninitialised values will contaminate horizontal averages. 
    e3t[:].mask = TSfield[:].mask

    field_min=None
    field_max=None
    ones3D = ma.ones(TSfield.shape)
    # apply 3D mask which will be inherited by 3D versions of e1t and e2t.
    ones3D.mask = TSfield.mask
    # broadcast 2D scale factor fields to 3D
    e1t_3d = e1t[:] * ones3D[:]
    e2t_3d = e2t[:] * ones3D[:]
    # xy_averaged_area is different on different depth levels because 
    # of the different mask at different depths. 
    xy_averaged_area = ma.average(ma.average(e1t_3d[:]*e2t_3d[:],axis=-1),axis=-1)
    x_averaged_e3t = ma.average(e3t[:,:,:],axis=-1)
    xy_averaged_e3t = ma.average(x_averaged_e3t[:,:],axis=-1)
    xyz_averaged_e3t = ma.average(xy_averaged_e3t[:])

    labels_for_plot=[]
    for fieldname in fields:

        print '>>> horizontally averaging '+fieldname
        x_averaged_field = ma.average(trd[fieldname][:,:,:]*e1t[:,:]*e2t[:,:]*e3t[:,:,:],axis=-1)
        xy_averaged_field = ma.average(x_averaged_field[:,:],axis=-1)    
        if fieldname == 'model_trend':
            # model_trend field already in units of change per year
            xy_averaged_field[:] = ( xy_averaged_field[:] / (xy_averaged_area[:] * xy_averaged_e3t[:]) )
        else:
            # scale from per second to per 365-day year:
            xy_averaged_field[:] = ( xy_averaged_field[:] / (xy_averaged_area[:] * xy_averaged_e3t[:]) ) * 31536000 

        if integralcheck is True:
            xyz_averaged_field = ma.average(xy_averaged_field[:]*xy_averaged_e3t[:])/xyz_averaged_e3t
            print '>>> domain average for '+fieldname+' : ',xyz_averaged_field

        if field_min is None:
            field_min = xy_averaged_field.min()
        else:
            field_min=min(xy_averaged_field.min(),field_min)

        if field_max is None:
            field_max = xy_averaged_field.max()
        else:
            field_max=max(xy_averaged_field.max(),field_max)

        plt.plot(xy_averaged_field[:],depth[:])
        labels_for_plot.append(field_labels[fieldname])
        
    if xmin is None:
        xmin = field_min
    if xmax is None:
        xmax = field_max
    plt.gca().set_xlim([xmin,xmax])

    if depthmax is not None:
        plt.gca().set_ylim([depthmax,0.0])
    else:
        plt.gca().invert_yaxis()

    plt.gca().set_xlabel(xlabel)
    plt.gca().set_ylabel('depth (m)')
    plt.legend(labels_for_plot,loc=3,fontsize='medium')

    if outfile is not None:
        matplotlib.rcParams['font.size'] =8
        plt.savefig(outfile,dpi=200)
    else:
        plt.show()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="name of input file")
    parser.add_argument("-c", "--coordsfile", action="store",dest="coordsfile",default=None,
                    help="name of coordinates file (compulsory)")
    parser.add_argument("-F", "--fieldtype", action="store",dest="fieldtype",default=None,
                    help="type of field: temperature (-F T) or salinity (-F S)")
    parser.add_argument("-f", "--fields", action="store",dest="fields",nargs='+',default=None,
                    help="names of fields to plot excluding leading 's' or 't' from ?trd_ type diagnostics")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",default=None,
                    help="plot to output file - filetype by extension. If unset plot to GUI.")
    parser.add_argument("-x", "--xmin", action="store",dest="xmin",type=float,default=None,
                    help="start of x-axis (min field value to plot)")
    parser.add_argument("-X", "--xmax", action="store",dest="xmax",type=float,default=None,
                    help="end of x-axis (max field value to plot)")
    parser.add_argument("-D", "--depthmax", action="store",dest="depthmax",type=float,default=None,
                    help="maximum depth for plot")
    parser.add_argument("-T", "--totaltrend", action="store",dest="totaltrend",default=None,
                    help="file containing the total model trend")
    parser.add_argument("-I", "--integralcheck", action="store_true",dest="integralcheck",
                    help="check volume integrals of trends")
    parser.add_argument("-a", "--diagnostics", action="store_true",dest="diagnostics",
                    help="write diagnostic fields for calculation of advective trends")
    parser.add_argument("-N", "--nomaskcheck", action="store_true",dest="nomaskcheck",
                    help="turn off check for masking, ie. allow averaging of unmasked fields")

    args = parser.parse_args()

    vert_tracer_budget( args.filename, coordsfile=args.coordsfile, fields=args.fields, outfile=args.outfile, 
                        xmin=args.xmin, xmax=args.xmax, depthmax=args.depthmax, fieldtype=args.fieldtype, 
                        totaltrend=args.totaltrend, integralcheck=args.integralcheck, diagnostics=args.diagnostics,
                        nomaskcheck=args.nomaskcheck )
