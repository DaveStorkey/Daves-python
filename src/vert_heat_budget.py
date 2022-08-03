#! /usr/bin/env python

'''
Routine to plot vertical heat budget.

NB. Assumes that the input file doesn't have multiple times.

To do: 
  1. Generalise it to do salinity trends as well as temperature trends. 
  2. Remove seasonal variability from the "eddy" trends.

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

def vert_heat_budget(filename,coordsfile=None,fields=None,outfile=None,xmin=None,xmax=None,depthmax=None,
                     totaltrend=None,integralcheck=None,diagnostics=None):

    # field labels

    field_labels = { 'ttrd_adv_mean'           : 'advection by mean',
                     'ttrd_adv_eddy'           : 'advection by eddies',
                     'ttrd_adv_GM'             : 'advection by Gent-McWilliams bolus velocity',
                     'ttrd_adv_mean_plus_eddy' : 'sum of advection by mean and eddies',
                     'ttrd_adv_mean_plus_eddy_plus_GM' : 'sum of advection by mean and eddies and Gent-McWilliams',
                     'ttrd_totad'              : 'total advection (model diagnostic)',
                     'ttrd_xad'                : 'x-component of advection',
                     'ttrd_yad'                : 'y-component of advection',
                     'ttrd_zad'                : 'z-component of advection',
                     'ttrd_ldf'                : 'lateral diffusion',
                     'ttrd_iso_x'              : 'x-component of isopycnal diffusion',
                     'ttrd_iso_y'              : 'y-component of isopycnal diffusion',
                     'ttrd_iso_z'              : 'z-component of isopycnal diffusion',
                     'ttrd_iso_z1'             : 'partial z-component of isopycnal diffusion',
                     'ttrd_zdf'                : 'total vertical diffusion',
                     'ttrd_zdfp'               : 'pure vertical diffusion with convection',
                     'ttrd_zdfp_minus_evd'     : 'pure vertical diffusion without convection',
                     'ttrd_evd'                : 'EVD convection',
                     'ttrd_qsr'                : 'penetrating solar heating',
                     'ttrd_sum'                : 'sum of model trends',
                     'ttrd_tot'                : 'total trend (model diagnostic)',
                     'model_trend'             : 'total trend (end field minus start field)'
                   }

    # reference density and heat capacity: NEMO 3.6 values!
    rau0 = 1026.0            
    rcp  = 3991.86795711963  

    if fields is None:
        fields=['','','','']

    coords_in= nc.Dataset(coordsfile,'r')
    infile_T = nc.Dataset(filename+'_grid_T.nc','r')
    depth=infile_T.variables['deptht']

    e1t = coords_in.variables['e1t'][:]
    e2t = coords_in.variables['e2t'][:]
    # unmask the e3 arrays
    e3t = infile_T.variables['e3t'][:]; e3t[:].mask = ma.nomask
    temp = infile_T.variables['votemper'][:]

    ######## Calculate advective trends from uT type diagnostics if required ########

    if 'ttrd_adv_mean' in fields or 'ttrd_adv_eddy' in fields or 'ttrd_adv_mean_plus_eddy' in fields or 'ttrd_adv_mean_plus_eddy_plus_GM' in fields: 

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
            for dim in infile_T.variables['votemper'].dimensions:
                outdata_T.createDimension(dim, len(indims_T[dim])) 

        # Create <U><T> term and <U'T'> fluxes. Need to horizontally and vertically interpolate T.
        # New arrays inherit mask from temp and are set to zero at unmasked points.
        temp_interp   = ma.zeros(temp.shape)
        meanU_meanT   = np.zeros(temp.shape) # Don't want these fields
        Uprime_Tprime = np.zeros(temp.shape) # to be masked.
        ttrd_adv_mean = ma.zeros(temp.shape)
        ttrd_adv_eddy = ma.zeros(temp.shape)
    
        print '>>> calculating x-fluxes and adding divergence to trend'
        Uvel = infile_U.variables['vozocrtx'][:]
        mean_UT = infile_U.variables['ut'][:]
        GM=True
        try: 
            ueiv_heattr = infile_U.variables['ueiv_heattr3d'][:]
            # unmask the array, putting zeroes at masked points
            ueiv_heattr[ueiv_heattr.mask] = 0.0
            ueiv_heattr.mask = ma.nomask
        except KeyError:
            GM=False
        else:
            print '>>>> found GM trends'
            ttrd_adv_GM = ma.zeros(temp.shape)
      
        # Calculate fluxes on internal interfaces between T-cells. 
        # Fluxes on external interfaces are zero.
        for ii in range(temp.shape[2]-3):
            temp_interp[:,:,ii+1]   = 0.5 * ( temp[:,:,ii+1] + temp[:,:,ii+2] )
            meanU_meanT[:,:,ii+1]   = Uvel[:,:,ii+1] * temp_interp[:,:,ii+1]
            Uprime_Tprime[:,:,ii+1] = mean_UT[:,:,ii+1] - meanU_meanT[:,:,ii+1]
        for ii in range(temp.shape[2]-2):
            ttrd_adv_mean[:,:,ii+1] =  ttrd_adv_mean[:,:,ii+1] + \
              ( e3u[:,:,ii]*meanU_meanT[:,:,ii] - e3u[:,:,ii+1]*meanU_meanT[:,:,ii+1] ) \
              / ( e1t[:,ii+1]*e3t[:,:,ii+1] )
            ttrd_adv_eddy[:,:,ii+1] =  ttrd_adv_eddy[:,:,ii+1] + \
              ( e3u[:,:,ii]*Uprime_Tprime[:,:,ii] - e3u[:,:,ii+1]*Uprime_Tprime[:,:,ii+1] ) \
              / ( e1t[:,ii+1]*e3t[:,:,ii+1] )
        if GM:
            for ii in range(temp.shape[2]-2):
                ttrd_adv_GM[:,:,ii+1] =  ttrd_adv_GM[:,:,ii+1] + \
                  ( ueiv_heattr[:,:,ii] - ueiv_heattr[:,:,ii+1] ) \
                  / ( rau0*rcp*e1t[:,ii+1]*e2t[:,ii+1]*e3t[:,:,ii+1] )

        print '>>> calculating y-fluxes and adding divergence to trend'
        # re-use arrays to save memory
        Uvel[:] = 0.0
        mean_UT[:] = 0.0
        Uvel = infile_V.variables['vomecrty'][:]
        mean_UT = infile_V.variables['vt'][:]
        if GM:
            ueiv_heattr[:] = 0.0
            ueiv_heattr = infile_V.variables['veiv_heattr3d'][:]
            # unmask the array, putting zeroes at masked points
            ueiv_heattr[ueiv_heattr.mask] = 0.0
            ueiv_heattr.mask = ma.nomask
        temp_interp[:] = 0.0
        meanU_meanT[:] = 0.0
        Uprime_Tprime[:] = 0.0
        # Calculate fluxes on internal interfaces between T-cells. 
        # Fluxes on external interfaces are zero.
        for ij in range(temp.shape[1]-3):
            temp_interp[:,ij+1,:] = 0.5 * ( temp[:,ij+1,:] + temp[:,ij+2,:] )
            meanU_meanT[:,ij+1,:] = Uvel[:,ij+1,:] * temp_interp[:,ij+1,:]
            Uprime_Tprime[:,ij+1,:] = mean_UT[:,ij+1,:] - meanU_meanT[:,ij+1,:]
        for ij in range(temp.shape[1]-2):
            ttrd_adv_mean[:,ij+1,:] =  ttrd_adv_mean[:,ij+1,:] + \
              ( e3v[:,ij,:]*meanU_meanT[:,ij,:] - e3v[:,ij+1,:]*meanU_meanT[:,ij+1,:] ) \
             / ( e2t[ij+1,:]*e3t[:,ij+1,:] )
            ttrd_adv_eddy[:,ij+1,:] =  ttrd_adv_eddy[:,ij+1,:] + \
              ( e3v[:,ij,:]*Uprime_Tprime[:,ij,:] - e3v[:,ij+1,:]*Uprime_Tprime[:,ij+1,:] ) \
              / ( e2t[ij+1,:]*e3t[:,ij+1,:] )
        if GM:
            for ij in range(temp.shape[1]-2):
                ttrd_adv_GM[:,ij+1,:] =  ttrd_adv_GM[:,ij+1,:] + \
                  ( ueiv_heattr[:,ij,:] - ueiv_heattr[:,ij+1,:] ) \
                  / ( rau0*rcp*e1t[ij+1,:]*e2t[ij+1,:]*e3t[:,ij+1,:] )

        print '>>> calculating z-fluxes'
        # re-use arrays to save memory
        Uvel[:] = 0.0
        mean_UT[:] = 0.0
        Uvel = infile_W.variables['vovecrtz'][:]
        mean_UT = infile_W.variables['wt'][:]
        if GM:
            ueiv_heattr[:] = 0.0
            ueiv_heattr = infile_W.variables['weiv_heattr3d'][:]
            # unmask the array, putting zeroes at masked points
            ueiv_heattr[ueiv_heattr.mask] = 0.0
            ueiv_heattr.mask = ma.nomask
        temp_interp[:] = 0.0
        meanU_meanT[:] = 0.0
        Uprime_Tprime[:] = 0.0
        # Calculate fluxes on internal interfaces between T-cells. 
        # Fluxes on external interfaces are zero.
        for ik in range(temp.shape[0]-3):
            temp_interp[ik+1,:,:] = 0.5 * ( temp[ik,:,:] + temp[ik+1,:,:] )
            meanU_meanT[ik+1,:,:] = Uvel[ik+1,:,:] * temp_interp[ik+1,:,:]
            Uprime_Tprime[ik+1,:,:] = mean_UT[ik+1,:,:] - meanU_meanT[ik+1,:,:]
        for ik in range(temp.shape[0]-1):
            ttrd_adv_mean[ik,:,:] =  ttrd_adv_mean[ik,:,:] + \
                        ( meanU_meanT[ik+1,:,:] - meanU_meanT[ik,:,:] ) / e3t[ik,:,:]
            ttrd_adv_eddy[ik,:,:] =  ttrd_adv_eddy[ik,:,:] + \
                        ( Uprime_Tprime[ik+1,:,:] - Uprime_Tprime[ik,:,:] ) / e3t[ik,:,:]
        if GM:
            for ik in range(temp.shape[0]-2):
                ttrd_adv_GM[ik,:,:] =  ttrd_adv_GM[ik,:,:] + \
                  ( ueiv_heattr[ik+1,:,:] - ueiv_heattr[ik,:,:] ) \
              / ( rau0*rcp*e1t[:,:]*e2t[:,:]*e3t[ik,:,:] )

        # mask the final 3D trends fields to ensure the area integrals are calculated correctly
        # (shouldn't include the zeroes in the outermost rows/columns in the integral).
        ttrd_adv_mean.mask = temp.mask
        ttrd_adv_eddy.mask = temp.mask
        if GM:
            ttrd_adv_GM.mask = temp.mask

        if diagnostics:
            outdata_T.createVariable('ttrd_adv_mean',datatype='f',dimensions=outdata_T.dimensions.keys(),fill_value=-1.0e+20)
            outdata_T.variables['ttrd_adv_mean'][:] = ttrd_adv_mean[:]
            outdata_T.createVariable('ttrd_adv_eddy',datatype='f',dimensions=outdata_T.dimensions.keys(),fill_value=-1.0e+20)
            outdata_T.variables['ttrd_adv_eddy'][:] = ttrd_adv_eddy[:]
            if GM:
                 outdata_T.createVariable('ttrd_adv_GM',datatype='f',dimensions=outdata_T.dimensions.keys(),fill_value=-1.0e+20)
                 outdata_T.variables['ttrd_adv_GM'][:] = ttrd_adv_GM[:]

    ######## Read in (or initialise) required fields ########

    ttrd = {} # will be dictionary of masked arrays
    for fieldname in fields:
        if fieldname == 'ttrd_adv_mean':
            ttrd[fieldname] = ttrd_adv_mean
        elif fieldname == 'ttrd_adv_eddy':
            ttrd[fieldname] = ttrd_adv_eddy
        elif fieldname == 'ttrd_adv_GM':
            ttrd[fieldname] = ttrd_adv_GM
        elif fieldname in ['ttrd_zdfp_minus_evd','ttrd_adv_mean_plus_eddy','ttrd_adv_mean_plus_eddy_plus_GM','ttrd_sum']:
            pass # calculate these later when all other fields read in
        else:   
            ttrd[fieldname] = infile_T.variables[fieldname][:]

    if 'ttrd_zdfp_minus_evd' in fields:
        for fieldname in ['ttrd_zdfp','ttrd_evd']:
            if fieldname not in fields:
                ttrd[fieldname] = infile_T.variables[fieldname][:]
        ttrd['ttrd_zdfp_minus_evd'] = ttrd['ttrd_zdfp'] - ttrd['ttrd_evd']

    if 'ttrd_adv_mean_plus_eddy' in fields:
        ttrd['ttrd_adv_mean_plus_eddy'] = ttrd['ttrd_adv_mean'] + ttrd['ttrd_adv_eddy']
        if diagnostics:
            outdata_T.createVariable('ttrd_adv_mean_plus_eddy',datatype='f',dimensions=outdata_T.dimensions.keys(),fill_value=-1.0e+20)
            outdata_T.variables['ttrd_adv_mean_plus_eddy'][:] = ttrd['ttrd_adv_mean_plus_eddy'][:]

    if 'ttrd_adv_mean_plus_eddy_plus_GM' in fields:
        if GM:
            ttrd['ttrd_adv_mean_plus_eddy_plus_GM'] = ttrd['ttrd_adv_mean'] + ttrd['ttrd_adv_eddy'] + ttrd['ttrd_adv_GM']
            if diagnostics:
                outdata_T.createVariable('ttrd_adv_mean_plus_eddy_plus_GM',datatype='f',dimensions=outdata_T.dimensions.keys(),fill_value=-1.0e+20)
                outdata_T.variables['ttrd_adv_mean_plus_eddy_plus_GM'][:] = ttrd['ttrd_adv_mean_plus_eddy_plus_GM'][:]
        else:
            print "Can't plot ttrd_adv_mean_plus_eddy_plus_GM because there aren't any GM trends in the file."

    # need to do this last of course
    if 'ttrd_sum' in fields:
        fields_to_sum = list(fields) # force it to copy the list not just point to the original
        for field_to_remove in ['ttrd_sum','ttrd_tot']:
            if field_to_remove in fields_to_sum:
                fields_to_sum.remove(field_to_remove)
        ttrd['ttrd_sum'] = ma.copy(ttrd[fields_to_sum[0]])
        for fieldname in fields_to_sum[1:]:
            ttrd['ttrd_sum'][:] = ttrd['ttrd_sum'][:] + ttrd[fieldname][:]

    if totaltrend is not None:
        totalfile = nc.Dataset(totaltrend,'r')
        for varname in ['votemper','tn']:
            try:
                ttrd['model_trend'] = totalfile.variables[varname][:]
            except KeyError:
                pass
            else:
                break
        else:
            raise Exception('Could not find temperature field in model trend file.')
        fields=fields+['model_trend']

    ######## Horizontal average trends and plot ########

    # Re-mask the e3 fields otherwise uninitialised values will contaminate horizontal averages. 
    e3t[:].mask = temp[:].mask

    field_min=None
    field_max=None
    ones3D = ma.ones(temp.shape)
    # apply 3D mask which will be inherited by 3D versions of e1t and e2t.
    ones3D.mask = temp.mask
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
        x_averaged_field = ma.average(ttrd[fieldname][:,:,:]*e1t[:,:]*e2t[:,:]*e3t[:,:,:],axis=-1)
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

    plt.gca().set_xlabel('temperature trend (degC/year)')
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
    parser.add_argument("-f", "--fields", action="store",dest="fields",nargs='+',default=None,
                    help="names of fields to plot")
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

    args = parser.parse_args()

    vert_heat_budget( args.filename, coordsfile=args.coordsfile, fields=args.fields, outfile=args.outfile, 
                      xmin=args.xmin, xmax=args.xmax, depthmax=args.depthmax, 
                      totaltrend=args.totaltrend, integralcheck=args.integralcheck, diagnostics=args.diagnostics )
