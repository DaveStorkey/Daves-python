#! /usr/bin/env python

'''
Routine to calculate and plot global mean drifts of T & S (or any other 3D field on the T-grid)
versus depth.

June 2022 : Upgrade to python 3. DS.

@author: Dave Storkey
'''

import socket
import matplotlib
if 'spice' in socket.gethostname():
    # Note this disables plt.show()
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as col
import netCDF4 as nc
import numpy as np
import numpy.ma as ma

def global_drifts(filenames=None,coordsname=None,tzfilename=None,field=None,thickname=None,
                  title=None,nlevs=None,mnfld=None,mxfld=None,xmin=None,xmax=None,ymin=None,ymax=None,
                  depthmax=None,outfile=None,time_offset=None,absolute_field=None,
                  west=None,east=None,south=None,north=None,maskfilename=None,maskfieldname=None,
                  noshow=False,nobar=False,calc_only=False):

###########################################
# 1. Read in the data and create t-z field
###########################################

    if tzfilename is None:
        coords = nc.Dataset(coordsname,'r')
        # meshmask and coordinates files are generally not masked. Force e1t and e2t
        # to be masked arrays here and set the mask from the input field later.
        e1t = ma.array(coords.variables['e1t'][:])
        e2t = ma.array(coords.variables['e2t'][:])
        selected_area = None
        mask_tag = None
        if west is not None or east is not None or south is not None or north is not None:
            mask_tag=''
            # read in lat/lon arrays and cope with east-west wrapping
            lon = coords.variables['glamt'][:]
            # restrict longitudes to the range [0,360] 
            lon = np.remainder(lon,360.0)
            lat = coords.variables['gphit'][:]
            if west is not None:
                west = west%360.0
            if east is not None:
                east = east%360.0
            # if we have chosen west and east limits such that the area crosses the
            # zero meridion we'll have to change to use [-180,180]
            if west is not None and east is not None and west > east:
                select_mask = np.where(lon > 180.0,1,0)
                lon = lon - 360.0*select_mask
                if west > 180.0:
                    west=west-360.0
                if east > 180.0:
                    east=east-360.0

            # create mask for specified meaning area
            # note that the mask has *True* at *masked* points.
            selected_area = ma.zeros(lon.shape)
            if west is not None:
                selected_area.mask = selected_area.mask | ( lon < west )
                mask_tag = mask_tag + '_'+'W'+str(west)
            if east is not None:
                selected_area.mask = selected_area.mask | ( lon > east )
                mask_tag = mask_tag + '_'+'E'+str(east)
            if south is not None:
                selected_area.mask = selected_area.mask | ( lat < south )
                mask_tag = mask_tag + '_'+'S'+str(south)
            if north is not None:
                selected_area.mask = selected_area.mask | ( lat > north )
                mask_tag = mask_tag + '_'+'N'+str(north)

        if maskfilename is not None:
            maskfile = nc.Dataset(maskfilename,'r')
            if maskfieldname is None:
                maskfieldname = 'mask'
            if mask_tag is None:
                mask_tag = '_'+'mask_'+maskfieldname
            else:
                mask_tag = mask_tag + '_'+'mask_'+maskfieldname
            maskfield = maskfile.variables[maskfieldname][:]
            if selected_area is None:
                selected_area = ma.zeros(e1t.shape)
                selected_area.mask = ( maskfield == 0 )
            else:
                selected_area.mask = selected_area.mask | ( maskfield == 0 ) | maskfield.mask

        tzfield=[]
        for filename in filenames:
            print("Working on file : "+filename)
            infile = nc.Dataset(filename,'r')
            infield = infile.variables[field]
            # Counting masked points is a more secure check than just checking for the existence
            # of the _FillValue attribute. Can have the attribute with no masked points. 
            # Only do the check for the first file because it is slow.
            if filename == filenames[0]:
                if ma.count_masked(infield[:]) == 0:
                    raise Exception("Error: field "+field+" in file "+filename+
                                    " is not masked, so global spatial average will be incorrect.")
                else: # create level-by-level average cell area
                    if len(infield.shape) == 4:
                        infield0 = infield[0]
                    else:
                        infield0 = infield
                    ones3D = ma.ones(infield0.shape)
                    # apply 3D mask which will be inherited by 3D versions of e1t and e2t.
                    ones3D.mask = infield0.mask
                    # broadcast 2D scale factor fields to 3D
                    e1t_3d = e1t[:] * ones3D[:]
                    e2t_3d = e2t[:] * ones3D[:]
                    average_cell_area = ma.average(ma.average(e1t_3d[:]*e2t_3d[:],axis=-1),axis=-1)

            # mask out required area if specified
            if selected_area is not None:
                # create a copy of the field and mask with the selected area mask
                infield_data = infield[:].copy()
                # note that if infield has a time dimension the bitwise "or" operator will broadcast 
                # the select_area mask along the time dimension
                infield_data.mask = infield_data.mask | selected_area.mask
                # write out the first masked field as a check that it is doing what you expect
                if filename == filenames[0]:
                    maskfile_out = nc.Dataset(filename.split("/")[-1].replace('.nc',mask_tag+'.nc'),'w')
                    for dim in infield.dimensions:
                        maskfile_out.createDimension(dim,len(infile.dimensions[dim]))
                    maskfile_out.createVariable(field,datatype="f",dimensions=infield.dimensions,fill_value=-1.0e+20)
                    maskfile_out.variables[field][:] = infield_data[:]
                    maskfile_out.close()
            else:
                # create a simple pointer to the field
                infield_data = infield[:] 

            # loop over the time dimension if there is one. Assume it is the first dimension.
            if infile.dimensions[infield.dimensions[0]].isunlimited():
                for rec in range(len(infile.dimensions[infield.dimensions[0]])):
                    x_averaged_field = ma.average(infield_data[rec,:,:,:]*e1t[:,:]*e2t[:,:],axis=-1)
                    tzfield.append(ma.average(x_averaged_field,axis=-1)/average_cell_area[:])
            else:
                x_averaged_field = ma.average(infield_data[:,:,:],axis=-1)
                tzfield.append(ma.average(x_averaged_field,axis=-1)/average_cell_area[:])

        # Convert to a drift from initial conditions:
        if not absolute_field:
            tzfield0 = np.copy(tzfield[0])
            for rec in range(len(tzfield)):
                tzfield[rec][:] = tzfield[rec][:] - tzfield0[:]

        ntimes = len(tzfield)
        ndepths = len(tzfield[0])
        tzfield = ma.concatenate(tzfield)
        tzfield = tzfield.reshape(ntimes,ndepths)
        tzfield = tzfield.transpose()    
        print('shape of tzfield : ',tzfield.shape)

        # Label the time axis properly eventually...
        time=np.array(range(ntimes))
        for depthname in ['depth', 'deptht', 'depthu', 'depthv', 'depthw']:
            try:
                depth = infile.variables[depthname]
            except KeyError:
                pass
            else:
                break
        else:
            raise Exception('Could not find depth coordinate variable in file.')

        # write out field for future reference
        keepfilename = filenames[0][0:-3].split("/")[-1]+"_"+filenames[-1][0:-3].split("/")[-1]+"_tz_"+field
        if mask_tag is not None:
            keepfilename = keepfilename + mask_tag
        keepfilename = keepfilename+".nc"
        keepfile = nc.Dataset(keepfilename, 'w')
        # Putting None for the length of the time dimension makes it unlimited.
        # BUT the unlimited dimension is the second dimension in tzfield
        # whereas it normally should be the first dimension. This confuses NCO 
        # operators - should fix it. 
        for (dim,dimlen) in zip(["time","depth"],[None,ndepths]):
            keepfile.createDimension(dim, dimlen)
        if selected_area is not None:
            (ny,nx) = selected_area.shape
            for (dim,dimlen) in zip(["x","y"],[nx,ny]):
                keepfile.createDimension(dim, dimlen)
        keepfile.createVariable("time",datatype="f",dimensions=("time"))
        keepfile.variables["time"][:] = time[:]
        keepfile.createVariable("depth",datatype="f",dimensions=("depth"))
        keepfile.variables["depth"][:] = depth[:]
        keepfile.createVariable(field,datatype="f",dimensions=("depth","time"),fill_value=-1.0e+20)
        keepfile.variables[field][:,:] = tzfield[:,:]
        keepfile.close()

########################################################
# 1a. ...or just read in previously-calculated t-z field
########################################################

    else:
        keepfile = nc.Dataset(tzfilename, 'r')
        time = ma.copy(keepfile.variables["time"][:])
        depth = ma.copy(keepfile.variables["depth"][:])
        tzfield = ma.copy(keepfile.variables[field][:,:])
        keepfile.close()

############################################################
# 2. Do the time-depth plot if no thickness file specified
############################################################

    if time_offset is not None:
        time[:] = time[:] + time_offset

    if xmin is None:
        xmin = time.min()
    if xmax is None:
        xmax = time.max()

    if thickname is None and not calc_only: 

        if mnfld is None:
            mnfld = tzfield.min()
        if mxfld is None:
            mxfld = tzfield.max()

        print("mnfld,mxfld: ",mnfld,mxfld)
        levs = np.linspace(mnfld,mxfld,nlevs+1) 
        print("levs: ",levs)
        cmap = getattr(matplotlib.cm,'RdBu_r')
#        ind = np.intp(np.rint(np.linspace(2,255,nlevs+1))).tolist()    
#    # ensure we have white space in the middle for an odd number of levels.
#        if 2*(nlevs/2) != nlevs:
#            ind[nlevs/2] = 127
#        print('ind',ind)
#        colorlist=cmap(ind)

        plt.gca().patch.set_facecolor('grey')
    #    np.set_printoptions(threshold='nan')

    #    clev=plt.contourf(time,depth,tzfield,colors=colorlist,levels=levs)
        clev=plt.contourf(time,depth,tzfield,cmap=cmap,levels=levs)
        if title is not None:
            plt.gca().set_title(title)
        plt.gca().set_xlim([xmin,xmax])
        if depthmax is not None:
            plt.gca().set_ylim([depthmax,0.0])
        else:
            plt.gca().invert_yaxis()
 
        # Label the time axis properly eventually...
        plt.gca().set_xlabel('year of spin up')
        plt.gca().set_ylabel('depth (m)')
        plt.gca().ticklabel_format(style='plain',axis='x',useOffset=False)

        if not nobar:
            cax = plt.colorbar(clev,orientation='horizontal')
            cax.set_ticks(levs)
            labels=["%1.2f" %lev for lev in levs]
            cax.ax.set_xticklabels(labels,rotation=45)
    
        if outfile is not None:
            matplotlib.rcParams['font.size'] =8
            plt.savefig(outfile,dpi=200)
        else:
            if not noshow:
                plt.show()

#####################################################################################
# 3. Do the time-series plot of volume-integrated drifts if thickness file specified
#####################################################################################

    if thickname is not None and not calc_only: 
        thickfile = nc.Dataset(thickname,'r')
        for var in ['e3t','e3t_0']:
            try:
                e3t = ma.array(thickfile.variables[var][:])
            except KeyError:
                pass
            else:
                break
        else:
            raise Exception('Could not find thickness variable in file.')
        e3t.mask = tzfield[:,:].transpose().mask
        tseries_field = ma.average(tzfield[:,:].transpose()*e3t[:],-1)/ma.average(e3t[:])
        print('tseries_field.shape : ',tseries_field.shape)
        if ymin is None:
            ymin = tseries_field.min()
        if ymax is None:
            ymax = tseries_field.max()
        plt.plot(time,tseries_field)
        plt.gca().set_xlim([xmin,xmax])
        plt.gca().set_ylim([ymin,ymax])
        plt.gca().ticklabel_format(style='plain',axis='x',useOffset=False)
        plt.gca().set_xlabel('year of spin up')
        if field == 'votemper' or field == 'temperature' or field == 'thetao':
            plt.gca().set_ylabel('temperature drift from initial condition (deg C)')
        elif field == 'vosaline' or field == 'salinity' or field == 'so':
            plt.gca().set_ylabel('salinity drift from initial condition (psu)')
        else:
            plt.gca().set_ylabel('field drift from initial condition (unknown units)')
        if outfile is not None:
            matplotlib.rcParams['font.size'] =8
            plt.savefig(outfile,dpi=200)
        else:
            if not noshow:
                plt.show()
    

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--filenames", action="store", dest="filenames", help="name of input file(s)",nargs="+")
    parser.add_argument("-c", "--coordsname", action="store", dest="coordsname",
                    help="name of coordinates file")
    parser.add_argument("-d", "--thickname", action="store",dest="thickname",default=None,
                    help="name of file with cell thicknesses: if specified will plot volume averaged timeseries")
    parser.add_argument("-I", "--tzfilename", action="store", dest="tzfilename", 
                    help="name of input file containing previously-calculated t-z field")
    parser.add_argument("-j", "--field", action="store", dest="field", help="name of field to plot")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",default=None,
                    help="plot to output file - filetype by extension. If unset plot to GUI.")
    parser.add_argument("-n", "--nlevs", action="store",dest="nlevs",type=int,default=15,
                    help="number of contour levels")
    parser.add_argument("-x", "--xmin", action="store",dest="xmin",type=float,default=None,
                    help="start of x-axis (record number)")
    parser.add_argument("-X", "--xmax", action="store",dest="xmax",type=float,default=None,
                    help="end of x-axis (record number)")
    parser.add_argument("-y", "--ymin", action="store",dest="ymin",type=float,default=None,
                    help="start of y-axis if timeseries plot")
    parser.add_argument("-Y", "--ymax", action="store",dest="ymax",type=float,default=None,
                    help="end of y-axis if timeseries plot")
    parser.add_argument("-f", "--mnfld", action="store",dest="mnfld",type=float,default=None,
                    help="minimum field value to plot")
    parser.add_argument("-F", "--mxfld", action="store",dest="mxfld",type=float,default=None,
                    help="maximum field value to plot")
    parser.add_argument("-D", "--depthmax", action="store",dest="depthmax",type=float,default=None,
                    help="maximum depth for plot")
    parser.add_argument("-T", "--time_offset", action="store",dest="time_offset",type=float,default=None,
                    help="offset to add to time axis")
    parser.add_argument("-t", "--title", action="store",dest="title",default=None,
                    help="title for plot")
    parser.add_argument("-A", "--absolute", action="store_true",dest="absolute_field",
                    help="plot absolute field instead of drift from initial field")
    parser.add_argument("-W", "--west", action="store",dest="west",type=float,default=None,
                    help="western limit of area to average")
    parser.add_argument("-E", "--east", action="store",dest="east",type=float,default=None,
                    help="eastern limit of area to average")
    parser.add_argument("-S", "--south", action="store",dest="south",type=float,default=None,
                    help="southern limit of area to average")
    parser.add_argument("-N", "--north", action="store",dest="north",type=float,default=None,
                    help="northern limit of area to average")
    parser.add_argument("-m", "--maskfile", action="store",dest="maskfilename",default=None,
                    help="name of file containing 1/0 mask")
    parser.add_argument("-M", "--maskfield", action="store",dest="maskfieldname",default=None,
                    help="name of mask field (default = 'mask')")
    parser.add_argument("-b", "--nobar", action="store_true",dest="nobar",
                    help="don't plot a colour bar")
    parser.add_argument("-s", "--noshow", action="store_true",dest="noshow",
                    help="don't run plt.show() even if output file unset")
    parser.add_argument("-C", "--calc_only", action="store_true",dest="calc_only",
                    help="Create tz file and don't plot anything.")

    args = parser.parse_args()

    global_drifts( filenames=args.filenames, coordsname=args.coordsname, tzfilename=args.tzfilename, 
    field=args.field, outfile=args.outfile, title=args.title, nlevs=args.nlevs,
    mnfld=args.mnfld, mxfld=args.mxfld,xmin=args.xmin,xmax=args.xmax,ymin=args.ymin,ymax=args.ymax, 
    depthmax=args.depthmax,time_offset=args.time_offset, thickname=args.thickname,absolute_field=args.absolute_field,
    west=args.west, east=args.east, south=args.south, north=args.north, maskfilename=args.maskfilename, 
    maskfieldname=args.maskfieldname, noshow=args.noshow, nobar=args.nobar, calc_only=args.calc_only )

