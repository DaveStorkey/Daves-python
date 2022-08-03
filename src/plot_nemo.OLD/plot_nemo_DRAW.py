#! /usr/bin/env python
'''
Wrapper for Tim Graham's nemo_iris_plot.contourf routine (originally). 

Created November 2013

August 2018 : This version completely rewritten:
              1. Does vector plotting. 
              2. Doesn't use nemo_iris_plot any more. 
              3. Everything done with cartopy now - no more Basemap.
              4. Cut out areas work with any projection.
              5. Can overplot a line contour or vector arrows on a colour-filled plot. 
DS.

Oct 2018 : Sorted out a more flexible way of specifying colours and colour maps. DS.

Oct 2018 : Added a key arrow for vector plots and labelled the colour bar with units. DS.

Nov 2018 : Added a "draw" option to draw line segments between points. DS.

To do:        1. Copy the block-plotting option (pcolormesh) from plot_nemo_section.py.


@author: Dave Storkey
'''

import argparse
import socket
import matplotlib
if 'spice' in socket.gethostname():
    # Note this disables plt.show()
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.figure as fig
import matplotlib.colors as mcolors
import matplotlib.cm
import iris
import cartopy.crs as ccrs
import numpy as np
import numpy.ma as ma
import monty
import textwrap
import U_at_T_points as UaT
import general_tools as gt

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def project_cube(cube_in, proj=None):

    cube_out, extent = iris.analysis.cartography.project(cube_in, proj, nx=cube_in.shape[-1], ny=cube_in.shape[-2])
    x = cube_out.coord('projection_x_coordinate').points
    y = cube_out.coord('projection_y_coordinate').points
    return x,y,cube_out

def get_coord(cube, coordname, filename):

#   Lats/lons will be recognised by Iris as auxillary coordinates if the "coordinates"
#   attribute of the field in the file has been set correctly. If not, the lat/lon fields
#   might still be in the file, so try reading them as simple fields rather than coordinates.

    try:
        coord = cube.coord(coordname).points
    except iris.exceptions.CoordinateNotFoundError:
        print 'reading '+coordname+' as field'
        coord=iris.load_cube(filename,coordname).data
    else:
        print 'reading '+coordname+' as auxillary coordinate'

    return coord

def plot_nemo(filenames=None,sca_names=None,vec_names=None,nlevs=13,mnfld=None,mxfld=None,title=None,rec=0,level=1,bottom=False,scientific=False,
              cmap=None,colors=None,reverse_colors=False,west=None,east=None,south=None,north=None,proj='pc0',maskfile=None,outfile=None,
              logscale=None,factor=None,lineplot=False,zeromean=False,vec_only=False,mag_only=False,arrows=None,whitebg=False,
              draw_points=None,draw_fmt=None):

    # short cuts
    projections = { 'pc0'      : ( ccrs.PlateCarree() , 'global' ),
                    'pc180'    : ( ccrs.PlateCarree(central_longitude=180.0), 'global' ), 
                    'merc'     : ( ccrs.Mercator(min_latitude=-86.0), 'global' ), 
                    'northps'  : ( ccrs.NorthPolarStereo(), [0,359.9,50,90] ), 
                    'southps'  : ( ccrs.SouthPolarStereo(), [0,359.9,-90,-50] ) 
    }   

    p = projections[proj]
        
    nscalar = 0
    if lineplot:
        splots=['l','c']
        sc = 1
    else:
        splots=['c','l']
        sc = 0
    nvector = 0
    regrid_vector = False
    filenames_to_read = []
    varnames_to_read = []
    if filenames is None or len(filenames) > 4:
        raise Exception('Error: specify between one and four filenames.')
    if sca_names is None and vec_names is None:
        raise Exception('Error: must specify at least one scalar or vector field.')
    if sca_names is not None:
        if len(sca_names) > 2: 
            raise Exception('Error: specify zero, one or two scalar fields')
        else:
            nscalar = len(sca_names)
            if len(filenames) < len(sca_names):
                filenames = [filenames[0],filenames[0]]
            for filename in filenames[0:len(sca_names)]:
                filenames_to_read.append(filename)
                filenames.remove(filename)
            varnames_to_read = varnames_to_read + sca_names
    if vec_names is not None:
        if len(vec_names) != 2: 
            raise Exception('Error: specify zero or two vector components')
        else:
            nvector = 1
            v1 = nscalar
            v2 = nscalar+1
            # If separate files specified for each vector component assume we need to regrid to the T-grid.
            # Otherwise both vector components in the same file and we don't need to regrid.
            if len(filenames) == 2:
                regrid_vector = True
            elif len(filenames) == 1:
                filenames = [filenames[0],filenames[0]] 
            else:
                raise Exception('Error: need to specify one or two filenames for the vector components.')
            filenames_to_read = filenames_to_read + filenames
            varnames_to_read = varnames_to_read + vec_names
 
    fld_in = []
    for filename,varname in zip(filenames_to_read,varnames_to_read):
        print 'Reading '+varname+' from '+filename
        fld_in.append(gt.read_cube(filename,varname))
        
    print 'number of scalar fields to plot : ',nscalar
    print 'number of vector fields to plot : ',nvector
    print 'len(fld_in) :', len(fld_in)
    print 'regrid_vector is : ',regrid_vector

# Organise control parameters that apply to individual fields.
# In general with multiple fields and a single value of the control
# parameter then that value gets applied to all fields, otherwise
# the number of values input has to match the number of fields, where
# the vector field counts as a single field. 

    nfields = nscalar + nvector
    if type(rec) is not list:
        rec = (nscalar + 2*nvector) * [rec]
    elif len(rec) != nfields:
        raise Exception("Error: number of values for REC (-r) doesn't match number of fields.")
    else:
        rec = rec + nvector*[rec[-1]]
    if type(level) is not list:
        level = (nscalar + 2*nvector) * [level]
    elif len(level) != nfields:
        raise Exception("Error: number of values for LEVEL (-l) doesn't match number of fields.")
    else:
        level = level + nvector*[level[-1]]
    if type(maskfile) is not list:
        maskfile = (nscalar + 2*nvector) * [maskfile]
    elif len(maskfile) != nfields:
        raise Exception("Error: number of values for MASKFILE (-m) doesn't match number of fields.")
    else:
        maskfile = maskfile + nvector*[maskfile[-1]]
    if type(mnfld) is not list:
        mnfld = nfields*[mnfld]
    elif len(mnfld) != nfields:
        raise Exception("Error: number of values for MNFLD (-f) doesn't match number of fields.")
    if type(mxfld) is not list:
        mxfld = nfields*[mxfld]
    elif len(mxfld) != nfields:
        raise Exception("Error: number of values for MXFLD (-F) doesn't match number of fields.")
    if type(bottom) is not list:
        bottom = (nscalar + 2*nvector) * [bottom]
    elif len(bottom) != nfields:
        raise Exception("Error: number of values for BOTTOM (-b) doesn't match number of fields.")
    else:
        bottom = bottom + nvector*[bottom[-1]]
    if type(nlevs) is not list:
        nlevs = nfields*[nlevs]
    elif len(nlevs) != nfields:
        raise Exception("Error: number of values for NLEVS (-n) doesn't match number of fields.")
    if type(factor) is not list:
        factor = (nscalar + 2*nvector) * [factor]
    elif len(factor) != nfields:
        raise Exception("Error: number of values for FACTOR (-x) doesn't match number of fields.")
    else:
        factor = factor + nvector*[factor[-1]]

    if nscalar > 1:
        if type(zeromean) is not list:
            zeromean = nscalar*[zeromean]
        elif len(zeromean) != nscalar:
            raise Exception("Error: number of values for ZEROMEAN (-x) doesn't match number of scalar fields.")

# NB. Assumes time is the first dimension.
    fld = []
    for fldi,rec_i in zip(fld_in,rec):
        if fldi.ndim == 4:
            fld.append(fldi[rec_i])
        else:
            fld.append(fldi)

    if maskfile is not None:
        for fldi,maskfile_i in zip(fld,maskfile):
# mask the field using the mask field in maskfile
            if maskfile_i is not None:
                for maskname in ['mask', 'bathy', 'bathymetry', 'Bathymetry']:
                    try:
                        mask = iris.load_cube(maskfile_i,iris.Constraint(cube_func=lambda nemoname: nemoname.var_name == maskname) ) 
                    except iris.exceptions.ConstraintMismatchError:
                        pass
                    else:
                        break
                else:
                    raise iris.exceptions.ConstraintMismatchError

                try:
                    maskvalue = fld._FillValue
                except AttributeError:
                    maskvalue = 0.0

# assume the mask field has the correct number of dimensions!
                fldi.data = ma.masked_where( mask.data == maskvalue, fldi.data )

    # interpolate to T-points. If we are plotting the vector field at the bottom then 
    # do this *before* selecting the bottom field:
    if regrid_vector and bottom[v1]:
        fld[v1],fld[v2] = UaT.U_at_T_points(fld[v1],fld[v2])

# extract bottom field if required using mask...
    fldslice=[]
    for fldi,bottom_i,level_i in zip(fld,bottom,level):
        if bottom_i and fldi.ndim == 3:
            fldslice.append(fldi[0].copy())
#            fldslice[-1][:] = 0.0
            if fldi.var_name == 'vovecrtz' or fldi.var_name == 'wo':
                # for vertical velocities take the level just above the bottom (because the bottom field is zero by construction).
                for lev in range(len(fldi.data[:,0,0])-1):
                    fldslice[-1].data = np.ma.where( (~fldi.data[lev,:,:].mask & fldi.data[lev+1,:,:].mask), 
                                             x=fldi.data[lev-1,:,:], y=fldslice.data )
            else:
                print 'fldi[0].data.mask.shape : ',fldi[0].data.mask.shape
                for lev in range(len(fldi.data[:,0,0])-1):
                    fldslice[-1].data[:,:] = np.ma.where( (~fldi.data.mask[lev,:,:] & fldi.data.mask[lev+1,:,:]), 
                                                  x=fldi.data[lev-1,:,:], y=fldslice[-1].data[:,:] )

# ... or select required level
        else:
            if fldi.ndim == 3:
               fldslice.append(fldi[level_i-1])
            else:
               fldslice.append(fldi)


        print 'min/max fldslice : ',fldslice[-1].data.min(),fldslice[-1].data.max()

    # Do interpolation and rotation of vector fields *after* we have 
    # selected a single record and single level - much more efficient!

    # interpolate to T-points:
    if regrid_vector and not bottom[v1]:
        fldslice[v1],fldslice[v2] = UaT.U_at_T_points(fldslice[v1],fldslice[v2])

    # rotate to be relative to geographical coordinates rather than model grid:
    # (unnecessary if we are only plotting the magnitude of the vector field).
    if nvector and not mag_only:
        fldslice[v1],fldslice[v2] = UaT.rotate(fldslice[v1],fldslice[v2])

    # read in lat/lon coords for each scalar field slice (including vector magnitude)
    lons=[]
    lats=[]
    for filenames_to_read_i,fldslice_i in zip(filenames_to_read,fldslice):
        lons.append(get_coord(fldslice_i,'longitude',filenames_to_read_i))
        lats.append(get_coord(fldslice_i,'latitude',filenames_to_read_i))

    # restrict longitudes to the range [0,360] 
    for lons_i in lons:
        lons_i[:] = np.remainder(lons_i[:],360.0)
    if west is not None:
        west = west%360.0
    if east is not None:
        east = east%360.0

    # if we have chosen west and east limits such that the area crosses the
    # zero meridion we'll have to change to use [-180,180]
    if west is not None and east is not None and west > east:
        for lons_i in lons:
            select_mask = np.where(lons_i[:] > 180.0,1,0)
            lons_i[:] = lons_i[:] - 360.0*select_mask
        if west > 180.0:
            west=west-360.0
        if east > 180.0:
            east=east-360.0
        print 'min/max lons_i after normalisation: ',lons_i.min(),lons_i.max() 

    if p[1] == 'global':
        (west_default,east_default,south_default,north_default) = (0.0,359.99,-90.0,90.0)
    else:
        (west_default,east_default,south_default,north_default) = p[1]

    if west is None and east is None and south is None and north is None:
        region = p[1]
    else:
        if west is None:
            west = max(west_default,min(lons_i.min() for lons_i in lons))
        if east is None:
            east = min(east_default,max(lons_i.max() for lons_i in lons))
        if south is None:
            south = max(south_default,min(lats_i.min() for lats_i in lats))
        if north is None:
            north = min(north_default,max(lats_i.max() for lats_i in lats))
        region = [west,east,south,north]

    print 'region : ',region

    # Mask out (or in the future cut out) bits of the field(s) to be contoured outwith the area being plotted 
    # so that the automatic min/max calculation works. Also so that the zero mean calculation is restricted to
    # the plotted area. Also so we can calculate the maximum vector length for the key arrow.

    if region != 'global':
        for lons_i, lats_i, fldslice_i in zip(lons[:],lats[:],fldslice[:]): 
            fldslice_i.data = ma.masked_where( ((lons_i < region[0]) | (lons_i > region[1])) | ((lats_i < region[2]) | (lats_i > region[3])) , fldslice_i.data )

    # Apply factor if required:

    if factor is not None:
        for fldslice_i,factor_i in zip(fldslice,factor):
            if factor_i is not None:
                fldslice_i.data = factor_i * fldslice_i.data

    # Zero out the spatial mean over the displayed area if required:

    if zeromean:
        if not nscalar:
            raise Exception("Error : zero mean option doesn't work for vector fields.")
        else:
            for filename_to_read,fldslice_i,zeromean_i in zip(filenames_to_read[:nscalar],fldslice[:nscalar],zeromean):
                print 'fldslice_i.standard_name, zeromean_i : ',fldslice_i.standard_name, zeromean_i
                if zeromean_i:
                    cell_area = iris.load_cube(filename_to_read,'cell_area')
                    area_mean = fldslice_i.collapsed(['longitude','latitude'],iris.analysis.MEAN,weights=cell_area.data)
                    fldslice_i.data = fldslice_i.data - area_mean.data

    # Calculate magnitude of vector field if there is one. If we are only plotting the vector field as arrows then 
    # this is used later to calculate the length of the key arrow. If we are also contouring 
    # the magnitude of the vector field (only done if no scalar fields to plot) then add the magnitude
    # field slice to the beginning of the fldslice list.
     
    if nvector:
        u = fldslice[v1].data
        v = fldslice[v2].data
        speed = ma.sqrt(u*u + v*v)
        if not nscalar and not vec_only:
            # Append vector magnitude to the beginning of the fldslice list.
            nscalar=1; v1=1; v2=2; sc=0;
            fldslice = [fldslice[0].copy()] + fldslice
            # hardwired for now - doesn't really matter
            fldslice[0].standard_name = "sea_water_speed"
            fldslice[0].data = speed[:]

    # Calculate contour levels for contour plots and colour map for colour-filled contouring.

    if nscalar or not vec_only:
        levs = []
        for mnfld_i,mxfld_i,nlevs_i,fldslice_i,splot_i in zip(mnfld[:nscalar],mxfld[:nscalar],nlevs[:nscalar],fldslice[:nscalar],splots[:nscalar]):
            if mnfld_i is None:
                mnfld_i = ma.min(fldslice_i.data)
            if mxfld_i is None:
                mxfld_i = ma.max(fldslice_i.data)

            print 'nlevs_i : ',nlevs_i
            if logscale and splot_i == 'c':
                tiny_factor=0.0001
                if mnfld_i == 0.0:
                    mnfld_i = tiny_factor * mxfld_i
                if mxfld_i == 0.0:
                    mxfld_i = tiny_factor * mnfld_i
                if mnfld_i < 0.0 and mxfld_i < 0.0:
                    levs_array = -np.logspace(np.log10(-mnfld_i),np.log10(-mxfld_i),nlevs_i+1) 
                if mnfld_i > 0.0 and mxfld_i > 0.0:
                    levs_array = np.logspace(np.log10(mnfld_i),np.log10(mxfld_i),nlevs_i+1) 
                if mnfld_i < 0.0 and mxfld_i > 0.0:
                    nlevs_neg = int(abs(mnfld_i)*nlevs_i/(mxfld_i-mnfld_i))
                    levs_neg = -np.logspace(np.log10(-mnfld_i),np.log10(-mnfld_i*tiny_factor),nlevs_i+1) 
                    nlevs_pos = int(mxfld_i*nlevs_i/(mxfld_i-mnfld_i))
                    levs_pos = np.logspace(np.log10(mxfld_i*tiny_factor),np.log10(mxfld_i),nlevs_i+1) 
                    levs_array = levs_neg + levs_pos
            else:
                levs_array = np.linspace(mnfld_i,mxfld_i,nlevs_i+1) 
            levs.append(levs_array)

            if colors is None and nscalar == 2:
                # for a line contour plot over a filled contour plot default to black lines
                colors = "black"

            if colors is not None and nscalar == 1:
                # colors takes precedence over cmap and matplotlib.pyplot.contour complains
                # if you give it colors and cmap both as not None.
                cmap = None

            # sort out colour map for filled or line contouring
            if cmap is not None:
                if cmap == 'Rainbow':
                    # shortcut for old IDL rainbow colour map:
                    cmap = monty.clr_cmap('/home/h04/frsy/IDL/rainbow_diff_nice.clr')
                elif '/' in cmap:
                    # this is pointing to a file so use monty.clr_cmap:                
                    print 'reading colour map using monty.clr_cmap'
                    cmap = monty.clr_cmap(cmap)
                else:
                    # assume that it refers to a matplotlib inbuilt colour map:
                    print 'reading inbuilt matplotlib colour map'
                    cmap = getattr(matplotlib.cm,cmap)
                if reverse_colors:
                    cmap = mcolors.ListedColormap(cmap(np.arange(255,-1,-1)))

#                ind = np.intp(np.rint(np.linspace(2,255,nlevs_i+1))).tolist()    
# ensure we have white space in the middle for an odd number of levels.
#                if 2*(nlevs_i/2) != nlevs_i:
#                    ind[nlevs_i/2] = 127
#                print 'ind',ind
#                cmap=cmap(ind)

        print 'len(levs) : ',len(levs)

    # Create new "projected" cubes. Shouldn't really have to do this but need
    # it a workaround for a cartopy bug. Note that the fields are all on the
    # same grid, so (x,y) is the same for each field.

    fldproj = []
    x=[]
    y=[]
    for fldslice_i in fldslice:
        projinfo = project_cube(fldslice_i, proj=p[0])
        x.append(projinfo[0])
        y.append(projinfo[1])
        fldproj.append(projinfo[2])

    plt.figure(1)

    # Create plot axes:

    ax = plt.axes(projection=p[0])
    ax.coastlines()
    if not whitebg:
        ax.stock_img()
    # set the plot extent
    if region == 'global':
        ax.set_global()
    else:
        ax.set_extent(region, crs=ccrs.PlateCarree())
        # ax.set_extent(region, crs=ccrs.Geodetic())

    # Contour plot(s) for the scalar field(s), including magnitude of vector field if plotting.

    if nscalar:
        for x_i,y_i,fldproj_i,splot_i,levs_i,fld_in_i in zip(x[:nscalar],y[:nscalar],fldproj[:nscalar],splots[:nscalar],levs[:nscalar],fld_in[:nscalar]):
            print 'levs_i : ',levs_i
            if splot_i == 'l':
                cs1=plt.contour(x_i,y_i,fldproj_i.data,levels=levs_i,mxfld=levs_i[-1],mnfld=levs_i[0],colorbar=None,colors=colors,cmap=cmap,linewidths=0.5)
            elif splot_i == 'c':
                # always use colour maps for the filled contouring
                cs1=plt.contourf(x_i,y_i,fldproj_i.data,levels=levs_i,mxfld=levs_i[-1],mnfld=levs_i[0],colorbar=None,cmap=cmap,transform=p[0])
                cax=plt.colorbar(orientation='horizontal')
                cax.set_ticks(levs_i)
                if scientific:
                    labels=["%1.2e" %lev for lev in levs_i]
                else:
                    labels=["%1.3f" %lev for lev in levs_i]
                cax.ax.set_xticklabels(labels,rotation=45)
                cax.ax.set_xlabel(str(fld_in_i.units))
            else:
                raise Exception("Error: Something has gone wrong. Splot should be 'l' or 'c'")

    # Arrows for vector fields

    if nvector and not mag_only:
        arrow_subsample = 10
        arrow_scale = 15
        if not nscalar and vec_only:
            arrow_colour = "black"
        else:
            arrow_colour = "white"
        key_arrow_length = 0.9 * ma.amax(speed[::arrow_subsample,::arrow_subsample])
        if arrows is not None:
            arrow_subsample = int(arrows[0])
            if len(arrows) > 1:
                arrow_scale = int(arrows[1])
            if len(arrows) > 2:
                arrow_colour=arrows[2]
            if len(arrows) > 3:
                key_arrow_length=arrows[3]
        cs2 = plt.quiver(x[v1][::arrow_subsample], y[v1][::arrow_subsample], 
                        fldproj[v1].data[::arrow_subsample,::arrow_subsample], fldproj[v2].data[::arrow_subsample,::arrow_subsample], 
                        color=arrow_colour, scale=arrow_scale, angles='xy' )
        plt.quiverkey(cs2, X=0.8, Y=-0.05, U=key_arrow_length, label='%1.1f'%key_arrow_length+str(fld_in[v1].units), labelpos='E')

    # Line segments

    if draw_points is not None:
        if len(draw_points)%4 != 0:
            raise Exception("Error: number of draw_points (-d) must be a multiple of 4: start_lon,start_lat,end_lon,end_lat")
        fmt = "k-"  # default to black solid lines        
        linewidth = 2
        if draw_fmt is not None:
            fmt = draw_fmt[0]
            if len(draw_fmt) == 2:
                linewidth = draw_fmt[1] 
        for ii in range(0,len(draw_points),4):
            # pyplot.plot takes all the lons first and the lats second...
            plot_lons = [draw_points[ii],draw_points[ii+2]]
            plot_lats = [draw_points[ii+1],draw_points[ii+3]]
            print 'plot_lons: ',plot_lons
            print 'plot_lats: ',plot_lats
            print 'fmt : ',fmt
            print 'linewidth : ',linewidth
            print 'min/max x[0] : ',np.amin(x[0]),np.amax(x[0])
            print 'min/max y[0] : ',np.amin(y[0]),np.amax(y[0])
            plt.plot([plot_lons[:]],[plot_lats[:]],fmt,linewidth=linewidth,transform=p[0])
#            plt.gca().plot(305,53,'ko',markersize=5)

    # Plot title

    if title is not None:
        title=textwrap.fill(title,70)
        plt.gcf().suptitle(title, fontsize=12, y=0.95)    
#        plt.gcf().text(x=0.5, y=0.92, text=title, fontsize=16 )    

    if outfile is not None:
        matplotlib.rcParams['font.size'] =8
        plt.savefig(outfile,dpi=200)
    else:
        plt.show()

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file",action="store",dest="filenames", nargs="+",
                    help="name(s) of input file(s) for scalar field 1, scalar field 2, vector component 1, vector component 2.")
    parser.add_argument("-s", "--scalar",action="store",dest="sca_names", nargs="+",
                    help="name(s) of one or two scalar field to plot. If two fields specified one will be plotted as colour filled contours and one as line contours.") 
    parser.add_argument("-v", "--vector",action="store",dest="vec_names", nargs="+",
                    help="names of two components of vector field to plot.") 
    parser.add_argument("-r", "--rec", action="store",dest="rec",type=int,default=0,nargs="+",
                    help="record number to plot (defaults to zero)")
    parser.add_argument("-l", "--level", action="store",dest="level",type=int,default=1,nargs="+",
                    help="model level to plot (for 3D fields)")
    parser.add_argument("-m", "--mask", action="store",dest="maskfile",default=None,nargs="+",
                    help="file containing mask array ('mask' or 'bathy' or 'bathymetry'; 0=land)")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",default=None,
                    help="plot to output file - filetype by extension. If unset plot to GUI.")
    parser.add_argument("-L", "--lineplot", action="store_true",dest="lineplot",
                    help="plot first scalar field with line contours instead of colour-filled contours.")
    parser.add_argument("-n", "--nlevs", action="store",dest="nlevs",type=float,default=15,nargs="+",
                    help="number of contour levels")
    parser.add_argument("-f", "--mnfld", action="store",dest="mnfld",type=float,default=None,nargs="+",
                    help="minimum field value to plot for colour filled contouring")
    parser.add_argument("-F", "--mxfld", action="store",dest="mxfld",type=float,default=None,nargs="+",
                    help="maximum field value to plot for colour filled contouring")
    parser.add_argument("-t", "--title", action="store",dest="title",default=None,
                    help="title for plot")
    parser.add_argument("-W", "--west", action="store",dest="west",type=float,default=None,
                    help="western limit of area to plot")
    parser.add_argument("-E", "--east", action="store",dest="east",type=float,default=None,
                    help="eastern limit of area to plot")
    parser.add_argument("-S", "--south", action="store",dest="south",type=float,default=None,
                    help="southern limit of area to plot")
    parser.add_argument("-N", "--north", action="store",dest="north",type=float,default=None,
                    help="northern limit of area to plot")
    parser.add_argument("-p", "--proj", action="store",dest="proj",default='pc0',
                    help="projection: pc0, pc180, mercator, northps, southps")
    parser.add_argument("-B", "--bottom", action="store",dest="bottom",type=str2bool,default=False,nargs="+",
                    help="plot bottom field")
    parser.add_argument("-e", "--scientific", action="store_true",dest="scientific",
                    help="scientific format for colour bar labels")
    parser.add_argument("-C", "--colors", action="store",dest="colors",default=None,nargs='+',
                    help="list of colors to use for line contouring - takes precedence over colour map in this case.")
    parser.add_argument("-c", "--cmap", action="store",dest="cmap",default=None,
                    help="colour map for contour plotting - inbuilt matplotlib colour map or external loaded with monty.clr_cmap")
    parser.add_argument("-R", "--reverse_colors", action="store_true",dest="reverse_colors",
                    help="reverse colour map table")
    parser.add_argument("-G", "--logscale", action="store_true",dest="logscale",
                    help="use colour logscale")
    parser.add_argument("-x", "--factor", action="store",dest="factor",type=float,default=None,nargs="+",
                    help="multiplicative factor to apply to plotted field")
    parser.add_argument("-Z", "--zeromean", action="store",dest="zeromean",type=str2bool,default=False,nargs="+",
                    help="Zero the spatial mean for the displayed area")
    parser.add_argument("-V", "--vec_only", action="store_true",dest="vec_only",
                    help="Only plot arrows for vector fields - don't contour the vector magnitude.")
    parser.add_argument("-M", "--mag_only", action="store_true",dest="mag_only",
                    help="Only plot magnitude of vector fields - don't plot arrows.")
    parser.add_argument("-A", "--arrows", action="store",dest="arrows",nargs="+",
                    help="arrow parameters : subsampling, scale, colour")
    parser.add_argument("-w", "--whitebg", action="store_true",dest="whitebg",
                    help="White background for plots. Don't apply stock image of land/ocean topography.")
    parser.add_argument("-d", "--draw_points", action="store",dest="draw_points",nargs="+",
                    help="list of points to draw line segments between in groups of four: start_lon,start_lat,end_lon,end_lat")
    parser.add_argument("-D", "--draw_fmt", action="store",dest="draw_fmt",nargs="+",
                    help="first argument is format for line segments as for fmt keyword for pyplot.plot; second optional argument is line thickness")

    args = parser.parse_args()

    plot_nemo(filenames=args.filenames,sca_names=args.sca_names,vec_names=args.vec_names,rec=args.rec,level=args.level,
              nlevs=args.nlevs,mnfld=args.mnfld,mxfld=args.mxfld,
              title=args.title,west=args.west,east=args.east,south=args.south,north=args.north,proj=args.proj,
              maskfile=args.maskfile,outfile=args.outfile,bottom=args.bottom,scientific=args.scientific,
              reverse_colors=args.reverse_colors,logscale=args.logscale,factor=args.factor,zeromean=args.zeromean,
              lineplot=args.lineplot,vec_only=args.vec_only,mag_only=args.mag_only,arrows=args.arrows,whitebg=args.whitebg,
              colors=args.colors,cmap=args.cmap,draw_points=args.draw_points,draw_fmt=args.draw_fmt)        
    
