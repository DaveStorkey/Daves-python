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

Nov 2018 : 1. Added block plotting option and rejigged control of plot types. DS.
           2. Added a proj="none" option. DS.

Feb 2019 : 1. Add option to do colour-filled and line plotting for the same field: 'cl' or 'lc'
           2. Limit the number of labels on the colour bar to max 10 (and only have 2 decimal places). 
           3. Use cartopy.feature.NaturalEarthFeature to draw the coastlines. 
           4. Add "extend='both'" to contourf calls. DS.

July 2019 : Conversion to python 3. DS.

Mar 2020  : New plot_type='lcf' or plot_type='clf' option to force lines for all contour levels in 'cl' plots.

Aug 2020  : New option "units" to override default. 

Nov 2020  : 'bl', 'lb', 'blf', 'lbf' options. And option for vertical color bar. DS.

Jan 2021  : 1. "nobar" option and return plot indices.
            2. More flexible treatment of figures (can adjust size) and axes (allow subplots).
            DS.

May 2021  : Add lat/lon labels and set facecolor for no projection option. DS. 

Aug 2021  : Bug fix for the case where you want to plot vector magnitude in colour + arrows. DS.

Feb 2022  : nocoast option. DS.
            fix limited area plotting for "none" projection. DS. 

To do    : 1. Add lat/lon tickmarks (look at Ocean Assess)
           2. Size the figure in a clever way depending on the aspect ratio of the input data. DS. 


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
import cartopy as ctpy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import numpy as np
import numpy.ma as ma
#import monty
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

def project_cube(cube_in, proj=None, bounds=False):

    # NB. iris.analysis.cartography.project doesn't actually do a projection (ho hum), it does an *interpolation*
    #     to a regular grid in the desired projection. (nx,ny) specify the dimensions of the *global* target grid. 
    #     We need to do some jiggery pokery to try to make sure the target grid is a similar resolution to our 
    #     input data since the input data might not be global. 

    lon_extent = np.amax(cube_in.coord('longitude').points) - np.amin(cube_in.coord('longitude').points)
    nx = np.int_(cube_in.shape[-1] * 360.0/lon_extent)
    lat_extent = np.amax(cube_in.coord('latitude').points) - np.amin(cube_in.coord('latitude').points)
    ny = np.int_(cube_in.shape[-2] * 180.0/lat_extent)

    print('lon_extent, cube_in.shape[-1], nx : ',lon_extent, cube_in.shape[-1], nx)
    print('lat_extent, cube_in.shape[-2], ny : ',lat_extent, cube_in.shape[-2], ny)

    cube_out, extent = iris.analysis.cartography.project(cube_in, proj, nx=nx, ny=ny)
    x = cube_out.coord('projection_x_coordinate').points
    y = cube_out.coord('projection_y_coordinate').points
    if bounds:
        # need to guess the coordinate bounds because Iris doesn't do bounds
        # for 2D auxillary coordinates (even though they are in the files).
        xbounds,ybounds = guess_bounds(x,y)
        return x,y,cube_out,xbounds,ybounds
    else:
        return x,y,cube_out

def get_coord(cube, coordname, filename):

#   Lats/lons will be recognised by Iris as auxillary coordinates if the "coordinates"
#   attribute of the field in the file has been set correctly. If not, the lat/lon fields
#   might still be in the file, so try reading them as simple fields rather than coordinates.
    try:
        coord = cube.coord(coordname).points.copy()
    except iris.exceptions.CoordinateNotFoundError:
        print('reading '+coordname+' as field')
        coord=iris.load_cube(filename,coordname).data
    else:
        print('reading '+coordname+' as auxillary coordinate')
    # force lat/lon coordinates to be 2D arrays
    if len(coord.shape) == 1:
        try:
            coord_out = coord * np.ones(cube.shape)
        except ValueError:
            coord_out = (coord * np.ones(cube.shape).transpose()).transpose()
    else:
        coord_out = coord
    return coord_out

def guess_bounds(x,y):

    # get approximate bounds on the x and y arrays by interpolating to midpoints and extrapolating at end points.
    if len(x.shape) == 1:
        xbounds = np.zeros(x.shape[0]+1)
        xbounds[1:-1] = 0.5 * ( x[:-1] + x[1:] )
        xbounds[0] = x[0] - (xbounds[1] - x[0])
        xbounds[-1] = x[-1] + (x[-1] - xbounds[-2])
    elif len(x.shape) == 2:
        xbounds = np.zeros((x.shape[0]+1,x.shape[1]+1))
        xbounds[1:,1:-1] = 0.5 * ( x[:,:-1] + x[:,1:] )
        xbounds[1:,0] = x[:,0] - (xbounds[1:,1] - x[:,0])
        xbounds[1:,-1] = x[:,-1] + (x[:,-1] - xbounds[1:,-2])
        xbounds[0,:] = xbounds[1,:]
    else:
        raise Exception("Error: Can only deal with 1D or 2D coordinates arrays.") 

    if len(y.shape) == 1:
        ybounds = np.zeros(y.shape[0]+1)
        ybounds[1:-1] = 0.5 * ( y[:-1] + y[1:] )
        ybounds[0] = y[0] - (ybounds[1] - y[0])
        ybounds[-1] = y[-1] + (y[-1] - ybounds[-2])
    elif len(y.shape) == 2:
        ybounds = np.zeros((y.shape[0]+1,y.shape[1]+1))
        ybounds[1:-1,1:] = 0.5 * ( y[:-1,:] + y[1:,:] )
        ybounds[0,1:] = y[0,:] - (ybounds[1,1:] - y[0,:])
        ybounds[-1,1:] = y[-1,:] + (y[-1,:] - ybounds[-2,1:])
        ybounds[:,0] = ybounds[:,1]
    else:
        raise Exception("Error: Can only deal with 1D or 2D coordinates arrays.")

    return xbounds,ybounds

def plot_nemo(filenames=None,sca_names=None,vec_names=None,nlevs=13,mnfld=None,mxfld=None,title=None,rec=0,level=1,bottom=False,
              scientific=False,cmap=None,colors=None,reverse_colors=False,glob=None,west=None,east=None,south=None,north=None,
              proj=None,maskfile=None,outfile=None,logscale=None,factor=None,plot_types=None,zeromean=False,arrows=None,
              whitebg=False,noshow=None,units=None,vertbar=None,nobar=None,figx=None,figy=None,subplot=None,nocoast=None,
              draw_points=None,draw_fmt=None,fontsizes=None):

    # short cuts
    projections = { 'none'          : ( None               , 'global' ),
                    'pc0'           : ( ccrs.PlateCarree() , 'global' ),
                    'pc180'         : ( ccrs.PlateCarree(central_longitude=180.0), 'global' ), 
                    'merc'          : ( ccrs.Mercator(min_latitude=-86.0), 'global' ), 
                    'northps'       : ( ccrs.NorthPolarStereo(), [0,359.9,50,90] ), 
#                    'northps'       : ( ccrs.NorthPolarStereo(), [0,359.9,60,90] ), 
                    'southps_large' : ( ccrs.SouthPolarStereo(), [0,359.9,-90,-35] ), 
                    'southps'       : ( ccrs.SouthPolarStereo(), [0,359.9,-90,-45] ), 
                    'southps_zoom'  : ( ccrs.SouthPolarStereo(), [0,359.9,-90,-60] ),
                    'south_stereo'  : ( ccrs.Stereographic(central_latitude=-90.0, central_longitude=0.0),[0,359,-90,-60] )
    }   

    if proj is None:
        proj = 'pc0'

    p = projections[proj]
        
    nscalar = 0
    nvector = 0
    regrid_vector = False
    filenames_to_read = []
    varnames_to_read = []
    if filenames is not None and type(filenames) is not list:
        filenames = [filenames]
    elif filenames is None or len(filenames) > 4:
        raise Exception('Error: specify between one and four filenames.')
    if sca_names is not None and type(sca_names) is not list:
        sca_names = [sca_names]
    if vec_names is not None and type(vec_names) is not list:
        vec_names = [vec_names]
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
        print('Reading '+varname+' from '+filename)
        fld_in.append(gt.read_cube(filename,varname))
        
    print('number of scalar fields to plot : ',nscalar)
    print('number of vector fields to plot : ',nvector)
    print('len(fld_in) :', len(fld_in))
    print('regrid_vector is : ',regrid_vector)

    fontsizes_list = [14,12,10]
    for count,fontsize in enumerate(fontsizes):
        if fontsize is not None:
            fontsizes_list[count]=fontsize
    fontsizes = fontsizes_list

# Organise control parameters that apply to individual fields.
# In general with multiple fields and a single value of the control
# parameter then that value gets applied to all fields, otherwise
# the number of values input has to match the number of fields, where
# the vector field counts as a single field. 

    nfields = nscalar + nvector
    if plot_types is None:
        plot_types=[]
        if nscalar == 1:
            plot_types = plot_types+['c']
        elif nscalar == 2:
            plot_types = plot_types+['c','l']
        if nvector:
            plot_types = plot_types+['a']
    else:
        if type(plot_types) is not list:
            plot_types = [plot_types]
        if len(plot_types) != nfields:
            raise Exception("Error: number of values for plot_types (-p) doesn't match number of fields.")
    
    vec_only = False
    mag_only = False
    if nvector:
        # need to make sure that plot_types has the same number of elements
        # as fldslice otherwise a zip goes wrong later on :)
        plot_types = plot_types + [plot_types[-1]]
        if plot_types[v1] == 'a':
            vec_only = True
        elif plot_types[v1] == 'c' or plot_types[v1] == 'b' or plot_types[v1] == 'l':
            mag_only = True

    if type(rec) is not list:
        rec = (nscalar + 2*nvector) * [rec]
    elif len(rec) != nfields:
        raise Exception("Error: number of values for REC (-r) doesn't match number of fields.")
    else:
        rec = rec + nvector*[rec[-1]]
    if type(level) is not list:
        level = (nscalar + 2*nvector) * [level]
    elif len(level) == 1:
        level = (nscalar + 2*nvector) * level
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
        mnfld = (nscalar + 2*nvector) * [mnfld]
    elif len(mnfld) == 1:
        mnfld = nfields*mnfld
    elif len(mnfld) != nfields:
        raise Exception("Error: number of values for MNFLD (-f) doesn't match number of fields.")
    if type(mxfld) is not list:
        mxfld = (nscalar + 2*nvector) * [mxfld]
    elif len(mxfld) == 1:
        mxfld = nfields*mxfld
    elif len(mxfld) != nfields:
        raise Exception("Error: number of values for MXFLD (-F) doesn't match number of fields.")
    if type(bottom) is not list:
        bottom = (nscalar + 2*nvector) * [bottom]
    elif len(bottom) != nfields:
        raise Exception("Error: number of values for BOTTOM (-b) doesn't match number of fields.")
    else:
        bottom = bottom + nvector*[bottom[-1]]
    if type(nlevs) is not list:
        nlevs = (nscalar + 2*nvector) * [nlevs]
    elif len(nlevs) == 1:
        nlevs = nfields*nlevs
    elif len(nlevs) != nfields:
        raise Exception("Error: number of values for NLEVS (-n) doesn't match number of fields.")
    if type(factor) is not list:
        factor = (nscalar + 2*nvector) * [factor]
    elif len(factor) != nfields:
        raise Exception("Error: number of values for FACTOR (-x) doesn't match number of fields.")
    else:
        factor = factor + nvector*[factor[-1]]
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
            print('Extracting bottom field for '+fldi.var_name)
            fldslice.append(fldi[0].copy())
#            fldslice[-1][:] = 0.0
            if fldi.var_name == 'vovecrtz' or fldi.var_name == 'wo':
                # for vertical velocities take the level just above the bottom (because the bottom field is zero by construction).
                for lev in range(len(fldi.data[:,0,0])-1):
                    fldslice[-1].data = np.ma.where( (~fldi.data[lev,:,:].mask & fldi.data[lev+1,:,:].mask), 
                                             x=fldi.data[lev-1,:,:], y=fldslice.data )
            else:
                print('fldi[0].data.mask.shape : ',fldi[0].data.mask.shape)
                for lev in range(len(fldi.data[:,0,0])-1):
                    fldslice[-1].data[:,:] = np.ma.where( (~fldi.data.mask[lev,:,:] & fldi.data.mask[lev+1,:,:]), 
                                                  x=fldi.data[lev-1,:,:], y=fldslice[-1].data[:,:] )

# ... or select required level
        else:
            if fldi.ndim == 3:
               print ('Extracting level ',level_i,' from field.')
               fldslice.append(fldi[level_i-1])
            else:
               fldslice.append(fldi)


        print('min/max fldslice : ',fldslice[-1].data.min(),fldslice[-1].data.max())

    # Do interpolation and rotation of vector fields *after* we have 
    # selected a single record and single level - much more efficient!

    # interpolate to T-points:
    if regrid_vector and not bottom[v1]:
        fldslice[v1],fldslice[v2] = UaT.U_at_T_points(fldslice[v1],fldslice[v2])

    # rotate to be relative to geographical coordinates rather than model grid:
    # (unnecessary if we are only plotting the magnitude of the vector field).
    if nvector and not mag_only and proj != 'none':
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
        print('min/max lons_i after normalisation: ',lons_i.min(),lons_i.max())

    if p[1] == 'global':
        (west_default,east_default,south_default,north_default) = (0.0,359.99,-90.0,90.0)
    else:
        (west_default,east_default,south_default,north_default) = p[1]

    if glob:
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

    print('region : ',region)

    # Mask out or cut out bits of the field(s) to be contoured outwith the area being plotted 
    # so that the automatic min/max calculation works. Also so that the zero mean calculation is restricted to
    # the plotted area. Also so we can calculate the maximum vector length for the key arrow.

    fldslice_cutout=[]
    if region != 'global':
        for lons_i, lats_i, fldslice_i in zip(lons[:],lats[:],fldslice[:]): 
            fldslice_i.data = ma.masked_where( ((lons_i < region[0]) | (lons_i > region[1])) , fldslice_i.data )
            fldslice_i.data = ma.masked_where( ((lats_i < region[2]) | (lats_i > region[3])) , fldslice_i.data )
            if proj == 'none':
                # cut out the field in this case because ax.set_extent doesn't work with no projection
                # Note that here we can't just redefine fldslice_i as the cutout field because fldslice_i
                # is only a view on the original array, which isn't modified in this case. Instead we have
                # to create a new list of arrays fldslice_cutout and make the substition after the end of the loop.
                ny,nx = fldslice_i.data.shape
                ii = ma.arange(nx) * ma.ones((ny,nx))
                jj = (ma.arange(ny) * ma.ones((nx,ny))).transpose()
                ii = ma.masked_where( ((lons_i < region[0]) | (lons_i > region[1])) , ii ) 
                ii = ma.masked_where( ((lats_i < region[2]) | (lats_i > region[3])) , ii ) 
                jj = ma.masked_where( ((lons_i < region[0]) | (lons_i > region[1])) , jj ) 
                jj = ma.masked_where( ((lats_i < region[2]) | (lats_i > region[3])) , jj ) 
                fldslice_cutout.append(fldslice_i[int(jj.min()):int(jj.max())+1,int(ii.min()):int(ii.max()+1)])
            else:
                fldslice_cutout.append(fldslice_i)

    fldslice[:] = fldslice_cutout[:]

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
                print('fldslice_i.standard_name, zeromean_i : ',fldslice_i.standard_name, zeromean_i)
                if zeromean_i:
                    print("Reading cell_area from "+filename_to_read)
                    try:
                        cell_area = iris.load_cube(filename_to_read,'cell_area')
                    except(iris.exceptions.ConstraintMismatchError):
                        raise Exception('Could not read cell_area from file.\n'+
                                        'For the zero mean option you need a field called cell_area \n'+
                                        'in the file with the correct coordinates attribute set.')
                    area_mean = fldslice_i.collapsed(['longitude','latitude'],iris.analysis.MEAN,weights=cell_area.data)
                    print('area_mean.data : ',area_mean.data)
                    fldslice_i.data = fldslice_i.data - area_mean.data

    # Calculate magnitude of vector field if there is one. If we are only plotting the vector field as arrows then 
    # this is used later to calculate the length of the key arrow. If we are also/instead contouring 
    # the magnitude of the vector field (only done if no scalar fields to plot) then add the magnitude
    # field slice to the beginning of the fldslice list.
     
    if nvector:
        u = fldslice[v1].data
        v = fldslice[v2].data
        speed = ma.sqrt(u*u + v*v)
        if not nscalar and not vec_only:
            # Append vector magnitude to the beginning of the fldslice list.
            nscalar=1; v1=1; v2=2;
            plot_types = [plot_types[v1][0]]+plot_types[:]
            fldslice = [fldslice[0].copy()] + fldslice
            # hardwired for now - doesn't really matter
            fldslice[0].standard_name = "sea_water_speed"
            fldslice[0].data = speed[:]

    # Calculate contour levels for contour plots and colour map for colour-filled contouring.

    if nscalar or not vec_only:
        levs = []
        for mnfld_i,mxfld_i,nlevs_i,fldslice_i,plot_type in zip(mnfld[:nscalar],mxfld[:nscalar],nlevs[:nscalar],fldslice[:nscalar],plot_types[:nscalar]):
            if mnfld_i is None or mnfld_i == "None":
                mnfld_i = ma.min(fldslice_i.data)
            else:
                mnfld_i = float(mnfld_i)
            if mxfld_i is None or mxfld_i == "None":
                mxfld_i = ma.max(fldslice_i.data)
            else:
                mxfld_i = float(mxfld_i)

            print('nlevs_i : ',nlevs_i)
            if logscale and ( 'c' in plot_type or 'b' in plot_type ):
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
            if cmap is not None and ('c' in plot_type or 'b' in plot_type):
#                if cmap == 'Rainbow':
#                    # shortcut for old IDL rainbow colour map:
#                    cmap = monty.clr_cmap('/home/h04/frsy/IDL/rainbow_diff_nice.clr')
#                    cmap.set_under(cmap(0.01))
#                    cmap.set_over(cmap(0.99))
#                elif '/' in cmap:
#                    # this is pointing to a file so use monty.clr_cmap:                
#                    print('reading colour map using monty.clr_cmap')
#                    cmap = monty.clr_cmap(cmap)
#                else:
                # assume that it refers to a matplotlib inbuilt colour map:
                print('reading inbuilt matplotlib colour map')
                cmap = getattr(matplotlib.cm,cmap)
                if reverse_colors:
                    cmap = mcolors.ListedColormap(cmap(np.arange(255,-1,-1)))

#                ind = np.intp(np.rint(np.linspace(2,255,nlevs_i+1))).tolist()    
# ensure we have white space in the middle for an odd number of levels.
#                if 2*(nlevs_i/2) != nlevs_i:
#                    ind[nlevs_i/2] = 127
#                print('ind',ind)
#                cmap=cmap(ind)

        print('len(levs) : ',len(levs))

    # Create new "projected" cubes. Shouldn't really have to do this but need
    # it a workaround for a cartopy bug. Note that the fields are all on the
    # same grid, so (x,y) is the same for each field.

    fldproj = []
    x=[]
    y=[]
    xbounds=[]
    ybounds=[]

    if p[0] is None:
        # no projection
        for fldslice_i, plot_type in zip(fldslice,plot_types):
            print('fldslice_i.shape : ',fldslice_i.shape)
            fldproj.append(fldslice_i)
            x.append(fldslice_i.coord('longitude').points)
            y.append(fldslice_i.coord('latitude').points)
            if 'b' in plot_type:
                xbounds_i,ybounds_i = guess_bounds(x[-1],y[-1])
                xbounds.append(xbounds_i)
                ybounds.append(ybounds_i)
            else:
                xbounds.append(None)
                ybounds.append(None)
            print('min/max x : ',np.amin(x),np.amax(x))
            print('min/max y : ',np.amin(y),np.amax(y))
    else:
        for fldslice_i, plot_type in zip(fldslice,plot_types):
            if 'b' in plot_type:
                projinfo = project_cube(fldslice_i, proj=p[0], bounds=True)
            else:
                projinfo = project_cube(fldslice_i, proj=p[0], bounds=False)
            x.append(projinfo[0])
            y.append(projinfo[1])
            fldproj.append(projinfo[2])
            if 'b' in plot_type:
                xbounds.append(projinfo[3])
                ybounds.append(projinfo[4])
            else:
                xbounds.append(None)
                ybounds.append(None)

    if figx is None:
        figx = matplotlib.rcParams["figure.figsize"][0]
    if figy is None:
        figy = matplotlib.rcParams["figure.figsize"][1]

    # If figure 1 already exists then plt.figure will 
    # just return a reference to it.
    fig = plt.figure(1, figsize=[figx,figy])

    # Create plot axes:

    if subplot is not None:
        if len(subplot) == 1:
            ax = fig.add_subplot(subplot,projection=p[0])
        elif len(subplot) == 3:
            ax = fig.add_subplot(subplot[0],subplot[1],subplot[2],projection=p[0])
        else:
            raise Exception("Error : len(subplot) must be 1 or 3.")
    else:
        ax = plt.axes(projection=p[0])

    if p[0] is None:
        # In this case it will label lat and lon values along the axes, so add axis labels
        ax.set_xlabel("longitude (deg east)")
        ax.set_ylabel("latitude (deg north)")
        if whitebg:
            ax.set_facecolor('1')
        else:
            # default to grey over land
            ax.set_facecolor('0.75')
    else:
        if whitebg:
            facecolor='1'
        else:
            # default to grey over land
            facecolor='0.75'
        if nocoast:
            ax.set_facecolor(facecolor)
        else:
            feature = ctpy.feature.NaturalEarthFeature('physical', 'coastline', '50m',facecolor=facecolor,edgecolor='k')
            ax.add_feature(feature,linewidth=0.5)
        # set the plot extent
        if region == 'global':
            ax.set_global()
        else:
            ax.set_extent(region, crs=ccrs.PlateCarree())
            # ax.set_extent(region, crs=ccrs.Geodetic())

        if p[0] == ccrs.PlateCarree() or p[0] == ccrs.PlateCarree(central_longitude=180.0):
            # Axis formatting - ONLY WORKS FOR PLATECARREE (or no projection)!!
            gl = plt.gca().gridlines(crs=p[0], draw_labels=True)
            gl.xlines = gl.ylines = gl.xlabels_top = gl.ylabels_right = False
            gl.xlabel_style = gl.ylabel_style = {'size':fontsizes[2], 'color':'gray'}
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER

    csline=None
    cscolor=None
    csarrows=None

    if nscalar:
        for x_i,xbounds_i,y_i,ybounds_i,fldproj_i,plot_type,levs_i,fld_in_i in \
          zip(x[:nscalar],xbounds[:nscalar],y[:nscalar],ybounds[:nscalar],fldproj[:nscalar],plot_types[:nscalar],levs[:nscalar],fld_in[:nscalar]):
            print('levs_i : ',levs_i)
            if plot_type == 'l':
                ### contour lines ###
                csline=plt.contour(x_i,y_i,fldproj_i.data,levels=levs_i,mxfld=levs_i[-1],mnfld=levs_i[0],colorbar=None,colors=colors,linewidths=0.5)
            elif plot_type == 'c':
                ### colour-filled contours ###
                cscolor=plt.contourf(x_i,y_i,fldproj_i.data,levels=levs_i,mxfld=levs_i[-1],mnfld=levs_i[0],colorbar=None,cmap=cmap,extend='both')
            elif 'cl' in plot_type or 'lc' in plot_type:
                ### colour-filled contours and lines for the same field ###
                cscolor=plt.contourf(x_i,y_i,fldproj_i.data,levels=levs_i,mxfld=levs_i[-1],mnfld=levs_i[0],colorbar=None,cmap=cmap,extend='both')
                if len(levs_i) > 10 and 'f' not in plot_type:
                    levs_line = levs_i[::2]
                else:
                    levs_line = levs_i
                csline=plt.contour(x_i,y_i,fldproj_i.data,levels=levs_line,mxfld=levs_i[-1],mnfld=levs_i[0],colorbar=None,colors="black",linewidths=0.5)
            elif plot_type == 'b':
                ### block plot ###
                cscolor = plt.pcolormesh(xbounds_i, ybounds_i, fldproj_i.data, cmap=cmap, vmin=levs_i[0], vmax=levs_i[-1])
            elif 'bl' in plot_type or 'lb' in plot_type:
                ### block plot and lines for the same field ###
                cscolor = plt.pcolormesh(xbounds_i, ybounds_i, fldproj_i.data, cmap=cmap, vmin=levs_i[0], vmax=levs_i[-1])
                if len(levs_i) > 10 and 'f' not in plot_type:
                    levs_line = levs_i[::2]
                else:
                    levs_line = levs_i
                csline=plt.contour(x_i,y_i,fldproj_i.data,levels=levs_line,mxfld=levs_i[-1],mnfld=levs_i[0],colorbar=None,colors="black",linewidths=0.5)
            else:
                raise Exception("Error: Something has gone wrong. Scalar plot_type should be 'l' or 'c' or 'b' or 'cl' or 'lc'")

            ### colour bar ###
            if ('c' in plot_type or 'b' in plot_type) and not nobar:
                if vertbar:
                    cax=plt.colorbar(cscolor,orientation='vertical')
                else:
                    cax=plt.colorbar(cscolor,orientation='horizontal')
                if len(levs_i) > 11:
                    levs_ticks = np.linspace(levs_i[0],levs_i[-1],11)
                else:
                    levs_ticks = levs_i
                cax.set_ticks(levs_ticks)
                if scientific:
                    labels=["%1.2e" %lev for lev in levs_ticks]
                else:
                    labels=["%1.2f" %lev for lev in levs_ticks]
                if vertbar:
                    cax.ax.set_yticklabels(labels,fontsize=fontsizes[2])
                else:
                    cax.ax.set_xticklabels(labels,fontsize=fontsizes[2],rotation=45)
                if units is not None:
                    # units keyword overrides units read in from file
                    cax.ax.set_xlabel(units,fontsize=fontsizes[1])
                elif str(fld_in_i.units) != "1":
                    cax.ax.set_xlabel(str(fld_in_i.units),fontsize=fontsizes[1])

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
                arrow_scale = float(arrows[1])
            if len(arrows) > 2:
                arrow_colour=arrows[2]
            if len(arrows) > 3:
                key_arrow_length=arrows[3]
        print("Plotting arrows")
        csarrows = plt.quiver(x[v1][::arrow_subsample], y[v1][::arrow_subsample], 
                   fldproj[v1].data[::arrow_subsample,::arrow_subsample], fldproj[v2].data[::arrow_subsample,::arrow_subsample], 
                   color=arrow_colour, scale=arrow_scale, scale_units='inches', angles='xy' )
        if scientific:
            arrow_label = '%1.1e'%key_arrow_length+str(fld_in[v1].units)
        else:        
            arrow_label = '%1.1f'%key_arrow_length+str(fld_in[v1].units)
        plt.quiverkey(csarrows, X=0.8, Y=-0.05, U=key_arrow_length, label=arrow_label, labelpos='E')

    # Line segments

    if draw_points is not None:
        if len(draw_points)%4 != 0:
            raise Exception("Error: number of draw_points (-d) must be a multiple of 4: start_x,start_depth,end_x,end_depth")
        fmt = "k-"  # default to black solid lines        
        linewidth = 2
        if draw_fmt is not None:
            fmt = draw_fmt[0]
            if len(draw_fmt) == 2:
                linewidth = draw_fmt[1] 
        for ii in range(0,len(draw_points),4):
            # pyplot.plot takes all the x-values first and the y-values second...
            plot_x = [draw_points[ii],draw_points[ii+2]]
            plot_y = [draw_points[ii+1],draw_points[ii+3]]
            print("Plotting : ",plot_x,plot_y," on projection ",p[0])
            plt.gca().plot(plot_x[:],plot_y[:],fmt,linewidth=linewidth,transform=p[0])
            

    # Plot title

    if title is not None:
        title=textwrap.fill(title,70)
        if subplot is not None:
            plt.gca().set_title(title, fontsize=fontsizes[0])    
        else:
            plt.gcf().suptitle(title, fontsize=fontsizes[0], y=0.95)    

    if outfile is not None:
        matplotlib.rcParams['font.size'] =8
        plt.savefig(outfile,dpi=200)
    elif not noshow:
        plt.show()

    # need this if we are calling plot_nemo multiple times from a python script
    # otherwise the second call will start working on the same figure as the first call etc.
    if not noshow:
        plt.close()

    # Return references to the plots for use in calling python routines. In general one or more of these might have None values. 
    return csline, cscolor, csarrows

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file",action="store",dest="filenames", nargs="+",
                    help="name(s) of input file(s) for scalar field 1, scalar field 2, vector component 1, vector component 2.")
    parser.add_argument("-s", "--scalar",action="store",dest="sca_names", nargs="+",
                    help="name(s) of one or two scalar field to plot.") 
    parser.add_argument("-v", "--vector",action="store",dest="vec_names", nargs="+",
                    help="names of two components of vector field to plot.") 
    parser.add_argument("-p", "--plot_types", action="store",dest="plot_types",nargs="+",
                    help="plot types for each field: 'l' = line contours, 'c' = colour-filled contours, 'b' = block plot, 'a' = arrows only, 'ca' = colour filled contours and arrows, 'ba' = block plot and arrows")
    parser.add_argument("-r", "--rec", action="store",dest="rec",type=int,default=0,nargs="+",
                    help="record number to plot (defaults to zero)")
    parser.add_argument("-l", "--level", action="store",dest="level",type=int,default=1,nargs="+",
                    help="model level to plot (for 3D fields). Top level = 1.")
    parser.add_argument("-m", "--mask", action="store",dest="maskfile",default=None,nargs="+",
                    help="file containing mask array ('mask' or 'bathy' or 'bathymetry'; 0=land)")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",default=None,
                    help="plot to output file - filetype by extension. If unset plot to GUI.")
    parser.add_argument("-n", "--nlevs", action="store",dest="nlevs",type=int,default=15,nargs="+",
                    help="number of contour levels")
    parser.add_argument("-f", "--mnfld", action="store",dest="mnfld",default=None,nargs="+",
                    help="minimum field value(s) to plot for colour filled contouring - specify None for default")
    parser.add_argument("-F", "--mxfld", action="store",dest="mxfld",default=None,nargs="+",
                    help="maximum field value(s) to plot for colour filled contouring - specify None for default")
    parser.add_argument("-t", "--title", action="store",dest="title",default=None,
                    help="title for plot")
    parser.add_argument("-Y", "--fontsizes", action="store",dest="fontsizes",nargs='+',
                    help="fontsizes for : plot title, axis labels, axis tick labels. Use None for default.")
    parser.add_argument("-G", "--global", action="store_true",dest="glob",default=False,
                    help="force global limits on plot even if data is limited area")
    parser.add_argument("-W", "--west", action="store",dest="west",type=float,default=None,
                    help="western limit of area to plot")
    parser.add_argument("-E", "--east", action="store",dest="east",type=float,default=None,
                    help="eastern limit of area to plot")
    parser.add_argument("-S", "--south", action="store",dest="south",type=float,default=None,
                    help="southern limit of area to plot")
    parser.add_argument("-N", "--north", action="store",dest="north",type=float,default=None,
                    help="northern limit of area to plot")
    parser.add_argument("-P", "--proj", action="store",dest="proj",default='pc0',
                    help="projection: pc0, pc180, mercator, northps, southps")
    parser.add_argument("-B", "--bottom", action="store",dest="bottom",type=str2bool,default=False,nargs="+",
                    help="plot bottom field")
    parser.add_argument("-e", "--scientific", action="store_true",dest="scientific",
                    help="scientific format for colour bar labels")
    parser.add_argument("-C", "--colors", action="store",dest="colors",default=None,nargs='+',
                    help="list of colors to use for line contouring.")
    parser.add_argument("-c", "--cmap", action="store",dest="cmap",default=None,
                    help="colour map for filled contour or block plotting - inbuilt matplotlib colour map or external loaded with monty.clr_cmap")
    parser.add_argument("-R", "--reverse_colors", action="store_true",dest="reverse_colors",
                    help="reverse colour map table")
    parser.add_argument("-g", "--logscale", action="store_true",dest="logscale",
                    help="use colour logscale")
    parser.add_argument("-V", "--vertbar", action="store_true",dest="vertbar",
                    help="use vertical color bar (default horizontal)")
    parser.add_argument("--nobar", action="store_true",dest="nobar",
                    help="Suppress colour bar.")
    parser.add_argument("-x", "--factor", action="store",dest="factor",type=float,default=None,nargs="+",
                    help="multiplicative factor to apply to plotted field")
    parser.add_argument("-Z", "--zeromean", action="store",dest="zeromean",type=str2bool,default=False,nargs="+",
                    help="Zero the spatial mean for the displayed area")
    parser.add_argument("-A", "--arrows", action="store",dest="arrows",nargs="+",
                    help="arrow parameters : subsampling, scale, colour")
    parser.add_argument("-u", "--units", action="store",dest="units",default=None,
                    help="units (label for colour bar)")
    parser.add_argument("-w", "--whitebg", action="store_true",dest="whitebg",
                    help="White background for plots. Don't apply stock image of land/ocean topography.")
    parser.add_argument("--noshow", action="store_true",dest="noshow",
                    help="Don't run plt.show() even if no output file specified (for calls from other python routines).")
    parser.add_argument("--nocoast", action="store_true",dest="nocoast",
                    help="Don't draw coastline.")
    parser.add_argument("--figx", action="store",dest="figx",default=None,
                    help="x-dimension of figure (in inches I think)")
    parser.add_argument("--figy", action="store",dest="figy",default=None,
                    help="y-dimension of figure (in inches I think)")
    parser.add_argument("--draw_points", action="store",dest="draw_points",nargs="+",
                    help="list of points to draw line segments between in groups of four: start_lon,start_lat,end_lon,end_lat")
    parser.add_argument("--draw_fmt", action="store",dest="draw_fmt",nargs="+",
                    help="first argument is format for line segments as for fmt keyword for pyplot.plot; second optional argument is line thickness")

    args = parser.parse_args()

    plot_nemo(filenames=args.filenames,sca_names=args.sca_names,vec_names=args.vec_names,rec=args.rec,level=args.level,
              nlevs=args.nlevs,mnfld=args.mnfld,mxfld=args.mxfld,
              title=args.title,glob=args.glob,west=args.west,east=args.east,south=args.south,north=args.north,proj=args.proj,
              maskfile=args.maskfile,outfile=args.outfile,bottom=args.bottom,scientific=args.scientific,
              reverse_colors=args.reverse_colors,logscale=args.logscale,factor=args.factor,zeromean=args.zeromean,
              plot_types=args.plot_types,arrows=args.arrows,whitebg=args.whitebg,vertbar=args.vertbar,
              colors=args.colors,cmap=args.cmap,noshow=args.noshow,units=args.units,nobar=args.nobar,fontsizes=args.fontsizes,
              figx=args.figx,figy=args.figy,nocoast=args.nocoast,draw_points=args.draw_points,draw_fmt=args.draw_fmt)        
    
