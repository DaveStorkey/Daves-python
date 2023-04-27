#! /usr/bin/env python

'''
Routine to plot the cross section of fields on
an irregular horizontal grid. 

July 2017 : Set all defaults in function definition and remove from argparse. Better 
            practice because it means it is safe to call the function from within
            python bypassing the argparse bit. Also add noshow and nobar flags. DS. 

May 2018  : Added option (-B) to do block plotting instead of filled contour plotting.
            Also option to plot a section along a gridline (method='gridline'). In this
            case "lat" is interpreted as the y-index or "lon" as the x-index for the section
            and xmin,xmax and the min/max x/y index along the section. DS.

Oct 2018  : 0. Defaults set in the body of the routine. Setting defaults in the "def" statement doesn't
               work if you use argparse because argparse arguments default to None and override the 
               defaults in the "def" statement. 
            1. Enable line contouring and overplotting a line contour and a colour-filled contour plot.
            2. Make sure a line coutour plot has land greyed out and fill the field prior to contouring
               to avoid gaps next to bathymetry.
            3. Sorted out masking properly, including the -z option to mask out zero values which 
               you get in NEMO U-grid and V-grid fields. 
            4. Got block plotting (pcolormesh) working properly using the coordinate bounds.
            5. For "gridline" method you can now specify things in terms of the (i,j) indices or in terms
               of (lon,lat). In the latter case it will find the nearest gridline.  

Nov 2018    1. Create the get_section routine so I can call it from other modules. After some iteration
               this returns indices that you can use to extract the section from the field. 
            2. Put in drawing option for lines and points. 
            3. Option to input xsec_indices rather than calculating them. 
            4. Rejig control of plot types slightly.

June 2018   Update to python 3.

Feb 2020  : 1. Enable plot_nemo_section to read in a section file created by cdf_xtrac_brokenline. (Eventually
               would be good if plot_nemo_section could call cdf_xtrac_brokenline to create the section). 
            2. Units on the color bar.

July 2020 : Hack to fill in latitude values for north-south sections that include the "patched" areas of the
            extended ORCA grids. DS.

Nov 2021 : "units" and "vertbar" options copied from plot_nemo.py. DS.

Nov 2021 : Fix for case where no masked values on a section. DS

Aug 2022 : Add subplot, figx, figy options and return the plot objects so can be called from a 
           wrapper routine to do multi-plots. DS.

To do : Currently we don't take account of partial cells at the bottom. In principle could do this by reading in the
        3D e3t field, but I think getting the masking to work properly could be quite tricky.  

@author: Dave Storkey and Pat Hyder
'''

import socket
import matplotlib
if 'spice' in socket.gethostname():
    # Note this disables plt.show()
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.cm as cm
import iris
import numpy as np
import numpy.ma as ma
import general_tools as gt
from pylab import *
import textwrap
# import glob

def get_coord(cube, coordname, filename):
    #   Lats/lons will be recognised by Iris as auxillary coordinates if the "coordinates"
    #   attribute of the field in the file has been set correctly. If not, the lat/lon fields
    #   might still be in the file, so try reading them as simple fields rather than coordinates.
    try:
        coord = cube.coord(coordname).points
        coord_bounds = cube.coord(coordname).bounds
    except iris.exceptions.CoordinateNotFoundError:
        print('reading '+coordname+' as field')
        coord = iris.load_cube(filename,coordname).data
        # No easy way to get the bounds in this case.
        # Would need to look at the "bounds" attribute and then
        # read this in using the netCDF4 library rather than Iris.
        coord_bounds = None
    else:
        print('reading '+coordname+' as auxillary coordinate')
    return coord, coord_bounds

def fill_in_south(xpoints):
    # Switchable hack to deal with the case that we want to do a section that includes some of the "patched"
    # areas in Antarctica for the extended ORCA grids. In this case there are areas over land where the lat/lon
    # values have been filled with zeroes which messes up the plotting. Just extrapolate south from where we have 
    # sensible values. 
    count = 0
    for xpoint in xpoints:
        if xpoint == 0.:
            count += 1
        else:
            break
    if count > 1:
        dx = xpoints[count+1] - xpoints[count]
        xpoints[count-1::-1] = [ xpoints[count] - jj*dx for jj in np.arange(count) ]
    return xpoints

def get_xpoints_bounds(xpoints):
    # get approximate bounds on the xpoints array by interpolating to midpoints and extrapolating at end points.
    xpoints_bounds = np.zeros(xpoints.shape[0]+1)
    print('xpoints.shape, xpoints_bounds.shape : ',xpoints.shape,xpoints_bounds.shape)
    xpoints_bounds[1:-1] = 0.5 * ( xpoints[:-1] + xpoints[1:] )
    xpoints_bounds[0] = xpoints[0] - (xpoints_bounds[1] - xpoints[0])
    xpoints_bounds[-1] = xpoints[-1] + (xpoints[-1] - xpoints_bounds[-2])
    return xpoints_bounds

def get_depth_bounds(depths):
    depth_bounds = np.zeros(len(depths)+1)
    depths_next = np.roll(depths,-1)
    half_thickness = depths[0]
    for count,(depth,depthnext) in enumerate(zip(depths,depths_next)):
        depth_bounds[count+1] = depth_bounds[count] + 2*half_thickness
        half_thickness = depthnext - depth - half_thickness
    return depth_bounds

def fill_field(fld_in):
    # fill out missing data values from neighbouring values by one point in each direction
    # otherwise the contouring stops before it gets to the mask. 
    
    fld_out = fld_in.copy()

    # fill to the left
    fld_shift = np.roll(fld_in,shift=-1,axis=0)
    fld_out[:-1,:] = ma.where( ( fld_out[:-1,:].mask & ~fld_shift[:-1,:].mask ), fld_shift[:-1,:], fld_out[:-1,:] )    

    # fill to the left
    fld_shift = np.roll(fld_in,shift=+1,axis=0)
    fld_out[1:,:] = ma.where( ( fld_out[1:,:].mask & ~fld_shift[1:,:].mask ), fld_shift[1:,:], fld_out[1:,:] )    

    # fill above (we might be in a cavity)
    fld_shift = np.roll(fld_in,shift=-1,axis=1)
    fld_out[:,:-1] = ma.where( ( fld_out[:,:-1].mask & ~fld_shift[:,:-1].mask ), fld_shift[:,:-1], fld_out[:,:-1] )    

    # fill below
    fld_shift = np.roll(fld_in,shift=+1,axis=1)
    fld_out[:,1:] = ma.where( ( fld_out[:,1:].mask & ~fld_shift[:,1:].mask ), fld_shift[:,1:], fld_out[:,1:] )    

    return fld_out

def find_section_indices(slon,slat,nav_lat,nav_lon):
    # routine adapted from code by Pat Hyder
    slon[slon<0.]=slon[slon<0.]+360.  
    latm=np.squeeze(nav_lat[:,1])
    indm=np.where(abs(latm-np.mean(slat))==min(abs((latm-np.mean(slat)))))
    indm=indm[0]
    print(indm)
    print(latm)
    dist_int=abs(latm[indm-1]-latm[indm])*60*1.852
    dist_int=10.0
    
    # set up the transect at the right distance internal

    print('size slon', size(slon))
    sdist=ones(size(slon))*nan
    sdist[0]=0.0
    # calculate a distance coordinate along the zig zag section
    for lp in range(1,size(slon)):
        sdistx=(slon[lp]-slon[lp-1])*1.852*60.*math.cos(pi/180.*(slat[lp]+slat[lp-1])/2)
        sdisty=(slat[lp]-slat[lp-1])*1.852*60.
        sdist[lp]=sdist[lp-1]+(sdistx**2+sdisty**2)**0.5
  
    # find the indices of the points at model res along the zig zag transect
    npt=int(max(sdist)/dist_int)+1
    print('sdist, dist_int, npt : ',sdist, dist_int, npt)
    sxdist=np.ones(npt)*np.nan
    sxlon=np.ones(npt)*np.nan
    sxlat=np.ones(npt)*[np.nan]
    sxdist[0]=0.0
    sxlon[0]=slon[0]
    sxlat[0]=slat[0]
    for lp in range(1,npt):
        print('Working on ',lp,' of ',npt,' points')
        sxdist_tmp=dist_int*float(lp)
        sind=np.nanmax(where(sdist-sxdist_tmp<=0))
        eind=np.nanmin(where(sdist-sxdist_tmp>0))
        dist=sxdist_tmp-sdist[sind]
        tdist=sdist[eind]-sdist[sind]
        sxlat[lp]=slat[sind]+dist/tdist*(slat[eind]-slat[sind])
        sxlon[lp]=slon[sind]+dist/tdist*(slon[eind]-slon[sind])
        sxdistx=(sxlon[lp]-sxlon[lp-1])*1.852*60.*math.cos(pi/180.*(sxlat[lp]+sxlat[lp-1])/2)
        sxdisty=(sxlat[lp]-sxlat[lp-1])*1.852*60.
        sxdist[lp]=sxdist[lp-1]+(sxdistx**2+sxdisty**2)**0.5
        #print(sxlat[lp],sxlon[lp],sxdist[lp])

    nav_lonm=nav_lon.copy()
    nav_lonm[nav_lonm<0]=nav_lonm[nav_lonm<0]+360.
    xind=npt*[np.nan]
    yind=npt*[np.nan]
    # find the x and y model indices for the transect
    for lp in range(0,npt):
        print('Working on ',lp,' of ',npt,' points')
        adistx=(nav_lonm-sxlon[lp])*1.852*60.*math.cos(pi/180.*sxlat[lp])
        adisty=(nav_lat-sxlat[lp])*1.852*60.
        adist=(adistx**2+adisty**2)**0.5
        [yindv,xindv]=where(adist == np.nanmin(adist))
        
        yind[lp]=yindv[0]
        xind[lp]=xindv[0]
        #print(xind,yind)
        #raw_input('here')

    return(xind,yind,sxlat,sxlon,sxdist)

def get_section(lons,lats,lons_bounds=None,lats_bounds=None,toldeg=None,method=None,
                lon=None,lat=None,index_i=None,index_j=None,xmin=None,xmax=None,fill_south=None):

    # Routine to return lists of indices and coordinates along a section using one of several methods
  
    print('get_section : lon, lat, index_i, index_j, xmin, xmax : ',lon, lat, index_i, index_j, xmin, xmax)

    # try to cope with different ranges for input longitudes:
    centred_0 = False
    centred_180 = False
    if lon is not None:
        if lon < 0.0:
            centred_0 = True
        elif lon > 180.0:
            centred_180 = True
    elif lat is not None:
        if xmin < 0.0:
            centred_0 = True
        elif xmax > 180.0:
            centred_180 = True
    if centred_0:        
        select_mask = np.where(lons > 180.0,1,0)
        lons = lons - 360.0*select_mask
        if lons_bounds is not None:
            select_mask = np.where(lons_bounds > 180.0,1,0)
            lons_bounds = lons_bounds - 360.0*select_mask
    if centred_180:        
        select_mask = np.where(lons < 0.0,1,0)
        lons = lons + 360.0*select_mask
        if lons_bounds is not None:
            select_mask = np.where(lons_bounds < 0.0,1,0)
            lons_bounds = lons_bounds + 360.0*select_mask

    if index_i is not None:
        xpoints_name='latitude (degrees north)'
        # interpret xmin/xmax as min/max j-index:
        if xmin is None:
            xmin = 0
        if xmax is None:
            xmax = lats.shape[0]-1
        jlist = np.int_(np.arange(xmax-xmin+1)+xmin)
        ilist = np.int_(np.ones((len(jlist)))*index_i)
        xsec_indices=(jlist,ilist)
        if len(lats.shape) == 1:
            xpoints = lats[xsec_indices[0]]
            if lats_bounds is not None:
                xpoints_bounds = np.concatenate((lats_bounds[xsec_indices[0],0],np.array([lats_bounds[xmax+1,1]])))
        else:
            xpoints = lats[xsec_indices]
            if lats_bounds is not None:
                xpoints_bounds = np.concatenate((lats_bounds[xsec_indices[0],xsec_indices[1],0],np.array([lats_bounds[xmax+1,index_i,1]])))
        if lats_bounds is None:
            xpoints_bounds = get_xpoints_bounds(xpoints)
        if fill_south:
            xpoints = fill_in_south(xpoints)
            xpoints_bounds = get_xpoints_bounds(xpoints)

    elif index_j is not None:
        xpoints_name='longitude (degrees east)'
        # interpret xmin/xmax as min/max i-index:
        if xmin is None:
            xmin = 0
        if xmax is None:
            xmax = lons.shape[-1]-1
        ilist = np.int_(np.arange(xmax-xmin+1)+xmin)
        jlist = np.int_(np.ones((len(ilist)))*index_j)
        xsec_indices=(jlist,ilist)
        if len(lons.shape) == 1:
            xpoints = lons[ilist]
            if lons_bounds is not None:
                xpoints_bounds = np.concatenate((lons_bounds[xsec_indices[1],0],np.array([lons_bounds[xmax+1,1]])))
        else:
            xpoints = lons[xsec_indices]
            if lons_bounds is not None:
                xpoints_bounds = np.concatenate((lons_bounds[xsec_indices[0],xsec_indices[1],0],np.array([lons_bounds[index_j,xmax+1,1]])))
        if lons_bounds is None:
            xpoints_bounds = get_xpoints_bounds(xpoints)
    elif lat is not None:
        xpoints_name='longitude (degrees east)'
        if xmin is None:
            xmin = lons.min()
        if xmax is None:
            xmax = lons.max()
        print('lons.min(),lons.max()',lons.min(),lons.max())
        print('xmin,xmax',xmin,xmax)
        if method == 'gridline':
            # find nearest gridline to specified lat/lon:
            if len(lons.shape) == 1:
                if xsec_indices is None:
                    j0 = ma.argmin((lats[:]-lat)*(lats[:]-lat))
                    i0 = ma.argmin((lons[:]-xmin)*(lons[:]-xmin)) 
                    i1 = ma.argmin((lons[:]-xmax)*(lons[:]-xmax))                
                    print('j0,[i0:i1] : ',j0,',[',i0,',',i1,']')
                    ilist = np.arange(i1-i0+1)+i0
                    jlist = np.int_(np.ones((len(ilist)))*j0)
                    xsec_indices=(jlist,ilist)
                xpoints = lons[xsec_indices[1]]
                if lons_bounds is not None:
                    xpoints_bounds = np.concatenate((lons_bounds[xsec_indices[1],0],np.array([lons_bounds[i1+1,1]])))
                else:
                    xpoints_bounds = get_xpoints_bounds(xpoints)
            else:
                lon_mid = 0.5 * (xmin + xmax)
                (j0,i_mid) = np.unravel_index(ma.argmin((lons[:]-lon_mid)*(lons[:]-lon_mid)+(lats[:]-lat)*(lats[:]-lat)),lats.shape)
                i0 = ma.argmin((lons[j0,:]-xmin)*(lons[j0,:]-xmin)) 
                i1 = ma.argmin((lons[j0,:]-xmax)*(lons[j0,:]-xmax))                
                print('j0,[i0:i1] : ',j0,',[',i0,',',i1,']')
                ilist = np.arange(i1-i0+1)+i0
                jlist = np.int_(np.ones((len(ilist)))*j0)
                xsec_indices=(jlist,ilist)
                xpoints = lons[xsec_indices]
                if lons_bounds is not None:
                    xpoints_bounds = np.concatenate((lons_bounds[jlist,ilist,0],np.array([lons_bounds[j0,i1+1,1]])))
                else:
                    xpoints_bounds = get_xpoints_bounds(xpoints)
        elif method == 'dave':
            xsec_indices_unsorted = np.where( ((lat-toldeg < lats) & (lats < lat+toldeg)) & ((lons > xmin) & (lons < xmax)) ) 
            xpoints = lons[xsec_indices_unsorted]
            sort_indices = np.argsort(xpoints)
            xsec_indices = (xsec_indices_unsorted[0][sort_indices],xsec_indices_unsorted[1][sort_indices])
            xpoints = xpoints[sort_indices]
            xpoints_bounds = get_xpoints_bounds(xpoints)
        elif method == 'pat':
            start_lat = lat
            start_lon = xmin
            end_lat = lat
            end_lon = xmax
            lona = np.array((start_lon,end_lon))
            lata = np.array((start_lat,end_lat))
            [sec_xind,sec_yind,sec_lat,sec_lon,sec_dist] = find_section_indices(lona,lata,lats,lons)
#            fld_xsec_list.append(np.squeeze(fld.data[:,sec_yind,sec_xind]).copy())
            xsec_indices=(sec_yind,sec_xind)
            xpoints = np.squeeze(lons[sec_yind,sec_xind])
#            xpoints_bounds = ??
    elif lon is not None:
        xpoints_name='latitude (degrees north)'
        if xmin is None:
            xmin = lats.min()
        if xmax is None:
            xmax = lats.max()
        print('lats.min(),lats.max()',lats.min(),lats.max())
        print('xmin,xmax',xmin,xmax)
        if method == 'gridline':
            # find nearest gridline to specified lon:
            if len(lons.shape) == 1:
                i0 = ma.argmin((lons[:]-lon)*(lons[:]-lon))
                j0 = ma.argmin((lats[:]-xmin)*(lats[:]-xmin)) 
                j1 = ma.argmin((lats[:]-xmax)*(lats[:]-xmax))                
                print('i0,[j0:j1] : ',i0,',[',j0,',',j1,']')
                jlist = np.arange(j1-j0+1)+j0
                ilist = np.int_(np.ones((len(jlist)))*i0)
                xsec_indices=(jlist,ilist)
                xpoints = lats[jlist]
                if lats_bounds is not None:
                    xpoints_bounds = np.concatenate((lats_bounds[jlist,0],np.array([lons_bounds[j1+1,1]])))
                else:
                    xpoints_bounds = get_xpoints_bounds(xpoints)
            else:
                lat_mid = 0.5 * (xmin + xmax)
                (j_mid,i0) = np.unravel_index(ma.argmin((lats[:]-lat_mid)*(lats[:]-lat_mid)+(lons[:]-lon)*(lons[:]-lon)),lons.shape)
                j0 = ma.argmin((lats[:,i0]-xmin)*(lats[:,i0]-xmin)) 
                j1 = ma.argmin((lats[:,i0]-xmax)*(lats[:,i0]-xmax))                
                print('i0,[j0:j1] : ',i0,',[',j0,',',j1,']')
                jlist = np.arange(j1-j0+1)+j0
                ilist = np.int_(np.ones((len(jlist)))*i0)
                xsec_indices=(jlist,ilist)
                xpoints = lats[xsec_indices]
                if lats_bounds is not None:
                    xpoints_bounds = np.concatenate((lats_bounds[jlist,ilist,0],np.array([lats_bounds[j1+1,i0,1]])))
                else:
                    xpoints_bounds = get_xpoints_bounds(xpoints)
        elif method == 'dave':
            xsec_indices_unsorted = np.where( ((lon-toldeg < lons) & (lons < lon+toldeg)) & ((lats > xmin) & (lats < xmax)) ) 
            xpoints = lats[xsec_indices_unsorted]
            sort_indices = np.argsort(xpoints)
            xsec_indices = (xsec_indices_unsorted[0][sort_indices],xsec_indices_unsorted[1][sort_indices])
            xpoints = xpoints[sort_indices]
            xpoints_bounds = get_xpoints_bounds(xpoints)
#        elif method == 'pat':
#            start_lat = lat
#            start_lon = xmin
#            end_lat = lat
#            end_lon = xmax
#            lona = np.array((start_lon,end_lon))
#            lata = np.array((start_lat,end_lat))
#            [sec_xind,sec_yind,sec_lat,sec_lon,sec_dist] = find_section_indices(lona,lata,lats,lons)
#            xsec_indices=(sec_yind,sec_xind)
#            xpoints = np.squeeze(lons[sec_yind,sec_xind])

    return xsec_indices,xpoints,xpoints_name,xpoints_bounds

########################################## Main routine ################################################

def plot_nemo_section(filenames=None,var_names=None,title=None,
                 xsec_indices=None,xpoints=None,xpoints_bounds=None,xpoints_name=None,
                 lat=None,lon=None,index_i=None,index_j=None,levs=None,nlevs=15,
                 rec=None,mnfld=None,mxfld=None,xmin=None,xmax=None,depthmax=None,method=None,logscale=None, 
                 reverseX=None,kilometres=None,outfile=None,toldeg=None,
                 noshow=None,nobar=None,scientific=None,figx=None,figy=None,subplot=None,
                 plot_types=None,colors=None,cmap=None,factor=None,fmt=None,fontsizes=None,
                 maskzero=None,draw_points=None,draw_fmt=None,text=None,textsize=None,
                 xtrac=None,label_lon=None,label_lat=None,fill_south=None,units=None,vertbar=None):

    print('>>> rec is : ',rec)
    print('>>> method is : ',method)

    matplotlib.rcParams.update({'font.size': 7})

    if plot_types is None:
        plot_types=['c','l']

#############################
# 1. Check control parameters
#############################

    if not xtrac:
        if xsec_indices is not None:
            if xpoints is None or xpoints_bounds is None:
                raise Exception('Error: if you specify xsec_bounds in the input, must also specify xpoints and xpoints_bounds')
            else:
                print("Plotting section using pre-calculated indices and along-section coordinates.")
                print("NB. If two fields specified the *same* section will be plotted for each field.")
        elif lat is not None:
            if lon is not None or index_i is not None or index_j is not None:
                raise Exception("Only specify one of lat (-p), lon (-l), index_i (-I), index_j (-J).")
            else:
                print("Plotting section at latitude ",lat)
        elif lon is not None:
            if index_i is not None or index_j is not None:
                raise Exception("Only specify one of lat (-p), lon (-l), index_i (-I), index_j (-J).")
            else:
                print("Plotting section at longitude ",lon)
        elif index_i is not None:
            if index_j is not None:
                raise Exception("Only specify one of lat (-p), lon (-l), index_i (-I), index_j (-J).")
            else:
                print("Plotting section at i-index = ",index_i)
        elif index_j is not None:
            print("Plotting section at j-index = ",index_j)
        else:
            raise Exception("Error: must specify one of lat (-P), lon (-L), index_i (-I), index_j (-J)")
        
    if var_names is None:
        raise Exception("Error: must specify at least one field to plot.")
    elif type(var_names) is not list:
        var_names=[var_names]
    elif len(var_names) > 2:
        raise Exception("Error: can specify at most two fields to plot.")
    nvar = len(var_names)

    if rec is None:
        rec = 0
    if type(rec) is not list:
        rec = nvar*[rec]
    elif len(rec) != nvar:
        raise Exception("Error: number of values for REC (-r) doesn't match number of fields (-v).")

    if method is None:
        method = 'dave'

    if nlevs is None:
        nlevs = 15
    if type(nlevs) is not list:
        nlevs = nvar*[nlevs]
    elif len(nlevs) != nvar:
        raise Exception("Error: number of values for NLEVS (-n) doesn't match number of fields (-v).")

    # apply maskzero to all fields for now
    if maskzero is None:
        maskzero=False

    if type(factor) is not list:
        factor = nvar * [factor]
    elif len(factor) != nvar:
        raise Exception("Error: number of values for FACTOR (-Z) doesn't match number of fields (-v).")

    if filenames is None:
        raise Exception("Error: must specify at least one input file.")
    elif type(filenames) is not list:
        filenames = [filenames]
    if len(filenames) == 1 and len(var_names) == 2:
        # if two fields and only one input file try to read both fields from the same file:
        filenames=[filenames,filenames] 

    if xtrac and label_lon and label_lat:
        raise Exception("Error: should only specify one (or neither) of label_lon and label_lat")

    fontsizes_list = [14,12,10,8]
    for count,fontsize in enumerate(fontsizes):
        if fontsize is not None:
            fontsizes_list[count]=fontsize
    fontsizes = fontsizes_list

    csline=None
    cscolor=None

#############################
# 2. Read in field(s)
#############################

    fld_list = []
    for filename, varname, rec_i in zip(filenames,var_names,rec):
        fld_read = gt.read_cube(filename,varname)

# Load required record from file.
# NB. Assuming the field has 3 spatial dimensions here.
        if fld_read.ndim == 4:
            fld_list.append(fld_read[rec_i].copy())
        else:
            fld_list.append(fld_read.copy())

#########################################
# 3. Extract cross section if necessary
#########################################

    fld_xsec_list = []
    xpoints_list = []
    xpoints_bounds_list = []

    if xtrac:
        # in this case we have already read in the section directly from the input file
        for fld,filename in zip(fld_list,filenames):
            if label_lon:
                xpoints_cube = gt.read_cube(filename,'longitude')
                xpoints = xpoints_cube.data[0,:-1]
                xpoints_name='longitude (degrees east)'
            elif label_lat:
                xpoints_cube = gt.read_cube(filename,'latitude')
                xpoints = xpoints_cube.data[0,:-1]
                xpoints_name='latitude (degrees north)'
            else:
                # if we don't want lon or lat labels on x-axis then revert to
                # distance along section in kilometres.
                e1v_cube = gt.read_cube(filename,'e1v')
                e1v = e1v_cube.data[0,0,:-1]
                xpoints = np.cumsum(e1v,axis=-1) * 0.001
                xpoints_name='distance along section (km)'
            #print('xpoints : ',xpoints)

            # xtrac fields always have a masked point at the end - don't read this. 
            flddata = fld.data[:,:,:-1].squeeze()
            # truncate data if xmin and/or xmax set
            if xmin is not None:
                indices = np.where(xpoints >= xmin)[0]
                xpoints = xpoints[indices]
                flddata = flddata[:,indices]
            if xmax is not None:
                indices = np.where(xpoints <= xmax)[0]
                xpoints = xpoints[indices]
                flddata = flddata[:,indices]

            # no "bounds" variables in xtrac files, so just use get_xpoints_bounds:
            xpoints_bounds = get_xpoints_bounds(xpoints)

            # if xpoints is distance in km, rescale to start at zero at start of section
            if not label_lon and not label_lat:
                xpoints[:] = xpoints[:] - xpoints_bounds[0]
                xpoints_bounds[:] = xpoints_bounds[:] - xpoints_bounds[0]

            xpoints_list.append(xpoints)
            xpoints_bounds_list.append(xpoints_bounds)
            fld_xsec_list.append(flddata)
            print('fld_xsec_list[-1].shape',fld_xsec_list[-1].shape)
    
    else:
        # if the input files have 3D fields then extract the required section
        if toldeg is None:
            toldeg = 0.1

        for fld,filename in zip(fld_list,filenames):

            if xsec_indices is None:
                lons, lons_bounds = get_coord(fld,'longitude',filename)
                lats, lats_bounds = get_coord(fld,'latitude',filename)

                if len(lons.shape) == 1:
                    if len(lats.shape) != 1:
                        print('len(lons.shape), len(lats.shape) : ',len(lons.shape),len(lats.shape))
                        raise Exception('Error - 1D longitude coordinate and > 1D latitude coordinate. Something is weird here.')
                    else:
                        print('1D lat and lon coordinates in file => defaulting to gridline method.')
                        method='gridline'

                (xsec_indices, xpoints, xpoints_name, xpoints_bounds) = get_section(lons,lats,
                                                   lons_bounds=lons_bounds,lats_bounds=lats_bounds, 
                                                   toldeg=toldeg,method=method,lon=lon,lat=lat,
                                                   index_i=index_i,index_j=index_j,xmin=xmin,xmax=xmax,
                                                   fill_south=fill_south)

#            print('xpoints, flddta : ')
#            for xpoint,flddta in zip(xpoints,fld.data[27,xsec_indices[0],xsec_indices[1]]):
#                print(xpoint,':',flddta)
            fld_xsec_list.append(fld.data[:,xsec_indices[0],xsec_indices[1]])
            print('fld_xsec_list[-1].shape : ',fld_xsec_list[-1].shape)
            if maskzero:
                fld_xsec_list[-1] = ma.masked_values(fld_xsec_list[-1],0.0)
            xpoints_list.append(xpoints.copy())
            xpoints_bounds_list.append(xpoints_bounds.copy())

#############################
# 4. Do the plot
#############################

#    =========  SANITY CHECK FOR PLOTTING ========
#    print('fld_xsec.shape : ',fld_xsec.shape)
#    nz = fld_xsec.shape[0]
#    nz3 = nz/3
#    nz23 = 2*nz/3
#    print('nz, nz3, nz23 : ',nz, nz3, nz23)
#    print('fld_xsec[:nz3,:].shape : ',fld_xsec[:nz3,:].shape)
#    fld_xsec[:nz3,:] = 3.0e-5
#    fld_xsec[nz3:nz23,:] = 0.0
#    fld_xsec[nz23:,:] = -3.0e-5
#    print('ma.min(fld_xsec), ma.max(fld_xsec) : ',ma.min(fld_xsec), ma.max(fld_xsec))
#    =============================================

    # This line disables "summary printing" of numpy arrays, ie. forces it to print the full array.
#    np.set_printoptions(threshold=np.nan)

    if type(mnfld) is not list:
        mnfld = nvar*[mnfld]
    elif len(mnfld) != nvar:
        raise Exception("Error: number of values for MNFLD (-f) doesn't match number of fields.")
    if type(mxfld) is not list:
        mxfld = nvar*[mxfld]
    elif len(mxfld) != nvar:
        raise Exception("Error: number of values for MXFLD (-F) doesn't match number of fields.")

    if colors is None and nvar == 2:
        # for a line contour plot over a filled contour plot default to black lines
        colors = "black"

    if colors is not None and nvar == 1:
        # colors takes precedence over cmap and matplotlib.pyplot.contour complains
        # if you give it colors and cmap both as not None.
        cmap = None

    # sort out colour map for filled or line contouring
    if cmap is not None:
        # assume that it refers to a matplotlib inbuilt colour map:
        print('reading inbuilt matplotlib colour map : ',cmap)
        cmap = getattr(matplotlib.cm,cmap)

    if colors is None and cmap is None:
        cmap = getattr(matplotlib.cm,'RdYlBu_r')

    # loop over fields to plot (2 max)
    count = 0
    for var_name,plot_type,xpoints,xpoints_bounds,fld,fld_xsec,fac,mnfld_i,mxfld_i,nlevs_i in \
       zip(var_names,plot_types[:nvar],xpoints_list,xpoints_bounds_list,fld_list,fld_xsec_list,factor,mnfld,mxfld,nlevs):

        #print("xpoints : ",xpoints)
        #print("xpoints_bounds : ",xpoints_bounds)

        count += 1

        if fac is not None:
            print('multiplying field ',var_name,' by ',fac)
            fld_xsec[:] = fac*fld_xsec[:]

        if mnfld_i is None or mnfld_i == "None":
             mnfld_i = fld_xsec.min()
        else:
             mnfld_i = float(mnfld_i)
        if mxfld_i is None or mxfld_i == "None":
            mxfld_i = fld_xsec.max()
        else:
            mxfld_i = float(mxfld_i)

        print('mnfld_i, mxfld_i : ',mnfld_i,mxfld_i)

        print('nlevs_i : ',nlevs_i)
        if plot_type == 'c' or levs is None:
            if logscale:
                tiny_factor=0.0001
                if mnfld_i == 0.0:
                    mnfld_i = tiny_factor * mxfld_i
                if mxfld_i == 0.0:
                    mxfld_i = tiny_factor * mnfld_i
                if mnfld_i < 0.0 and mxfld_i < 0.0:
                    levs_i = -np.logspace(np.log10(-mnfld_i),np.log10(-mxfld_i),nlevs_i+1) 
                if mnfld_i > 0.0 and mxfld_i > 0.0:
                    levs_i = np.logspace(np.log10(mnfld_i),np.log10(mxfld_i),nlevs_i+1) 
                if mnfld_i < 0.0 and mxfld_i > 0.0:
                    nlevs_neg = int(abs(mnfld_i)*nlevs_i/(mxfld_i-mnfld_i))
                    levs_neg = -np.logspace(np.log10(-mnfld_i),np.log10(-mnfld_i*tiny_factor),nlevs_i+1) 
                    nlevs_pos = int(mxfld_i*nlevs_i/(mxfld_i-mnfld_i))
                    levs_pos = np.logspace(np.log10(mxfld_i*tiny_factor),np.log10(mxfld_i),nlevs_i+1) 
                    levs_i = levs_neg + levs_pos
            else:
                levs_i = np.linspace(mnfld_i,mxfld_i,nlevs_i+1) 
        else:
            levs_i = levs
        print('levs_i : ',levs_i)

        for depthname in ['depth', 'deptht', 'depthu', 'depthv', 'depthw', 'lev']:
            try:
                depth = fld.coord( var_name = depthname ) 
            except iris.exceptions.CoordinateNotFoundError:
                pass
            else:
                break
        else:
            raise Exception("Error: could not find depth coordinate.")

        try:
            depth_bounds = np.concatenate((depth.bounds[:,0],np.array([depth.bounds[-1,1]])))
        except(TypeError):
            depth_bounds = get_depth_bounds(depth.points)
        print('depth_bounds : ',depth_bounds)

        if figx is None:
            figx = matplotlib.rcParams["figure.figsize"][0]
        if figy is None:
            figy = matplotlib.rcParams["figure.figsize"][1]
        print('figx,figy : ',figx,figy)

        # If figure 1 already exists then plt.figure will 
        # just return a reference to it.
        fig = plt.figure(1, figsize=[figx,figy])

        # Create plot axes:

        if subplot is not None:
            if len(subplot) == 1:
                ax = fig.add_subplot(subplot)
            elif len(subplot) == 3:
                ax = fig.add_subplot(subplot[0],subplot[1],subplot[2])
            else:
                raise Exception("Error : len(subplot) must be 1 or 3.")
        else:
            ax = plt.gca()

        if count == 1: 
            # Plot the mask - only for the first field or we overwrite previously-plotted fields!
            # Keep the original mask before we fill fld_xsec prior to contouring.
            mask_to_plot = fld_xsec.mask.copy()
            # Make a colormap with grey on the end instead of black
            # If no masked values on the section then mask=False which the plotting routine
            # doesn't like so check that mask_to_plot is iterable first.
            try:
                iterator = iter(mask_to_plot)
            except TypeError:
                pass
            else:
                cmap_mask = matplotlib.colors.ListedColormap(cm.binary(np.arange(128)))
                cscolor_mask = plt.pcolormesh(xpoints_bounds, depth_bounds, mask_to_plot, cmap=cmap_mask)
 

        if plot_type == 'l':
            ### contour lines ###
            fld_xsec_filled = fill_field(fld_xsec)
            print("Plotting line contours.")
            if colors is not None and cmap is not None:
                csline = plt.contour(xpoints,depth.points,fld_xsec_filled, colors=colors, levels=levs_i,linewidths=0.5)
            else:
                csline = plt.contour(xpoints,depth.points,fld_xsec_filled, colors=colors, cmap=cmap, 
                                     levels=levs_i,linewidths=0.5)
#            manual_locations=[(35,1000)]
            if fmt is None:
                fmt='%1.1f'
            plt.clabel(csline,csline.levels[::2],inline=True,manual=False,fontsize=fontsizes[3],fmt=fmt)
        elif plot_type == 'c':
            ### colour-filled contours ###
            print("Plotting colour-filled contours.")
            fld_xsec_filled = fill_field(fld_xsec)
            cscolor = plt.contourf(xpoints,depth.points,fld_xsec,cmap=cmap,levels=levs_i)
        elif plot_type == 'b':
            ### colour-filled blocks ###
            if colors is not None and nvar == 1:
                (cmap_block,norm_block) = matplotlib.colors.from_levels_and_colors(levs_i,colors)
            else:
                cmap_block = cmap
            cscolor = plt.pcolormesh(xpoints_bounds, depth_bounds, fld_xsec, cmap=cmap_block, vmin=mnfld_i, vmax=mxfld_i)
        elif plot_type == 'bl':
            ### block plot and lines for the same field ###
            if colors is not None and nvar == 1:
                (cmap_block,norm_block) = matplotlib.colors.from_levels_and_colors(levs_i,colors)
            else:
                cmap_block = cmap
            cscolor = plt.pcolormesh(xpoints_bounds, depth_bounds, fld_xsec, cmap=cmap_block, vmin=mnfld_i, vmax=mxfld_i)
            fld_xsec_filled = fill_field(fld_xsec)
            csline = plt.contour(xpoints,depth.points,fld_xsec_filled, colors="black", levels=levs_i,linewidths=0.5)
            if fmt is None:
                fmt='%1.1f'
            plt.clabel(csline,inline=True,manual=False,fontsize=7,fmt=fmt)
        else:
            raise Exception("Error: unknown plot_type: "+plot_type)

        if not nobar and plot_type != 'l':
            if vertbar:
                cax = plt.colorbar(cscolor,orientation='vertical')
            else:
                cax = plt.colorbar(cscolor,orientation='horizontal')
            if plot_type == 'b':
#                print('norm_block(levs) : ',norm_block(levs))
#                cax.set_ticks(norm_block(levs))
                print('levs_i : ',levs_i)
                cax.set_ticks(levs_i)
            else:
                print('levs_i : ',levs_i)
                cax.set_ticks(levs_i)
            if scientific:
                labels=["%1.2e" %lev for lev in levs_i]
            else:
                labels=["%1.2f" %lev for lev in levs_i]
            if vertbar:
                cax.ax.set_yticklabels(labels,fontsize=fontsizes[2])
            else:
                cax.ax.set_xticklabels(labels,rotation=45,fontsize=fontsizes[2])
            if len(labels) > 11:
                if vertbar:
                    for label in cax.ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)           
                else:
                    for label in cax.ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)           
            if units is not None:
                cax.ax.set_xlabel(str(units),fontsize='x-large')
            elif str(fld.units) != "1" and str(fld.units) != "unknown":
                cax.ax.set_xlabel(str(fld.units),fontsize='x-large')
    
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
            plot_depth = [draw_points[ii+1],draw_points[ii+3]]
            plt.plot(plot_x[:],plot_depth[:],fmt,linewidth=linewidth)

    if title is not None:
        ax.set_title(textwrap.fill(title,60),fontsize=fontsizes[0])
    if text is not None:
        if type(text) is list:
            if len(text)%3 !=0:
                raise Exception('Error: must specify sets of 3 arguments for text: x, y, s. As for matplotlib text.')
        else:
            raise Exception('Error: text should consist of multiples of 3 arguments.')
        if textsize is None:
            textsize=12
        for ii in range(0,len(text),3):
            ax.text(text[ii],text[ii+1],text[ii+2],fontsize=textsize)
#    if xtrac and (xmin is not None):
#        ax.set_xlim(xmin,ax.get_xlim()[-1])
#    if xtrac and (xmax is not None):
#        ax.set_xlim(ax.get_xlim()[0],xmax)
    if reverseX:
        ax.set_xlim(ax.get_xlim()[::-1])
    if depthmax is not None:
        ax.set_ylim([depthmax,0.0])
    else:
        ax.invert_yaxis()

    ax.set_ylabel('depth (m)',fontsize=fontsizes[1])
    if kilometres:
        if lon is not None or 'lat' in xpoints_name:
            ax.set_xlabel('south-north distance (km)',fontsize=fontsizes[1])
        elif lat is not None or 'lon' in xpoints_name:
            ax.set_xlabel('west-east distance (km)',fontsize=fontsizes[1])
    else:
        if xpoints_name is not None:
            ax.set_xlabel(xpoints_name,fontsize=fontsizes[1])

    ax.xaxis.set_tick_params(labelsize=fontsizes[2])
    ax.yaxis.set_tick_params(labelsize=fontsizes[2])

    if outfile is not None:
        plt.draw()
        plt.savefig(outfile,dpi=200)
    elif not noshow:
        plt.show()

    # Return references to the plots for use in calling python routines. 
    # In general one or more of these might have None values. 
    return csline, cscolor

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file",action="store",dest="filenames", nargs="+",
                    help="name(s) of input file(s) for scalar field 1 and scalar field 2")
    parser.add_argument("--xtrac", action="store_true",dest="xtrac",
                    help="input file(s) are section files created by cdf_xtrac_brokenline")
    parser.add_argument("-v", "--var",action="store",dest="var_names", nargs="+",
                    help="name(s) of field(s) to plot.") 
    parser.add_argument("-Z", "--factor", action="store",dest="factor",type=float,nargs="+",
                    help="multiplicative factor(s) to apply to plotted field(s)")
    parser.add_argument("-z", "--maskzero", action="store_true",dest="maskzero",
                    help="mask zero values in field")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",
                    help="plot to output file - filetype by extension. If unset plot to GUI.")
    parser.add_argument("-p", "--plot_types", action="store",dest="plot_types",nargs="+",
                    help="type of plot for each input field: 'l' = contour lines, 'c' = colour-filled contours, 'b' = block plot")
    parser.add_argument("-P", "--phi", action="store",dest="lat",type=float,
                    help="latitude of cross section")
    parser.add_argument("-L", "--lambda", action="store",dest="lon",type=float,
                    help="longitude of cross section")
    parser.add_argument("-I", "--index_i", action="store",dest="index_i",type=int,
                    help="i-index of cross section")
    parser.add_argument("--fill_south", action="store_true",dest="fill_south",
                    help="fill in zero latitude values for extended ORCA grids")
    parser.add_argument("-J", "--index_j", action="store",dest="index_j",type=int,
                    help="j-index of cross section")
    parser.add_argument("-m", "--levs", action="store",dest="levs",type=float,nargs='+',
                    help="specified contour levels for line plots")
    parser.add_argument("-n", "--nlevs", action="store",dest="nlevs",type=int,nargs='+',
                    help="number of contour levels for each field")
    parser.add_argument("-q", "--fmt", action="store",dest="fmt",
                    help="format string for contour labels (for line contours)")
    parser.add_argument("-r", "--rec", action="store",dest="rec",type=int,nargs='+',
                    help="record number to plot for each field (defaults to zero)")
    parser.add_argument("-x", "--xmin", action="store",dest="xmin",type=float,
                    help="start of x-axis (min lat or min lon as applicable)")
    parser.add_argument("-X", "--xmax", action="store",dest="xmax",type=float,
                    help="end of x-axis (max lat or max lon as applicable)")
    parser.add_argument("-R", "--reverseX", action="store_true",dest="reverseX",
                    help="reverse the x axis")
    parser.add_argument("-f", "--mnfld", action="store",dest="mnfld",nargs="+",
                    help="minimum field value to plot for colour filled contouring")
    parser.add_argument("-F", "--mxfld", action="store",dest="mxfld",nargs="+",
                    help="maximum field value to plot for colour filled contouring")
    parser.add_argument("-D", "--depthmax", action="store",dest="depthmax",type=float,
                    help="maximum depth for plot")
    parser.add_argument("-t", "--title", action="store",dest="title",
                    help="title for plot")
    parser.add_argument("-A", "--fontsizes", action="store",dest="fontsizes",nargs='+',
                    help="fontsizes for : plot title, axis labels, axis tick labels, contour labels. Use None for default.")
    parser.add_argument("-M", "--method", action="store",dest="method",
                    help="method for extracting section ('dave' or 'pat' or 'gridline')")
    parser.add_argument("-G", "--logscale", action="store_true",dest="logscale",
                    help="use colour logscale")
    parser.add_argument("-C", "--colors", action="store",dest="colors",nargs='+',
                    help="list of colors to use - takes precedence over colour map.")
    parser.add_argument("-c", "--cmap", action="store",dest="cmap",
                    help="colour map for contour plotting - inbuilt matplotlib colour map or external loaded with monty.clr_cmap")
    parser.add_argument("-N", "--noshow", action="store_true",dest="noshow",
                    help="Do not show GUI window even if output filename unset.")
    parser.add_argument("-b", "--nobar", action="store_true",dest="nobar",
                    help="Do not plot colour bar.")
    parser.add_argument("--label_lon", action="store_true",dest="label_lon",
                    help="x-axis to be labelled with longitudes (with --xtrac option)")
    parser.add_argument("--label_lat", action="store_true",dest="label_lat",
                    help="x-axis to be labelled with latitudes (with --xtrac option)")
    parser.add_argument("-k", "--kilometres", action="store_true",dest="kilometres",
                    help="x-axis to be labelled as kilometres rather than degrees (for f-plane models etc)")
    parser.add_argument("-Q", "--toldeg", action="store",dest="toldeg",type=float,
                    help="tolerance in degrees when calculating cross section (default=0.1)")
    parser.add_argument("-e", "--scientific", action="store_true",dest="scientific",
                    help="scientific format for colour bar labels")
    parser.add_argument("-d", "--draw_points", action="store",dest="draw_points",nargs="+",
                    help="list of points to draw line segments between in groups of four: start_x,start_depth,end_x,end_depth")
    parser.add_argument("-w", "--draw_fmt", action="store",dest="draw_fmt",nargs="+",
                    help="first argument is format for line segments as for fmt keyword for pyplot.plot; second optional argument is line thickness")
    parser.add_argument("-T", "--text", action="store",dest="text",nargs='+',
                    help="text to add to plot - multiples of 3 arguments: x, y, s as for matplotlib text")
    parser.add_argument("-S", "--textsize", action="store",dest="textsize",
                    help="fontsize of added text")
    parser.add_argument("-u", "--units", action="store",dest="units",default=None,
                    help="units (label for colour bar)")
    parser.add_argument("-V", "--vertbar", action="store_true",dest="vertbar",
                    help="use vertical color bar (default horizontal)")
    parser.add_argument("--figx", action="store",dest="figx",default=None,type=float,
                    help="x-dimension of figure (in inches I think)")
    parser.add_argument("--figy", action="store",dest="figy",default=None,type=float,
                    help="y-dimension of figure (in inches I think)")

    args = parser.parse_args()

    plot_nemo_section(filenames=args.filenames,var_names=args.var_names, outfile=args.outfile, title=args.title, fontsizes=args.fontsizes,
    lat=args.lat, lon=args.lon, index_i=args.index_i, index_j=args.index_j,levs=args.levs, nlevs=args.nlevs, rec=args.rec, 
    xmin=args.xmin, xmax=args.xmax,mnfld=args.mnfld, mxfld=args.mxfld, depthmax=args.depthmax, method=args.method, 
    logscale=args.logscale, reverseX=args.reverseX, kilometres=args.kilometres, toldeg=args.toldeg, 
    noshow=args.noshow, nobar=args.nobar, scientific=args.scientific, plot_types=args.plot_types,
    colors=args.colors, cmap=args.cmap, factor=args.factor,fmt=args.fmt,maskzero=args.maskzero,draw_points=args.draw_points,
    draw_fmt=args.draw_fmt,text=args.text,textsize=args.textsize,xtrac=args.xtrac,
    label_lon=args.label_lon,label_lat=args.label_lat,fill_south=args.fill_south,
    units=args.units,vertbar=args.vertbar,figx=args.figx,figy=args.figy)
       

