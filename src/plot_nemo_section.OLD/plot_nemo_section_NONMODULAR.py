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
            6. More compact version of the "dave" method which does exactly the same thing. 

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
import monty
import numpy as np
import numpy.ma as ma
import general_tools as gt
from pylab import *
import textwrap
# import glob

def get_xpoints_bounds(xpoints):
    # get approximate bounds on the xpoints array by interpolating to midpoints and extrapolating at end points.
    xpoints_bounds = np.zeros(xpoints.shape[0]+1)
    print 'xpoints.shape, xpoints_bounds.shape : ',xpoints.shape,xpoints_bounds.shape
    xpoints_bounds[1:-1] = 0.5 * ( xpoints[:-1] + xpoints[1:] )
    xpoints_bounds[0] = xpoints[0] - (xpoints_bounds[1] - xpoints[0])
    xpoints_bounds[-1] = xpoints[-1] + (xpoints[-1] - xpoints_bounds[-2])
    return xpoints_bounds

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
    print indm
    print latm
    dist_int=abs(latm[indm-1]-latm[indm])*60*1.852
    dist_int=10.0
    
    # set up the transect at the right distance internal

    print 'size slon', size(slon)
    sdist=ones(size(slon))*nan
    sdist[0]=0.0
    # calculate a distance coordinate along the zig zag section
    for lp in range(1,size(slon)):
        sdistx=(slon[lp]-slon[lp-1])*1.852*60.*math.cos(pi/180.*(slat[lp]+slat[lp-1])/2)
        sdisty=(slat[lp]-slat[lp-1])*1.852*60.
        sdist[lp]=sdist[lp-1]+(sdistx**2+sdisty**2)**0.5
  
    # find the indices of the points at model res along the zig zag transect
    npt=int(max(sdist)/dist_int)+1
    print 'sdist, dist_int, npt : ',sdist, dist_int, npt    
    sxdist=np.ones(npt)*np.nan
    sxlon=np.ones(npt)*np.nan
    sxlat=np.ones(npt)*[np.nan]
    sxdist[0]=0.0
    sxlon[0]=slon[0]
    sxlat[0]=slat[0]
    for lp in range(1,npt):
        print 'Working on ',lp,' of ',npt,' points'
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
        #print sxlat[lp],sxlon[lp],sxdist[lp]

    nav_lonm=nav_lon.copy()
    nav_lonm[nav_lonm<0]=nav_lonm[nav_lonm<0]+360.
    xind=npt*[np.nan]
    yind=npt*[np.nan]
    # find the x and y model indices for the transect
    for lp in range(0,npt):
        print 'Working on ',lp,' of ',npt,' points'
        adistx=(nav_lonm-sxlon[lp])*1.852*60.*math.cos(pi/180.*sxlat[lp])
        adisty=(nav_lat-sxlat[lp])*1.852*60.
        adist=(adistx**2+adisty**2)**0.5
        [yindv,xindv]=where(adist == np.nanmin(adist))
        
        yind[lp]=yindv[0]
        xind[lp]=xindv[0]
        #print xind,yind
        #raw_input('here')

    return(xind,yind,sxlat,sxlon,sxdist)


def plot_nemo_section(filenames=None,var_names=None,title=None,
                 lat=None,lon=None,index_i=None,index_j=None,levs=None,nlevs=15,
                 rec=0,mnfld=None,mxfld=None,
                 xmin=None,xmax=None,depthmax=None,method=None,logscale=None, 
                 reverseX=None,kilometres=None,outfile=None,toldeg=None,
                 noshow=None,nobar=None,scientific=None,block=None,
                 lineplot=None,colors=None,cmap=None,factor=None,fmt=None,
                 maskzero=None):

    print '>>> rec is : ',rec 
    print '>>> method is : ',method

    if lineplot:
        plot_types=['l','c']
    else:
        plot_types=['c','l']

#############################
# 1. Check control parameters
#############################

    if lat is not None:
        if lon is not None or index_i is not None or index_j is not None:
            raise Exception("Only specify one of lat (-p), lon (-l), index_i (-I), index_j (-J).")
        else:
            print "Plotting section at latitude ",lat
    elif lon is not None:
        if index_i is not None or index_j is not None:
            raise Exception("Only specify one of lat (-p), lon (-l), index_i (-I), index_j (-J).")
        else:
            print "Plotting section at longitude ",lon
    elif index_i is not None:
        if index_j is not None:
            raise Exception("Only specify one of lat (-p), lon (-l), index_i (-I), index_j (-J).")
        else:
            print "Plotting section at i-index = ",index_i
    elif index_j is not None:
        print "Plotting section at j-index = ",index_j
    else:
        raise Exception("Error: must specify one of lat (-p), lon (-l), index_i (-I), index_j (-J)")
        
    if var_names is None:
        raise Exception("Error: must specify at least one field to plot.")
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
    elif len(filenames) == 1 and len(var_names) == 2:
        # if two fields and only one input file try to read both fields from the same file:
        filenames=[filenames,filenames] 

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

#############################
# 2. Extract cross section
#############################

    if toldeg is None:
        toldeg = 0.1

    fld_xsec_list = []
    xpoints_list = []
    for fld,filename in zip(fld_list,filenames):

#   Lats/lons will be recognised by Iris as auxillary coordinates if the "coordinates"
#   attribute of the temperature field is set correctly. If not, the lat/lon fields
#   might still be in the file, so try reading them as simple fields rather than coordinates.
        try:
            lats=fld.coord('latitude').points
            lats_bounds=fld.coord('latitude').bounds
        except iris.exceptions.CoordinateNotFoundError:
            print 'reading latitudes as field'
            lats=iris.load_cube(filename,'latitude').data
        else:
            print 'reading latitudes as auxillary coordinate'

        try:
            lons=fld.coord('longitude').points
            lons_bounds=fld.coord('longitude').bounds
        except iris.exceptions.CoordinateNotFoundError:
            print 'reading longitudes as field'
            lons=iris.load_cube(filename,'longitude').data
        else:
            print 'reading longitudes as auxillary coordinate'

        select_mask = np.where(lons < 0.0,1,0)
        lons = lons + 360.0*select_mask
        select_mask = np.where(lons_bounds < 0.0,1,0)
        if lons_bounds is not None:
            lons_bounds = lons_bounds + 360.0*select_mask
            print 'lons_bounds.shape : ',lons_bounds.shape

        if len(lons.shape) == 1:
            if len(lats.shape) != 1:
                print 'len(lons.shape), len(lats.shape) : ',len(lons.shape),len(lats.shape)
                raise Exception('Error - 1D longitude coordinate and > 1D latitude coordinate. Something is weird here.')
            else:
                print '1D lat and lon coordinates in file => defaulting to gridline method.'
                method='gridline'

        if index_i is not None:
            # interpret xmin/xmax as min/max j-index:
            if xmin is None:
                xmin = 0
            if xmax is None:
                xmax = fld.data.shape[1]
            fld_xsec = fld.data[:,xmin:xmax+1,index_i]
            if len(lats.shape) == 1:
                xpoints = lats[xmin:xmax+1]
                xpoints_bounds = np.concatenate((lats_bounds[xmin:xmax+1,0],np.array([lats_bounds[xmax+1,1]])))
            else:
                xpoints = lats[xmin:xmax+1,index_i]
                xpoints_bounds = np.concatenate((lats_bounds[xmin:xmax+1,index_i,0],np.array([lats_bounds[xmax+1,index_i,1]])))
        elif index_j is not None:
            # interpret xmin/xmax as min/max i-index:
            if xmin is None:
                xmin = 0
            if xmax is None:
                xmax = fld.data.shape[2]
            fld_xsec = fld.data[:,index_j,xmin:xmax+1]
            if len(lons.shape) == 1:
                xpoints = lons[xmin:xmax+1]
                xpoints_bounds = np.concatenate((lons_bounds[xmin:xmax+1,0],np.array([lons_bounds[xmax+1,1]])))
            else:
                xpoints = lons[index_j,xmin:xmax+1]
                xpoints_bounds = np.concatenate((lons_bounds[index_j,xmin:xmax+1,0],np.array([lons_bounds[index_j,xmax+1,1]])))
        elif lat is not None:
            if xmin is None:
                xmin = lons.min()
            if xmax is None:
                xmax = lons.max()
            print 'lons.min(),lons.max()',lons.min(),lons.max()
            print 'xmin,xmax',xmin,xmax
            if method == 'gridline':
                # find nearest gridline to specified lat/lon:
                if len(lons.shape) == 1:
                    j0 = ma.argmin((lats[:]-lat)*(lats[:]-lat))
                    i0 = ma.argmin((lons[:]-xmin)*(lons[:]-xmin)) 
                    i1 = ma.argmin((lons[:]-xmax)*(lons[:]-xmax))                
                    print 'j0,[i0:i1] : ',j0,',[',i0,',',i1,']'
                    fld_xsec = fld.data[:,j0,i0:i1+1]
                    xpoints = lons[i0:i1+1]
                    xpoints_bounds = np.concatenate((lons_bounds[i0:i1+1,0],np.array([lons_bounds[i1+1,1]])))
                else:
                    lon_mid = 0.5 * (xmin + xmax)
                    (j0,i_mid) = np.unravel_index(ma.argmin((lons[:]-lon_mid)*(lons[:]-lon_mid)+(lats[:]-lat)*(lats[:]-lat)),lats.shape)
                    i0 = ma.argmin((lons[j0,:]-xmin)*(lons[j0,:]-xmin)) 
                    i1 = ma.argmin((lons[j0,:]-xmax)*(lons[j0,:]-xmax))                
                    print 'j0,[i0:i1] : ',j0,',[',i0,',',i1,']'
                    fld_xsec = fld.data[:,j0,i0:i1+1]
                    xpoints = lons[j0,i0:i1+1]
                    xpoints_bounds = np.concatenate((lons_bounds[j0,i0:i1+1,0],np.array([lons_bounds[j0,i1+1,1]])))
            elif method == 'dave':
                xsec_indices = np.where( ((lat-toldeg < lats) & (lats < lat+toldeg)) & ((lons > xmin) & (lons < xmax)) ) 
                print 'xsec_indices : ',zip(xsec_indices[1],xsec_indices[0])
                fld_xsec = fld.data[:,xsec_indices[0],xsec_indices[1]]
                xpoints = lons[xsec_indices]
                sort_indices = np.argsort(xpoints)
                xpoints = xpoints[sort_indices]
                xpoints_bounds = get_xpoints_bounds(xpoints)
                fld_xsec = fld_xsec[:,sort_indices]
            elif method == 'pat':
                start_lat = lat
                start_lon = xmin
                end_lat = lat
                end_lon = xmax
                lona = np.array((start_lon,end_lon))
                lata = np.array((start_lat,end_lat))
                [sec_xind,sec_yind,sec_lat,sec_lon,sec_dist] = find_section_indices(lona,lata,lats,lons)
                fld_xsec = np.squeeze(fld.data[:,sec_yind,sec_xind])
                xpoints = np.squeeze(lons[sec_yind,sec_xind])
                xpoints_bounds = get_xpoints_bounds(xpoints)
        elif lon is not None:
            if xmin is None:
                xmin = lats.min()
            if xmax is None:
                xmax = lats.max()
            print 'lats.min(),lats.max()',lats.min(),lats.max()
            print 'xmin,xmax',xmin,xmax
            if method == 'gridline':
                # find nearest gridline to specified lon:
                if len(lons.shape) == 1:
                    i0 = ma.argmin((lons[:]-lon)*(lons[:]-lon))
                    j0 = ma.argmin((lats[:]-xmin)*(lats[:]-xmin)) 
                    j1 = ma.argmin((lats[:]-xmax)*(lats[:]-xmax))                
                    print 'i0,[j0:j1] : ',i0,',[',j0,',',j1,']'
                    fld_xsec = fld.data[:,j0:j1+1,i0]
                    xpoints = lats[j0:j1+1]
                    xpoints_bounds = np.concatenate((lats_bounds[j0:j1+1,0],np.array([lons_bounds[j1+1,1]])))
                else:
                    lat_mid = 0.5 * (xmin + xmax)
                    (j_mid,i0) = np.unravel_index(ma.argmin((lats[:]-lat_mid)*(lats[:]-lat_mid)+(lons[:]-lon)*(lons[:]-lon)),lons.shape)
                    j0 = ma.argmin((lats[:,i0]-xmin)*(lats[:,i0]-xmin)) 
                    j1 = ma.argmin((lats[:,i0]-xmax)*(lats[:,i0]-xmax))                
                    print 'i0,[j0:j1] : ',i0,',[',j0,',',j1,']'
                    fld_xsec = fld.data[:,j0:j1+1,i0]
                    xpoints = lats[j0:j1+1,i0]
                    xpoints_bounds = np.concatenate((lats_bounds[j0:j1+1,i0,0],np.array([lats_bounds[j1+1,i0,1]])))
            elif method == 'dave':
                xsec_indices = np.where( ((lon-toldeg < lons) & (lons < lon+toldeg)) & ((lats > xmin) & (lats < xmax)) ) 
                print 'xsec_indices : ',zip(xsec_indices[0],xsec_indices[1])
                fld_xsec = fld.data[:,xsec_indices[0],xsec_indices[1]]
                xpoints = lats[xsec_indices]
                sort_indices = np.argsort(xpoints)
                xpoints = xpoints[sort_indices]
                xpoints_bounds = get_xpoints_bounds(xpoints)
                fld_xsec = fld_xsec[:,sort_indices]
#            elif method == 'pat':
#                start_lat = lat
#                start_lon = xmin
#                end_lat = lat
#                end_lon = xmax
#                lona = np.array((start_lon,end_lon))
#                lata = np.array((start_lat,end_lat))
#                [sec_xind,sec_yind,sec_lat,sec_lon,sec_dist] = find_section_indices(lona,lata,lats,lons)
#                fld_xsec = np.squeeze(fld.data[:,sec_yind,sec_xind])
#                xpoints = np.squeeze(lons[sec_yind,sec_xind])

#        print 'xpoints : ',xpoints[:]
#        print 'xpoints_bounds : ',xpoints_bounds[:]
#        print 'fld_xsec[47.:] : ',fld_xsec[47,:]
#        print 'fld_xsec.mask[47.:] : ',fld_xsec.mask[47,:]

        fld_xsec_list.append(fld_xsec.copy())
        if maskzero:
            fld_xsec_list[-1] = ma.masked_values(fld_xsec_list[-1],0.0)
        xpoints_list.append(xpoints.copy())

#############################
# 3. Do the plot
#############################

#    =========  SANITY CHECK FOR PLOTTING ========
#    print 'fld_xsec.shape : ',fld_xsec.shape
#    nz = fld_xsec.shape[0]
#    nz3 = nz/3
#    nz23 = 2*nz/3
#    print 'nz, nz3, nz23 : ',nz, nz3, nz23
#    print 'fld_xsec[:nz3,:].shape : ',fld_xsec[:nz3,:].shape
#    fld_xsec[:nz3,:] = 3.0e-5
#    fld_xsec[nz3:nz23,:] = 0.0
#    fld_xsec[nz23:,:] = -3.0e-5
#    print 'ma.min(fld_xsec), ma.max(fld_xsec) : ',ma.min(fld_xsec), ma.max(fld_xsec)
#    =============================================

    np.set_printoptions(threshold='nan')

    if type(mnfld) is not list:
        mnfld = nvar*[mnfld]
    elif len(mnfld) != nvar:
        raise Exception("Error: number of values for MNFLD (-f) doesn't match number of fields.")
    if type(mxfld) is not list:
        mxfld = nvar*[mxfld]
    elif len(mxfld) != nvar:
        raise Exception("Error: number of values for MXFLD (-F) doesn't match number of fields.")

    # loop over fields to plot (2 max)
    for var_name,plot_type,xpoints,fld_xsec,fac,mnfld_i,mxfld_i,nlevs_i in zip(var_names,plot_types[:nvar],xpoints_list,fld_xsec_list,factor,mnfld,mxfld,nlevs):

        if fac is not None:
            print 'multiplying field ',var_name,' by ',fac
            fld_xsec[:] = fac*fld_xsec[:]

        if mnfld_i is None:
            mnfld_i = fld_xsec.min()
        if mxfld_i is None:
            mxfld_i = fld_xsec.max()

        print 'mnfld_i, mxfld_i : ',mnfld_i,mxfld_i

        print 'nlevs_i : ',nlevs_i
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
        print 'levs_i : ',levs_i

        if colors is None and nvar == 2:
            # for a line contour plot over a filled contour plot default to black lines
            colors = "black"

        if colors is not None and nvar == 1:
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
#            if reverse_colors:
#                cmap = mcolors.ListedColormap(cmap(np.arange(255,-1,-1)))

        if colors is None and cmap is None:
            cmap = monty.clr_cmap('/home/h04/frsy/IDL/rainbow_diff_nice.clr')

        for depthname in ['depth', 'deptht', 'depthu', 'depthv', 'depthw']:
            try:
                depth = fld.coord( var_name = depthname ) 
                depth_bounds = np.concatenate((depth.bounds[:,0],np.array([depth.bounds[-1,1]])))
            except iris.exceptions.CoordinateNotFoundError:
                pass
            else:
                break
        else:
            raise iris.exceptions.CoordinateNotFoundError

        # Keep the original mask before we fill fld_xsec prior to contouring.
        mask_to_plot = fld_xsec.mask.copy()
        # Plot the mask.
        # Make a colormap with grey on the end instead of black
        cmap_mask = matplotlib.colors.ListedColormap(cm.binary(np.arange(128)))
        clev_mask = plt.pcolormesh(xpoints_bounds, depth_bounds, mask_to_plot, cmap=cmap_mask)

        if plot_type == 'l':
            fld_xsec_filled = fill_field(fld_xsec)
            llev = plt.contour(xpoints,depth.points,fld_xsec_filled, colors=colors, cmap=cmap, levels=levs_i)
#            manual_locations=[(35,1000)]
            if fmt is None:
                fmt='%1.1f'
            plt.clabel(llev,inline=True,manual=False,fontsize=7,fmt=fmt)
        elif plot_type == 'c':
            if block:
                if colors is not None:
                    (cmap_block,norm_block) = matplotlib.colors.from_levels_and_colors(levs,colors)
                else:
                    cmap_block = cmap
                clev = plt.pcolormesh(xpoints_bounds, depth_bounds, fld_xsec, cmap=cmap_block, vmin=mnfld_i, vmax=mxfld_i)
            else:
                fld_xsec_filled = fill_field(fld_xsec)
                clev=plt.contourf(xpoints,depth.points,fld_xsec_filled,cmap=cmap,levels=levs_i)
        else:
            raise Exception("Error: unknown plot_type: "+plot_type)

        if not nobar and plot_type == 'c':
            cax = plt.colorbar(clev,orientation='horizontal')
            if block:
#                print 'norm_block(levs) : ',norm_block(levs)
#                cax.set_ticks(norm_block(levs))
                print 'levs_i : ',levs_i
                cax.set_ticks(levs_i)
            else:
                cax.set_ticks(levs_i)
            if scientific:
                labels=["%1.2e" %lev for lev in levs_i]
            else:
                labels=["%1.3f" %lev for lev in levs_i]
            cax.ax.set_xticklabels(labels,rotation=45)
    
    if title is not None:
        plt.gca().set_title(textwrap.fill(title,60))
    if reverseX:
        gca().set_xlim(gca().get_xlim()[::-1])
    if depthmax is not None:
        plt.gca().set_ylim([depthmax,0.0])
    else:
        plt.gca().invert_yaxis()

    plt.gca().set_ylabel('depth (m)')
    if lon is not None:
        if kilometres:
            plt.gca().set_xlabel('south-north distance (km)')
        else:
            plt.gca().set_xlabel('latitude (degrees north)')
    elif lat is not None:
        if kilometres:
            plt.gca().set_xlabel('west-east distance (km)')
        else:
            plt.gca().set_xlabel('longitude (degrees east)')

    if outfile is not None:
        matplotlib.rcParams['font.size'] =8
        plt.savefig(outfile,dpi=200)
    elif not noshow:
        plt.show()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file",action="store",dest="filenames", nargs="+",
                    help="name(s) of input file(s) for scalar field 1 and scalar field 2")
    parser.add_argument("-v", "--var",action="store",dest="var_names", nargs="+",
                    help="name(s) of field(s) to plot.") 
    parser.add_argument("-Z", "--factor", action="store",dest="factor",type=float,nargs="+",
                    help="multiplicative factor(s) to apply to plotted field(s)")
    parser.add_argument("-z", "--maskzero", action="store_true",dest="maskzero",
                    help="mask zero values in field")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",
                    help="plot to output file - filetype by extension. If unset plot to GUI.")
    parser.add_argument("-L", "--lineplot", action="store_true",dest="lineplot",
                    help="plot line contours instead of colour-filled contours for *first* field.")
    parser.add_argument("-p", "--phi", action="store",dest="lat",type=float,
                    help="latitude of cross section")
    parser.add_argument("-l", "--lambda", action="store",dest="lon",type=float,
                    help="longitude of cross section")
    parser.add_argument("-I", "--index_i", action="store",dest="index_i",type=int,
                    help="i-index of cross section")
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
    parser.add_argument("-f", "--mnfld", action="store",dest="mnfld",type=float,nargs="+",
                    help="minimum field value to plot for colour filled contouring")
    parser.add_argument("-F", "--mxfld", action="store",dest="mxfld",type=float,nargs="+",
                    help="maximum field value to plot for colour filled contouring")
    parser.add_argument("-D", "--depthmax", action="store",dest="depthmax",type=float,
                    help="maximum depth for plot")
    parser.add_argument("-t", "--title", action="store",dest="title",
                    help="title for plot")
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
    parser.add_argument("-k", "--kilometres", action="store_true",dest="kilometres",
                    help="x-axis to be labelled as kilometres rather than degrees (for f-plane models etc)")
    parser.add_argument("-Q", "--toldeg", action="store",dest="toldeg",type=float,
                    help="tolerance in degrees when calculating cross section (default=0.1)")
    parser.add_argument("-e", "--scientific", action="store_true",dest="scientific",
                    help="scientific format for colour bar labels")
    parser.add_argument("-B", "--block", action="store_true",dest="block",
                    help="block plotting rather than contouring")

    args = parser.parse_args()

    plot_nemo_section(filenames=args.filenames,var_names=args.var_names, outfile=args.outfile, title=args.title, 
    lat=args.lat, lon=args.lon, index_i=args.index_i, index_j=args.index_j,levs=args.levs, nlevs=args.nlevs, rec=args.rec, 
    xmin=args.xmin, xmax=args.xmax,mnfld=args.mnfld, mxfld=args.mxfld, depthmax=args.depthmax, method=args.method, 
    logscale=args.logscale, reverseX=args.reverseX, kilometres=args.kilometres, toldeg=args.toldeg, 
    noshow=args.noshow, nobar=args.nobar, scientific=args.scientific, block=args.block, lineplot=args.lineplot, 
    colors=args.colors, cmap=args.cmap, factor=args.factor,fmt=args.fmt,maskzero=args.maskzero)
       

