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

@author: Dave Storkey and Pat Hyder
'''

import socket
import matplotlib
if 'spice' in socket.gethostname():
    # Note this disables plt.show()
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as col
import iris
import iris.plot as iplt
import monty
import numpy as np
import general_tools as gt
from pylab import *
import textwrap
# import glob

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
                 lat=None,lon=None,nlevs=15,rec=0,mnfld=None,mxfld=None,
                 xmin=None,xmax=None,depthmax=None,method='dave',logscale=False, 
                 reverseX=False,kilometres=None,outfile=None,toldeg=None,
                 noshow=False,nobar=False,scientific=False,block=False,
                 lineplot=False):

    print '>>> rec is : ',rec 
    print '>>> method is : ',method

    if lineplot:
        plot_types=['l','c']
    else:
        plot_types=['c','l']

#############################
# 1. Read in the data
#############################

# HARDWIRED FOR ORCA025 FOR NOW!!
    gridfile='/data/cr1/hadtd/CMA/HadGEM3/ORCA025/mesh_mask_GO5.nc'
    maskfile='/project/ujcc/CDFTOOLS/mesh_ORCA025L75/subbasins_orca025_070909_rename.nc'
    maskvar='glbmsk'

    mask=(gt.read_nc_var(maskfile,maskvar)*gt.read_nc_var(gridfile,'tmaskutil'))

    if var_names is None:
        raise Exception("Error: must specify at least one field to plot.")
    elif len(var_names) > 2:
        raise Exception("Error: can specify at most two fields to plot.")
    nvar = len(var_names)

    if filenames is None:
        raise Exception("Error: must specify at least one input file.")
    elif len(filenames) == 1 and len(var_names) == 2:
        # if two fields and only one input file try to read both fields from the same file:
        filenames=[filenames,filenames] 

    fld_list = []
    for filename, varname in zip(filenames,var_names):
        fld_read = gt.read_cube(filename,varname)

# Load required record from file.
# NB. Assuming the field has 3 spatial dimensions here.
        if fld_read.ndim == 4:
            fld_list.append(fld_read[rec].copy())
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
        except iris.exceptions.CoordinateNotFoundError:
            print 'reading latitudes as field'
            lats=iris.load_cube(filename,'latitude').data
        else:
            print 'reading latitudes as auxillary coordinate'

        try:
            lons=fld.coord('longitude').points
        except iris.exceptions.CoordinateNotFoundError:
            print 'reading longitudes as field'
            lons=iris.load_cube(filename,'longitude').data
        else:
            print 'reading longitudes as auxillary coordinate'

        select_mask = np.where(lons < 0.0,1,0)
        lons = lons + 360.0*select_mask

        if lat is not None:
            if method == 'gridline':
                # interpret "lat" as y-index and xmin/xmax as min/max x-index
                if xmin is None:
                    xmin = 0
                if xmax is None:
                    xmax = fld.data.shape[2]
                fld_xsec = fld.data[:,lat,xmin:xmax+1]
                xpoints = lons[lat,xmin:xmax+1]
            elif method == 'dave':
                if xmin is None:
                    xmin = lons.min()
                if xmax is None:
                    xmax = lons.max()
                print 'lons.min(),lons.max()',lons.min(),lons.max()
                print 'xmin,xmax',xmin,xmax
                flat_indices = np.flatnonzero( np.where(lat-toldeg < lats,1,0) * np.where(lats < lat+toldeg,1,0) *
                                               np.where(lons > xmin,1,0) * np.where(lons < xmax,1,0) )
                print "len(flat_indices)",len(flat_indices)
                fld_flatij = np.reshape(fld.data,[fld.shape[0],fld.shape[1]*fld.shape[2]])
                fld_xsec = fld_flatij[:,flat_indices]
                print 'lons.shape',lons.shape
                xpoints_all = np.reshape(lons,lons.shape[0]*lons.shape[1])
                xpoints = xpoints_all[flat_indices]
                sort_indices = np.argsort(xpoints)
                xpoints = xpoints[sort_indices]
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
    
        if lon is not None:
            if method == 'gridline':
                # interpret "lon" as x-index and xmin/xmax as min/max y-index
                if xmin is None:
                    xmin = 0
                if xmax is None:
                    xmax = fld.data.shape[1]
                fld_xsec = fld.data[:,xmin:xmax+1,lon]
                xpoints = lats[xmin:xmax+1,lon]
            elif method == 'dave':
                if xmin is None:
                    xmin = lats.min()
                if xmax is None:
                    xmax = lats.max()
                print 'lons.min(),lons.max()',lons.min(),lons.max()
                print 'lats.min(),lats.max()',lats.min(),lats.max()
                print 'xmin,xmax',xmin,xmax
                print 'lon : ',lon
#                flat_indices = np.flatnonzero( np.where(lon-toldeg < lons,1,0) )
#                flat_indices = np.flatnonzero( np.where(lons < lon+toldeg,1,0) )
                flat_indices = np.flatnonzero( np.where(lon-toldeg < lons,1,0) * np.where(lons < lon+toldeg,1,0) *
                                               np.where(lats > xmin,1,0) * np.where(lats < xmax,1,0) )
                print "len(flat_indices)",len(flat_indices)
                fld_flatij = np.reshape(fld.data,[fld.shape[0],fld.shape[1]*fld.shape[2]])
                fld_xsec = fld_flatij[:,flat_indices]
                print 'lats.shape',lats.shape
                xpoints_all = np.reshape(lats,lats.shape[0]*lats.shape[1])
                xpoints = xpoints_all[flat_indices]
                sort_indices = np.argsort(xpoints)
                xpoints = xpoints[sort_indices]
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

            fld_xsec_list.append(fld_xsec.copy())
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

    for plot_type,xpoints,fld_xsec,mnfld_i,mxfld_i in zip(plot_types[:nvar],xpoints_list,fld_xsec_list,mnfld,mxfld):

        if mnfld_i is None:
            mnfld_i = fld_xsec.min()
        if mxfld_i is None:
            mxfld_i = fld_xsec.max()

        print 'nlevs : ',nlevs
        if logscale:
            tiny_factor=0.0001
            if mnfld_i == 0.0:
                mnfld_i = tiny_factor * mxfld_i
            if mxfld_i == 0.0:
                mxfld_i = tiny_factor * mnfld_i
            if mnfld_i < 0.0 and mxfld_i < 0.0:
                levs = -np.logspace(np.log10(-mnfld_i),np.log10(-mxfld_i),nlevs+1) 
            if mnfld_i > 0.0 and mxfld_i > 0.0:
                levs = np.logspace(np.log10(mnfld_i),np.log10(mxfld_i),nlevs+1) 
            if mnfld_i < 0.0 and mxfld_i > 0.0:
                nlevs_neg = int(abs(mnfld_i)*nlevs/(mxfld_i-mnfld_i))
                levs_neg = -np.logspace(np.log10(-mnfld_i),np.log10(-mnfld_i*tiny_factor),nlevs+1) 
                nlevs_pos = int(mxfld_i*nlevs/(mxfld_i-mnfld_i))
                levs_pos = np.logspace(np.log10(mxfld_i*tiny_factor),np.log10(mxfld_i),nlevs+1) 
                levs = levs_neg + levs_pos
        else:
            levs = np.linspace(mnfld_i,mxfld_i,nlevs+1) 

        if plot_type == 'c':
            cmap = monty.clr_cmap('/home/h05/hadtd/IDL/rainbow_diff_nice.clr')
            ind = np.intp(np.rint(np.linspace(2,255,nlevs+1))).tolist()    
# ensure we have white space in the middle for an odd number of levels.
            if 2*(nlevs/2) != nlevs:
                ind[nlevs/2] = 127
            print 'ind',ind
            colorlist=cmap(ind)
            print 'levs : ',levs
#        my_cmap = matplotlib.colors.ListedColormap(colorlist, 'my_cmap')
    
        if plot_type == 'c':
            plt.gca().patch.set_facecolor('grey')
        for depthname in ['depth', 'deptht', 'depthu', 'depthv', 'depthw']:
            try:
                depth = fld.coord( var_name = depthname ) 
            except iris.exceptions.CoordinateNotFoundError:
                pass
            else:
                break
        else:
            raise iris.exceptions.CoordinateNotFoundError

        if plot_type == 'l':
            print 'len(xpoints),len(depth.points),fld_xsec.shape : ',len(xpoints.data),len(depth.points),fld_xsec.shape
            llev = plt.contour(xpoints,depth.points,fld_xsec, colors='black',levels=levs)
    #        plt.clabel(llev,inline=False,manual=True)
        elif plot_type == 'c':
            if block:
                (cmap_block,norm_block) = matplotlib.colors.from_levels_and_colors(levs,colorlist[0:-1])
                clev = plt.pcolormesh(xpoints, depth.points, fld_xsec,  cmap=cmap_block, norm=norm_block)
            else:
                clev=plt.contourf(xpoints,depth.points,fld_xsec,colors=colorlist[:-1],levels=levs)
        else:
            raise Exception("Error: unknown plot_type: "+plot_type)

        if not nobar and plot_type == 'c':
            cax = plt.colorbar(clev,orientation='horizontal')
            if block:
#                print 'norm_block(levs) : ',norm_block(levs)
#                cax.set_ticks(norm_block(levs))
                cax.set_ticks(levs)
            else:
                cax.set_ticks(levs)
            if scientific:
                labels=["%1.2e" %lev for lev in levs]
            else:
                labels=["%1.3f" %lev for lev in levs]
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
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",
                    help="plot to output file - filetype by extension. If unset plot to GUI.")
    parser.add_argument("-C", "--lineplot", action="store_true",dest="lineplot",
                    help="plot line contours instead of colour-filled contours for *first* field.")
    parser.add_argument("-l", "--latitude", action="store",dest="lat",type=float,
                    help="latitude of cross section")
    parser.add_argument("-L", "--longitude", action="store",dest="lon",type=float,
                    help="longitude of cross section")
    parser.add_argument("-n", "--nlevs", action="store",dest="nlevs",type=int,default=15,
                    help="number of contour levels")
    parser.add_argument("-r", "--rec", action="store",dest="rec",type=int,default=0,
                    help="record number to plot (defaults to zero)")
    parser.add_argument("-x", "--xmin", action="store",dest="xmin",type=float,
                    help="start of x-axis (min lat or min lon as applicable)")
    parser.add_argument("-X", "--xmax", action="store",dest="xmax",type=float,
                    help="end of x-axis (max lat or max lon as applicable)")
    parser.add_argument("-f", "--mnfld", action="store",dest="mnfld",type=float,default=None,nargs="+",
                    help="minimum field value to plot for colour filled contouring")
    parser.add_argument("-F", "--mxfld", action="store",dest="mxfld",type=float,default=None,nargs="+",
                    help="maximum field value to plot for colour filled contouring")
    parser.add_argument("-D", "--depthmax", action="store",dest="depthmax",type=float,
                    help="maximum depth for plot")
    parser.add_argument("-t", "--title", action="store",dest="title",
                    help="title for plot")
    parser.add_argument("-m", "--method", action="store",dest="method",default="dave",
                    help="method for extracting section ('dave' or 'pat' or 'gridline')")
    parser.add_argument("-G", "--logscale", action="store_true",dest="logscale",default=False,
                    help="use colour logscale")
    parser.add_argument("-R", "--reverseX", action="store_true",dest="reverseX",default=False,
                    help="reverse the x axis")
    parser.add_argument("-N", "--noshow", action="store_true",dest="noshow",default=False,
                    help="Do not show GUI window even if output filename unset.")
    parser.add_argument("-b", "--nobar", action="store_true",dest="nobar",default=False,
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
    lat=args.lat, lon=args.lon, nlevs=args.nlevs, rec=args.rec, xmin=args.xmin, xmax=args.xmax, 
    mnfld=args.mnfld, mxfld=args.mxfld, depthmax=args.depthmax, method=args.method, logscale=args.logscale, 
    reverseX=args.reverseX, kilometres=args.kilometres, toldeg=args.toldeg, noshow=args.noshow, nobar=args.nobar,
    scientific=args.scientific, block=args.block, lineplot=args.lineplot)
       

