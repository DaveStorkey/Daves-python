#! /usr/bin/env python
'''
Routine to visualise cross sections of the partial steps bathymetry
superimposed on the model grid (x-z grid or y-z grid). 

Requires a bathymetry file and a meshmask file for the model grid. 

@author: Dave Storkey
@date: June 2016
'''

import netCDF4 as nc
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import textwrap

def partial_steps(bathyfile=None, meshfile=None, outfile=None, zco=False, zps=False, u_slopes=False,
                  xmin=None, xmax=None, ymin=None, ymax=None, e3zps_min=None, e3zps_rat=None,
                  depthmin=None, depthmax=None, title=None ):

    bathydata = nc.Dataset(bathyfile,mode="r")
    bathy = np.copy(bathydata.variables["Bathymetry"][:])
    bathydata.close()

    if meshfile is not None:
        meshdata = nc.Dataset(meshfile,mode="r")
        gdepw_0 = np.copy(meshdata.variables["gdepw_1d"][:])
        e3t_0 = np.copy(meshdata.variables["e3t_1d"][:])
        # create an extended array to make the projection onto the model vertical grid work below...
        gdepw_x = np.append(gdepw_0,gdepw_0[-1]+1000.0)

        outdata = nc.Dataset("test_out.nc","w")
        for dim in meshdata.dimensions.keys():
            outdata.createDimension(dim, len(meshdata.dimensions[dim])) 

        meshdata.close()

    bathy_to_plot = bathy

    if zco:
        # Create a lego z-level bathymetry from the bathymetry dataset read in. 
        bathy_zco = np.zeros(bathy.shape)
        mbathy = np.searchsorted(gdepw_0, bathy[:], side="left")
        bathy_zco = np.where( (mbathy > 0) 
         & ( (mbathy == len(gdepw_0)) | (np.abs(bathy-gdepw_x[mbathy-1]) < np.abs(bathy-gdepw_x[mbathy]) ) ),
         gdepw_x[mbathy-1], gdepw_x[mbathy] )
        bathy_to_plot = bathy_zco

    if zps:
        if e3zps_min is None:
            print "WARNING: e3zps_min unset. Using GO6 standard value of 25.0m"
            e3zps_min=25.0
        if e3zps_rat is None:
            print "WARNING: e3zps_rat unset. Using GO6 standard value of 0.2"
            e3zps_rat=0.2
        # Create a partial-steps z-level bathymetry from the bathymetry dataset read in.
        # This follows zgr_zps in domzgr.F90.
        bathy_zps = np.zeros(bathy.shape)
        zdepth = gdepw_0[:] + np.minimum( e3zps_min, e3t_0[:]*e3zps_rat )
        outdata.createVariable('zdepth',datatype='f',dimensions=("z"))
        outdata.variables["zdepth"][:] = zdepth[:]
        mbathy = np.searchsorted(zdepth, bathy[:], side="left")
        outdata.createVariable("mbathy",datatype="f",dimensions=("y","x"))
        outdata.variables["mbathy"][:] = mbathy[:]
        outdata.close()
        bathy_zps = np.where( (mbathy > 0) 
         & ( (mbathy == len(gdepw_0)) | ( bathy-gdepw_x[mbathy] <= 0.0 ) ),
         bathy, gdepw_x[mbathy] )
        bathy_to_plot = bathy_zps

    x_section=None
    y_section=None
    if xmin is not None and xmax is not None:
        if ymin is not None and ymax is not None:
            raise Exception("Don't know what to do. You specified all of xmin,xmax,ymin,ymax.")
        elif ymin is not None:
            y_section = ymin
        elif ymax is not None:
            y_section = ymax
        else:
            raise Exception("You need to specify the y-value for the cross section using ymin or ymax")
    if ymin is not None and ymax is not None:
        if xmin is not None:
            x_section = xmin
        elif xmax is not None:
            x_section = xmax
        else:
            raise Exception("You need to specify the x-value for the cross section using xmin or xmax")

    if x_section:
        bathyslice_t = bathy_to_plot[ymin:ymax,x_section]
        if u_slopes:
            mbathyslice = mbathy[ymin:ymax,x_section]
            xvals = np.zeros(3*len(bathyslice_t)-2)
            xvals[0::3] = np.arange(ymin,ymax) 
            xvals[1::3] = np.arange(ymin,ymax-1)+0.4999 
            xvals[2::3] = np.arange(ymin,ymax-1)+0.5001 
        else:
            xvals = np.arange(ymin,ymax)
        xlabel="y index"
        if title is None:
            title="section along x = "+str(x_section)
        else:
            title=title+": section along x = "+str(x_section)
    elif y_section:
        bathyslice_t = bathy_to_plot[y_section,xmin:xmax]
        if u_slopes:
            mbathyslice = mbathy[y_section,xmin:xmax]
            xvals = np.zeros(3*len(bathyslice_t)-2)
            xvals[0::3] = np.arange(xmin,xmax) 
            xvals[1::3] = np.arange(xmin,xmax-1)+0.4999 
            xvals[2::3] = np.arange(xmin,xmax-1)+0.5001 
        else:
            xvals = np.arange(xmin,xmax)
        xlabel="x index"
        if title is None:
            title="section along y = "+str(y_section)
        else:
            title=title+": section along y = "+str(y_section)
   
    if u_slopes:
        bathyslice_u1 = np.zeros(len(bathyslice_t)-1)
        bathyslice_u2 = np.zeros(len(bathyslice_t)-1)
        jump1 = bathyslice_t[1:]-bathyslice_t[0:-1]
        jump1 = np.minimum( jump1, gdepw_x[mbathyslice[0:-1]]-bathyslice_t[0:-1] )
        jump1 = np.maximum( jump1, gdepw_x[mbathyslice[0:-1]-1]-bathyslice_t[0:-1] )
        bathyslice_u1[:] = bathyslice_t[0:-1] + 0.5*jump1
        jump2 = bathyslice_t[0:-1]-bathyslice_t[1:]
        jump2 = np.minimum( jump2, gdepw_x[mbathyslice[1:]]-bathyslice_t[1:] )
        jump2 = np.maximum( jump2, gdepw_x[mbathyslice[1:]-1]-bathyslice_t[1:] )
        bathyslice_u2[:] = bathyslice_t[1:] + 0.5*jump2
        bathyslice = np.zeros(len(bathyslice_t)+len(bathyslice_u1)+len(bathyslice_u2))
        bathyslice[0::3] = bathyslice_t
        bathyslice[1::3] = bathyslice_u1
        bathyslice[2::3] = bathyslice_u2
    else:
        bathyslice = bathyslice_t

    if depthmin is None:
        depthmin=0.0
    if depthmax is None:
        depthmax=bathyslice.max()

    print 'xvals : ',xvals

    if u_slopes:
        plt.plot(xvals,bathyslice)
    else:
        plt.step(xvals,bathyslice,where="mid")
    plt.gca().set_xlabel(xlabel)
    plt.gca().set_ylim([depthmax,depthmin])
    plt.gca().set_ylabel('depth (m)')
    if meshfile is not None:
        yticks = gdepw_0[np.where( (depthmin < gdepw_0) & (gdepw_0 < depthmax) )]
        plt.gca().set_yticks(yticks)
    plt.grid()
    plt.gca().set_title(textwrap.fill(title,60))

    if outfile is not None:
        matplotlib.rcParams['font.size'] =8
        plt.savefig(outfile,dpi=200)
    else:
        plt.show()


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bathy", action="store",dest="bathyfile",default=None,
                    help="Name of bathymetry file. (Assume field is called Bathymetry)")
    parser.add_argument("-m", "--mesh", action="store",dest="meshfile",default=None,
                    help="Name of meshmask file.")
    parser.add_argument("-Z", "--zco", action="store_true",dest="zco", default=False,
                    help="calculate and plot lego z-level bathymetry equivalent to input bathymetry")
    parser.add_argument("-P", "--zps", action="store_true",dest="zps", default=False,
                    help="calculate and plot partial-steps bathymetry equivalent to input bathymetry")
    parser.add_argument("-U", "--u_slopes", action="store_true",dest="u_slopes", default=False,
                    help="plot cross section with sloping bathymetry at U- and V- points as per FVPS scheme")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",default=None,
                    help="plot to output file - filetype by extension. If unset plot to GUI.")
    parser.add_argument("-x", "--xmin", action="store",dest="xmin",type=float,default=None,
                    help="start of x-axis or x-value of cross section ")
    parser.add_argument("-X", "--xmax", action="store",dest="xmax",type=float,default=None,
                    help="end of x-axis  or x-value of cross section")
    parser.add_argument("-y", "--ymin", action="store",dest="ymin",type=float,default=None,
                    help="start of y-axis or y-value of cross section ")
    parser.add_argument("-Y", "--ymax", action="store",dest="ymax",type=float,default=None,
                    help="end of y-axis  or y-value of cross section")
    parser.add_argument("-d", "--depthmin", action="store",dest="depthmin",type=float,default=None,
                    help="minimum depth for plot")
    parser.add_argument("-D", "--depthmax", action="store",dest="depthmax",type=float,default=None,
                    help="maximum depth for plot")
    parser.add_argument("-p", "--e3zps_min", action="store",dest="e3zps_min",type=float,default=None,
                    help="minimum thickness for partial cell (metres)")
    parser.add_argument("-q", "--e3zps_rat", action="store",dest="e3zps_rat",type=float,default=None,
                    help="minimum ratio of partial cell thickness to full cell thickness")
    parser.add_argument("-t", "--title", action="store",dest="title",default=None,
                    help="title for plot")
 
    args = parser.parse_args()

    partial_steps(bathyfile=args.bathyfile, meshfile=args.meshfile, zco=args.zco, zps=args.zps, u_slopes=args.u_slopes,
                  outfile=args.outfile, depthmin=args.depthmin, depthmax=args.depthmax,
                  xmin=args.xmin, xmax=args.xmax, ymin=args.ymin, ymax=args.ymax, title=args.title,
                  e3zps_min=args.e3zps_min, e3zps_rat=args.e3zps_rat )

