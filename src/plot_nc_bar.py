#! /usr/bin/env python
"""
    Script to plot out a bar chart from a 1D vertical netcdf variable.
    Options for plottype keyword:
        'levels'     : (default) Plot a bar chart of values of field level by level.
        'depths'     : Plot a bar chart with depths as the vertical coordinate
                       and the area of the bars proportional to the values in 
                       field. Note that this requires deptht and e3t variables
                       in the input file. 
        'cumulative' : Plot the cumulative value of field summed from the surface 
                       with level index on the y-axis.
"""

__author__ = 'Dave Storkey (dave.storkey@metoffice.gov.uk)'
__date__ = '$Date: 2013/3/5 $'

import sys
import numpy as np
import matplotlib.pyplot as plt
import netCDF4

def plot_nc_bar(datafile,field,plottype='levels',maxdepth=None,
                xmin=None,xmax=None,ymin=None,ymax=None):

    dataid = netCDF4.Dataset(datafile, mode='r')
    fld = np.squeeze(np.copy(dataid.variables[field]))
    imax = len(fld) - 1
    if plottype == 'depths' or maxdepth is not None:
        thickness = np.copy(dataid.variables['e3t'])
        depths = np.zeros(len(thickness))
        for i in range(len(thickness)):
           depths[i] = thickness[0:i+1].sum()

        if maxdepth is not None:
            for i in range(len(thickness)):
                if depths[i] > maxdepth:
                    imax = i
                    break

    if plottype == 'levels':
        plt.xlabel('heat content change (W/m2 equivalent)')
        plt.xlim(-15.0,15.0)
    elif plottype == 'depths':
        fld = fld/thickness
        plt.xlabel('heat content change per unit depth (W/m3 equivalent)')
    elif plottype == 'cumulative':
        fld = np.cumsum(fld)
        plt.xlabel('heat content change (W/m2 equivalent) \n cumulative value summed from surface')
        plt.xlim(-160.0,160.0)

    if plottype == 'levels' or plottype == 'cumulative':
        pos = [-x for x in range(len(fld))]
        height = np.ones(len(fld))
        plt.ylabel('model level index')
    else:
        pos = [-x for x in depths]
        height = thickness
        plt.ylabel('depth (metres)')

    if xmin is not None:
        plt.gca().set_xlim(left=xmin)
    if xmax is not None:
        plt.gca().set_xlim(right=xmax)
    if ymin is not None:
        plt.gca().set_ylim(left=ymin)
    if ymax is not None:
        plt.gca().set_ylim(right=ymax)

    left=0
    plt.bar(left, height[0:imax], width=fld[0:imax], bottom=pos[0:imax], color='r',orientation='horizontal')
    plt.show()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="name of input file")
    parser.add_argument("field", help="name of field to plot")
    parser.add_argument("-l", "--levels", action="store_const",dest="plottype",const="levels",
                    default="depths",help="depth level to plot")
    parser.add_argument("-d", "--depths", action="store_const",dest="plottype",const="depths",
                    default="depths",help="depth level to plot")
    parser.add_argument("-c", "--cumulative", action="store_const",dest="plottype",const="cumulative",
                    default="depths",help="depth level to plot")
    parser.add_argument("-D", "--maxdepth", action="store",dest="maxdepth",type=int,default=None,
                    help="time level to plot")
    parser.add_argument("-x", "--xmin", action="store",dest="xmin",type=float,default=None,
                    help="start of x-axis")
    parser.add_argument("-X", "--xmax", action="store",dest="xmax",type=float,default=None,
                    help="end of x-axis")
    parser.add_argument("-y", "--ymin", action="store",dest="ymin",type=float,default=None,
                    help="start of y-axis")
    parser.add_argument("-Y", "--ymax", action="store",dest="ymax",type=float,default=None,
                    help="end of y-axis")
    args = parser.parse_args()

    plot_nc_bar(args.filename,args.field,plottype=args.plottype,maxdepth=args.maxdepth,
                xmin=args.xmin,xmax=args.xmax,ymin=args.ymin,ymax=args.ymax)        

#if __name__ == "__main__":
#    if len(sys.argv) == 3:
#        plot_nc_bar(datafile=sys.argv[1],field=sys.argv[2])
#    elif len(sys.argv) == 4:
#        plot_nc_bar(datafile=sys.argv[1],field=sys.argv[2],plottype=sys.argv[3])
