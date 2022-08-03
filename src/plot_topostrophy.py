#! /usr/bin/env python
'''
Routine to plot topostrophy in the style of 
Penduff et al (2007) Figure 7.

Dave Storkey
Aug 2018
'''

import argparse
import netCDF4 as nc
import numpy.ma as ma
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator
import monty

def plot_topostrophy(infile,outfile=None,nlevs=None,mnfld=None,mxfld=None,ratio=False):

    with nc.Dataset(infile,'r') as topo_in:
        depths=topo_in.variables['depthx'][:]
        topostr=topo_in.variables['topostrophy'][:]

    # depths in km:
    depths[:] = depths[:]*0.001

    if nlevs is None:
        nlevs = 6
    if mnfld is None:
        mnfld = ma.min(topostr)
    if mxfld is None:
        mxfld = ma.max(topostr)

    locator = MaxNLocator(nlevs+1)
    locator.create_dummy_axis()
    locator.set_bounds(mnfld, mxfld)
    levels = locator()

    plt.contour(depths,depths,topostr,colors='black',linewidths=0.5,levels=levels)
    if ratio:
        cmap_in = monty.clr_cmap('/home/h05/hadtd/IDL/rainbow_diff_nice.clr')
    else:
        cmap_in = cm.binary
    plt.contourf(depths,depths,topostr,cmap=cmap_in,levels=levels)

    plt.gca().invert_yaxis()
    plt.gca().set_ylabel("local depth (km)")
    plt.gca().set_xlabel("bottom depth (km)")

    cax = plt.colorbar(orientation='horizontal')

    if outfile is not None:
        matplotlib.rcParams['font.size'] =8
        plt.savefig(outfile,dpi=200)
    else:
        plt.show()

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile",help="input file name")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",default=None,
                    help="plot to output file - filetype by extension. If unset plot to GUI.")
    parser.add_argument("-n", "--nlevs", action="store",dest="nlevs",type=int,default=15,
                    help="number of contour levels")
    parser.add_argument("-f", "--mnfld", action="store",dest="mnfld",type=float,default=None,
                    help="minimum field value to plot for colour filled contouring")
    parser.add_argument("-F", "--mxfld", action="store",dest="mxfld",type=float,default=None,
                    help="maximum field value to plot for colour filled contouring")
    parser.add_argument("-R", "--ratio", action="store_true",dest="ratio",
                    help="ratio of two fields: use rainbow colours")

    args = parser.parse_args()

    plot_topostrophy(args.infile,outfile=args.outfile,nlevs=args.nlevs,mnfld=args.mnfld,mxfld=args.mxfld,ratio=args.ratio)

