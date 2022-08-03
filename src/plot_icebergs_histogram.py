#! /usr/bin/env python
'''
Routine to plot histogram of iceberg variables 
from an iceberg restart file. 

Created Nov 2020

@author: Dave Storkey
'''

import netCDF4 as nc
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import textwrap

def plot_icebergs_histogram(infiles=None,outfile=None,var=None,rec=None,nbins=None,logy=None,
                            title=None,legend=None):

    if rec is None:
        rec=0

    if nbins is None:
        nbins=50

    for filename in infiles:

        with nc.Dataset(filename,'r') as file_in:
            uvel = file_in.variables['uvel'][:]
            vvel = file_in.variables['vvel'][:]
            ntraj = uvel.shape
            print ('number of trajectories: ',ntraj)

        speed = np.sqrt(uvel*uvel + vvel*vvel)

        hist, bins = np.histogram(speed, bins=nbins)
        width = 0.8 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        plt.bar(center, hist, align='center', width=width)
    
    if logy:
        plt.yscale('log')

    plt.gca().set_xlabel("iceberg speed (m/s)")

    if legend:
        plt.gca().legend(legend)

    if title is not None:
        title=textwrap.fill(title,70)
        plt.gcf().suptitle(title, fontsize=12, y=0.95)    

    if outfile is not None:
        matplotlib.rcParams['font.size'] =8
        plt.savefig(outfile,dpi=200)
    else:
        plt.show()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infiles", action="store", dest="infiles", nargs="+",
                        help="name of input file(s)")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",default=None,
                    help="plot to output file - filetype by extension. If unset plot to GUI.")
    parser.add_argument("-v", "--var", action="store",dest="var",default=None,
                    help="variable to plot")
    parser.add_argument("-r", "--rec", action="store",dest="rec",type=int,default=0,
                    help="record number to plot (defaults to zero)")
    parser.add_argument("-n", "--nbins", action="store",dest="nbins",type=int,default=None,
                    help="number of bins for histogram (default 50)")
    parser.add_argument("-G", "--logy", action="store_true",dest="logy",default=False,
                    help="log scale for y-axis")
    parser.add_argument("-L", "--legend", action="store",dest="legend",nargs="+",default=False,
                    help="legend labels")
    parser.add_argument("-t", "--title", action="store",dest="title",default=None,
                    help="title for plot")

    args = parser.parse_args()

    plot_icebergs_histogram(infiles=args.infiles,var=args.var,outfile=args.outfile,nbins=args.nbins,logy=args.logy,
                            rec=args.rec,title=args.title,legend=args.legend)        
