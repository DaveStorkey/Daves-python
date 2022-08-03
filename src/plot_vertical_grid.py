#! /usr/bin/env python
'''
Plot vertical grid(s) from domcfg file(s)
as depth versus depth level number.

Created June 2019

@author: Dave Storkey
'''

import netCDF4 as nc
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def plot_vertical_grid(infiles,outfile=None):

    for infile in infiles:
        indata = nc.Dataset(infile,'r')
        e3t_1d = indata.variables['e3t_1d'][:]
        depths = np.cumsum(e3t_1d)
        plt.plot(depths)      

    plt.gca().set_xlabel('level index')
    plt.gca().set_ylabel('depth (m)')
    
    if outfile is not None:
        matplotlib.rcParams['font.size'] =8
        plt.savefig(outfile,dpi=200)
    else:
        plt.show()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", action="store",dest="infiles",nargs="+",
                    help="input domain cfg file(s)")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",default=None,
                    help="plot to output file - filetype by extension. If unset plot to GUI.")

    args = parser.parse_args()

    plot_vertical_grid(args.infiles,outfile=args.outfile)
