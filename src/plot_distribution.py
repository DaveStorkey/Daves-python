#! /usr/bin/env python
"""
 Script to blah
"""

__author__ = 'Dave Storkey (dave.storkey@metoffice.gov.uk)'
__date__ = '$Date: Oct 2015 $'

import netCDF4 as nc
import numpy.ma as ma
import matplotlib as mpl
import matplotlib.pyplot as plt

def plot_distribution(datafile,field,minfld=None,maxfld=None,bins=None,rec=None,
                xmin=None,xmax=None,ymin=None,ymax=None,outfile=None,title=None):

    dataset = nc.Dataset(datafile)
    if rec is not None:
        field_in = dataset.variables[field][rec].flatten()
    else:
        field_in = dataset.variables[field][:].flatten()
    field_in = field_in[~field_in.mask]
    if minfld is not None and maxfld is not None:
        field_to_plot = field_in[ma.where((field_in >= minfld) & (field_in <= maxfld))]
    elif minfld is not None:
        field_to_plot = field_in[ma.where(field_in >= minfld)]
    elif maxfld is not None:
        field_to_plot = field_in[ma.where(field_in >= minfld)]
    else:
        field_to_plot = field_in

    if bins is not None:
        plt.hist(field_to_plot,bins=bins)
    else:
        plt.hist(field_to_plot)

    if xmin is not None:
        plt.gca().set_xlim(left=xmin)
    if xmax is not None:
        plt.gca().set_xlim(right=xmax)
    if ymin is not None:
        plt.gca().set_ylim(bottom=ymin)
    if ymax is not None:
        plt.gca().set_ylim(top=ymax)
    
    if title is not None:
        plt.gcf().suptitle(title, fontsize=16, y=0.92)    

    if outfile is not None:
        mpl.rcParams['font.size'] =10
        plt.savefig(outfile,dpi=200)
    else:
        plt.show()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="name of input file")
    parser.add_argument("field", help="name of field to plot")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",default=None,
                    help="plot to output file - filetype by extension. If unset plot to GUI.")
    parser.add_argument("-t", "--title", action="store",dest="title",default=None,
                    help="title for plot")
    parser.add_argument("-b", "--bins", action="store",dest="bins",type=float,nargs='+',default=None,
                    help="minimum field value to plot")
    parser.add_argument("-r", "--rec", action="store",dest="rec",type=int,default=None,
                    help="record number to read in file (counting from zero)")
    parser.add_argument("-f", "--mnfld", action="store",dest="minfld",type=float,default=None,
                    help="minimum field value to plot")
    parser.add_argument("-F", "--mxfld", action="store",dest="maxfld",type=float,default=None,
                    help="maximum field value to plot")
#    parser.add_argument("-l", "--levels", action="store_const",dest="plottype",const="levels",
#                    default="depths",help="depth level to plot")
#    parser.add_argument("-d", "--depths", action="store_const",dest="plottype",const="depths",
#                    default="depths",help="depth level to plot")
#    parser.add_argument("-c", "--cumulative", action="store_const",dest="plottype",const="cumulative",
#                    default="depths",help="depth level to plot")
#    parser.add_argument("-D", "--maxdepth", action="store",dest="maxdepth",type=int,default=None,
#                    help="time level to plot")
    parser.add_argument("-x", "--xmin", action="store",dest="xmin",type=float,default=None,
                    help="start of x-axis")
    parser.add_argument("-X", "--xmax", action="store",dest="xmax",type=float,default=None,
                    help="end of x-axis")
    parser.add_argument("-y", "--ymin", action="store",dest="ymin",type=float,default=None,
                    help="start of y-axis")
    parser.add_argument("-Y", "--ymax", action="store",dest="ymax",type=float,default=None,
                    help="end of y-axis")
    args = parser.parse_args()

    plot_distribution(args.filename,args.field,bins=args.bins,minfld=args.minfld,maxfld=args.maxfld,
              xmin=args.xmin,xmax=args.xmax,ymin=args.ymin,ymax=args.ymax,outfile=args.outfile,title=args.title,
              rec=args.rec)        

