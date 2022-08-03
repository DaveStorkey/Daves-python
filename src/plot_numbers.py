#! /usr/bin/env python

'''
Plot numbers read from a text file.

Created Nov 2016.

@author: Dave Storkey
'''

import matplotlib
import matplotlib.pyplot as plt

def plot_numbers(infilenames=None, outfilename=None, x_offset=None, scale_factor=None,
                 ymin=None, ymax=None, xtitle=None, ytitle=None, legend=None):

    if x_offset is None:
        x_offset=0.0

    if scale_factor is None:
        scale_factor=1.0

    for infilename in infilenames:

        with open(infilename,'r') as f:
            read_data = f.read()
        float_data=[float(val)*scale_factor for val in read_data.splitlines()]

        if ymin is None:
            ymin = min(val for val in float_data)
        if ymax is None:
            ymax = max(val for val in float_data)

        xaxis = [val+x_offset for val in range(len(float_data))]

        plt.plot(xaxis, float_data)

    plt.gca().set_ylim([ymin,ymax])

    if xtitle is not None:
        plt.gca().set_xlabel(xtitle)
    if ytitle is not None:
        plt.gca().set_ylabel(ytitle)

    if legend is not None:
        plt.legend(legend,loc=1,fontsize='medium')

    if outfilename is not None:
        matplotlib.rcParams['font.size'] =8
        plt.savefig(outfilename,dpi=200)
    else:
        plt.show()
        
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", action="store", dest="infilenames", nargs="+", default=None,
                    help="name of input file")
    parser.add_argument("-o", "--outfile", action="store",dest="outfilename",default=None,
                    help="name of output file")
    parser.add_argument("-O", "--x_offset", action="store",dest="x_offset",type=float,default=None,
                    help="start of y-axis")
    parser.add_argument("-s", "--scale_factor", action="store",dest="scale_factor",type=float,default=None,
                    help="scale factor to multiply data")
    parser.add_argument("-y", "--ymin", action="store",dest="ymin",type=float,default=None,
                    help="start of y-axis")
    parser.add_argument("-Y", "--ymax", action="store",dest="ymax",type=float,default=None,
                    help="end of y-axis")
    parser.add_argument("-p", "--xtitle", action="store",dest="xtitle",default=None,
                    help="title for x-axis")
    parser.add_argument("-q", "--ytitle", action="store",dest="ytitle",default=None,
                    help="title for y-axis")
    parser.add_argument("-L", "--legend", action="store", dest="legend", nargs="+", default=None,
                    help="legend labels (one for each input file)")

    args = parser.parse_args()

    plot_numbers(infilenames=args.infilenames,outfilename=args.outfilename,
                 x_offset=args.x_offset,scale_factor=args.scale_factor,ymin=args.ymin,ymax=args.ymax,
                 xtitle=args.xtitle,ytitle=args.ytitle,legend=args.legend)
