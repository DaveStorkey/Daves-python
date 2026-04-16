#! /usr/bin/env python

'''
Script to plot timeseries at specific model points extracted by
calc_rms_series.py

@author: Dave Storkey
@date: Apr 2026
'''

import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt

def plot_series(infile=None, outfile=None, deffile=None, zero_mean=None,
                abs_diff=None, area_list=None, exclude=None, ylabel=None, title=None):

    legend=[]
    indices=[]
    if deffile is not None:
        labels=[]
        with open(deffile, "r") as deffile:
            # get rid of any blank lines.
            lines = [line.rstrip() for line in deffile]
            lines = [line for line in lines if line]
            for line in lines:
                print("line : ",line)
                labels.append(line.split("#")[1])
        if area_list is not None:
            indices_select=[]
            for ii, label in enumerate(labels):
                for area in area_list:
                    if area in label:
                        indices_select.append(ii)
                        break
            # remove duplicate values:
            indices_select=list(set(indices_select))

    dates=[]
    with open(infile,"r") as f:
        # f is an iterator. Converting it to a list means the whole thing is
        # stored in memory but makes it easier to loop over it more than once.
        flist = list(f)
        points=[[] for _ in range(len(flist[0].split(":")[1].split(",")))]
        if area_list is None:
            indices=[ii for ii in range(len(flist[0].split(":")[1].split(",")))]
        elif exclude:
            indices=[ii for ii in range(len(flist[0].split(":")[1].split(","))) if ii not in indices_select]
        else:
            indices=indices_select
        for ii,line in enumerate(flist):
            date_str=line.split(":")[0]
            dates.append(datetime.strptime(date_str,'%Y%m%d'))
            for ii, value in enumerate(line.split(":")[1].split(",")):
                points[ii].append(float(value))

    labels_to_plot=[]
    for ii, tseries in enumerate(points):
        if ii in indices:        
            labels_to_plot.append(labels[ii])
            tseries_to_plot=None
            if zero_mean:
                tseries=np.array(tseries)
                tseries_to_plot[:]=tseries[:]-np.mean(tseries[:])
                dates_to_plot=dates
            if abs_diff:
                tseries_roll=np.roll(tseries,1)
                tseries_to_plot=np.abs(tseries[1:]-tseries_roll[1:])
                dates_to_plot=dates[1:]
            if tseries_to_plot is None:
                tseries_to_plot=tseries
                dates_to_plot=dates
            plt.plot(dates_to_plot,tseries_to_plot)
            
    plt.legend(labels_to_plot)
    plt.xlabel("calendar year")
    if ylabel is not None:
        plt.ylabel(ylabel)
    if title is not None:
        plt.title(title)
    plt.savefig(outfile)
    
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", action="store",dest="infile",
                         help="input data file")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",
                         help="output plot file")
    parser.add_argument("-d", "--deffile", action="store",dest="deffile",
                         help="points definitions file for legend labels")
    parser.add_argument("-y", "--ylabel", action="store",dest="ylabel",
                         help="y-axis label - details of field plotted")
    parser.add_argument("-t", "--title", action="store",dest="title",
                         help="title for plot")
    parser.add_argument("-Z", "--zero_mean", action="store_true",dest="zero_mean",
                         help="subtract mean from each timeseries")
    parser.add_argument("-D", "--abs_diff", action="store_true",dest="abs_diff",
                         help="plot absolute sequential differences")
    parser.add_argument("-A", "--area_list", action="store",dest="area_list",nargs="+",
                         help="list of areas to plot based on points definition file.")
    parser.add_argument("-X", "--exclude", action="store_true",dest="exclude",
                         help="exclude selected areas from plot")

    args = parser.parse_args()
    plot_series(infile=args.infile, outfile=args.outfile, deffile=args.deffile, zero_mean=args.zero_mean,
                abs_diff=args.abs_diff,area_list=args.area_list, exclude=args.exclude, ylabel=args.ylabel, title=args.title)
