#! /usr/bin/env python

'''
Script to plot timeseries of specific model points extracted by
calc_rms_series.py

@author: Dave Storkey
@date: Apr 2026
'''

import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt

def plot_series(infile=None, outfile=None, deffile=None, zero_mean=None,
                area_list=None, ylabel=None, title=None):

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
            for ii, label in enumerate(labels):
                for area in area_list:
                    if area in label:
                        indices.append(ii)
                        break
            # remove duplicate values:
            indices=list(set(indices))

    dates=[]
    with open(infile,"r") as f:
        flist = list(f)
        points=[[] for _ in range(len(flist[0].split(":")[1].split(",")))]
        if len(indices) == 0:
            indices=[ii for ii in range(len(flist[0].split(":")[1].split(",")))]
        for ii,line in enumerate(flist):
            date_str=line.split(":")[0]
            dates.append(datetime.strptime(date_str,'%Y%m%d'))
            for ii, value in enumerate(line.split(":")[1].split(",")):
                points[ii].append(float(value))

    for ii, tseries in enumerate(points):
        if ii in indices:        
            if zero_mean:
                tseries=np.array(tseries)
                tseries[:]=tseries[:]-np.mean(tseries[:])
            plt.plot(dates,tseries)
        else:
            del_label = labels.pop(ii)
            
    plt.legend(labels)
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
    parser.add_argument("-A", "--area_list", action="store",dest="area_list",nargs="+",
                         help="list of areas to plot based on points definition file.")

    args = parser.parse_args()
    plot_series(infile=args.infile, outfile=args.outfile, deffile=args.deffile, zero_mean=args.zero_mean,
                area_list=args.area_list, ylabel=args.ylabel, title=args.title)
