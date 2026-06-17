#! /usr/bin/env python

'''
Script to extract residual norm values from "aa_restart.h5" files
from an Anderson acceleration calculation.

@author: Dave Storkey
@date: June 2026
'''

import h5py
import numpy as np

def extract_res_norm(files_in=None, file_out=None, group=None, var=None,nx=None):

    if files_in is None:
        raise Exception("Error : must specify at least one input file.")
        
    if group is None:
        group="x"
        
    if var is None:
        var="res_norm"

    dates=[]
    values=[]
    for file in files_in:
        if "_out" in file:
            dates.append(file[15:23])
        else:
            dates.append(file[11:19])
        with h5py.File(file) as f:
            if group not in f.keys():
                raise Exception("Error: could not find group "+group+" in file.")
            elif var not in f[group].keys():
                raise Exception("Error: could not find var "+var+" in group "+group+".")
            else:
                values.append(f[group][var][()])

    if nx is not None:
        print("Converting values to RMS. NX is "+str(nx))
        values_rms=[]
        for value in values:
            values_rms.append(np.sqrt((value*value)/float(nx)))
        values=values_rms
            
    with open(file_out,"w") as f:
        for date,value in zip(dates,values):
            f.write(date+f" : {value}\n")

            
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--files_in", action="store",dest="files_in",nargs="+",
                         help="list of input files")
    parser.add_argument("-v", "--var", action="store",dest="var",
                         help="name of field to extract")
    parser.add_argument("-o", "--file_out", action="store",dest="file_out",
                         help="name of output file")
    parser.add_argument("-N", action="store",dest="nx",type=int,
                         help="size of model state vector. If set convert from sqrt(x.x) to RMS.")

    args = parser.parse_args()

    extract_res_norm(files_in=args.files_in,var=args.var,file_out=args.file_out,nx=args.nx)
