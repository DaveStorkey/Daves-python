#! /usr/bin/env python

'''
Routine to difference two fields that have slightly different
formats, meta-data, with or without halo points etc.

@author: Dave Storkey
@date: Nov 2021
'''

import xarray as xr
import numpy as np

def diff_fields(infiles=None, invars=None, outfile=None, outformat=None):

    if outformat is None:
        outformat = 0

    try:
        iter1 = iter(infiles)
    except TypeError:
        raise Exception("Error: must specify two input files.")
    else:
        if len(infiles) != 2: 
            raise Exception("Error: must specify two input files.")

    try:
        iter1 = iter(invars)
    except TypeError:
        raise Exception("Error: must specify two input variables.")
    else:
        if len(infiles) != 2: 
            raise Exception("Error: must specify two input variables.")

    field = []
    for file, var in zip(infiles,invars):
        with xr.open_dataset(file) as indata:
            field.append(indata.get(var).squeeze())

    print("field[0].shape : ",field[0].shape)
    print("field[1].shape : ",field[1].shape)

    xdiff = field[1].shape[-1] - field[0].shape[-1]
    ydiff = field[1].shape[-2] - field[0].shape[-2]

    if xdiff == 0 and ydiff == 0:
        diff_field = field[0].values - field[1].values                          
    elif xdiff == 2 and ydiff == 2:
        if outformat == 0:
            if len(field[0].shape) == 3:
                diff_field = field[0].values - field[1].values[:,1:-1,1:-1]
            else:
                diff_field = field[0].values - field[1].values[1:-1,1:-1]
        else:
            diff_field = np.empty(field[1].shape)
            diff_field[:] = np.nan
            if len(field[0].shape) == 3:
                diff_field[:,1:-1,1:-1] = field[0].values - field[1].values[:,1:-1,1:-1]
            else:
                diff_field[1:-1,1:-1] = field[0].values - field[1].values[1:-1,1:-1]
    elif xdiff == -2 and ydiff == -2:
        if outformat == 0:
            diff_field = np.empty(field[0].shape)
            diff_field[:] = np.nan
            if len(field[0].shape) == 3:
                diff_field[:,1:-1,1:-1] = field[0].values[:,1:-1,1:-1] - field[1].values
            else:
                diff_field[1:-1,1:-1] = field[0].values[1:-1,1:-1] - field[1].values
        else:
            if len(field[0].shape) == 3:
                diff_field = field[0].values[:,1:-1,1:-1] - field[1].values
            else:
                diff_field = field[0].values[1:-1,1:-1] - field[1].values
    else:
        raise Exception("Error: can only deal with 1-point haloes.")

    # deep copy by default
    outfield = field[outformat].copy(data=diff_field)    

    outdata = outfield.to_dataset()    
    outdata.to_netcdf(outfile)
                
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", action="store",dest="infiles",nargs="+",
                         help="name of file with closea_mask in it")
    parser.add_argument("-v", "--invars", action="store",dest="invars",nargs="+",
                         help="source bathymetry to copy from")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",
                         help="name of output file")
    parser.add_argument("-F", "--outformat", action="store",dest="outformat",type=int,
                         help="which input field format to use for output field (0 or 1, default 0)")

    args = parser.parse_args()

    diff_fields(infiles=args.infiles,invars=args.invars,outfile=args.outfile,
                outformat=args.outformat)

