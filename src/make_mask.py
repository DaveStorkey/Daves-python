#! /usr/bin/env python

'''
Script to generate a 2D mask based on the properties of
an input field, eg. mask based on a maximum depth of a 
bathymetry. 

@author: Dave Storkey
@date: May 2022
'''

import xarray as xr
import numpy as np
import re

def fill_field(field_in, fill_value):

    field_out = np.where( (np.roll(field_in, 1,axis=0) == fill_value) | \
                          (np.roll(field_in,-1,axis=0) == fill_value) | \
                          (np.roll(field_in,-1,axis=1) == fill_value) | \
                          (np.roll(field_in,-1,axis=1) == fill_value), fill_value,field_in)

    return field_out

def make_mask(files_in=None, file_out=None, varnames=None, 
              min_thresholds=None, max_thresholds=None,
              north=None, south=None, east=None, west=None,
              fill_antarc_shelf=None):

    if len(files_in) != len(varnames):
        if len(files_in) == 1:
            files_in = [files_in]*len(varnames)
        else:
            raise Exception("Error: number of input files must equal number of varnames (or be 1)")

    if min_thresholds is not None:
        if len(min_thresholds) != len(varnames):
            raise Exception("Error: number of min_thresholds must equal number of varnames")
    else:
        min_thresholds = [None]*len(varnames)

    if max_thresholds is not None:
        if len(max_thresholds) != len(varnames):
            raise Exception("Error: number of max_thresholds must equal number of varnames")
    else:
        max_thresholds = [None]*len(varnames)

    with xr.open_dataset(files_in[0]) as indata:
        try:
            nav_lat = indata.nav_lat
        except(AttributeError):
            nav_lat = None
            if north is not None or south is not None:
                raise Exception("Error: can't find nav_lat in input file.")
        try:
            nav_lon = indata.nav_lon
        except(AttributeError):
            nav_lon = None
            if east is not None or west is not None:
                raise Exception("Error: can't find nav_lon in input file.")

    fields=[]
    for file_in, varname in zip(files_in,varnames):
        with xr.open_dataset(file_in) as indata:
            fields.append(getattr(indata,varname).squeeze())

    mask = fields[0].copy()
    mask[:] = 1.0
    mask.name = varnames[0]+"_mask"

    for field, min_threshold, max_threshold in zip(fields,min_thresholds,max_thresholds):
        if min_threshold is not None:
            mask = mask.where(field > min_threshold, other=0.0)
        if max_threshold is not None:
            mask = mask.where(field < max_threshold, other=0.0)
    if north is not None:
        mask = mask.where(nav_lat < north, other=0.0)
    if south is not None:
        mask = mask.where(nav_lat > south, other=0.0)
    if east is not None:
        mask = mask.where(nav_lon < east, other=0.0)
    if west is not None:
        mask = mask.where(nav_lon > west, other=0.0)

    if fill_antarc_shelf:
        field.values[0] = -1
        field.values = fill_field(field,-1)
        for jj in range(mask.shape[0]-1,1,-1):
            mask.values[jj-1,:] = np.maximum(mask.values[jj,:],mask.values[jj-1,:]) \
                                * np.sign(field.values[jj,:])       
            
    outdata = mask.to_dataset()    
    if nav_lat is not None and nav_lon is not None:
        outdata.update({'nav_lat':nav_lon ,
                        'nav_lat':nav_lon })

    outdata.to_netcdf(file_out)
                
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--files_in", action="store",dest="files_in",nargs="+",
                         help="input file(s)")
    parser.add_argument("-o", "--file_out", action="store",dest="file_out",
                         help="output file")
    parser.add_argument("-v", "--varnames", action="store",dest="varnames",nargs="+",
                         help="name of variable(s) to edit")
    parser.add_argument("-m", "--min_thresholds", action="store",dest="min_thresholds",type=float,nargs="+",
                         help="minimum threshold(s) for mask")
    parser.add_argument("-M", "--max_thresholds", action="store",dest="max_thresholds",type=float,nargs="+",
                         help="maximum threshold(s) for mask")
    parser.add_argument("-N", "--north", action="store",dest="north",type=float,
                         help="northern latitude limit for mask")
    parser.add_argument("-S", "--south", action="store",dest="south",type=float,
                         help="southern latitude limit for mask")
    parser.add_argument("-E", "--east", action="store",dest="east",type=float,
                         help="eastern longitude limit for mask")
    parser.add_argument("-W", "--west", action="store",dest="west",type=float,
                         help="western longitude limit for mask")
    parser.add_argument("-F", "--fill_shelf", action="store_true",dest="fill_antarc_shelf",
                         help="special function to help create a mask for Antarctic shelf")

    args = parser.parse_args()

    make_mask(files_in=args.files_in,file_out=args.file_out,varnames=args.varnames,
              min_thresholds=args.min_thresholds, max_thresholds=args.max_thresholds,
              north=args.north,south=args.south,east=args.east,west=args.west,
              fill_antarc_shelf=args.fill_antarc_shelf )

