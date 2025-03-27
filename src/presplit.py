#! /usr/bin/env python
'''
Split up input files for some processing task so that it can be submitted
to multiple Spice processors in parallel. Splits along the "x" dimension
by default.

@author: Dave Storkey
@date: Dec 2021
'''

import xarray as xr
import numpy as np

def presplit(infile=None, outfile=None, nsplit=None, extra_cols=None, split_dim=None):

    if nsplit is None:
        nsplit = 10

    if split_dim is None:
        split_dim = 'x'

    if extra_cols is None:
        extra_cols = 0

    dataset_in = xr.open_dataset(infile)
    dimlen = dataset_in.dims[split_dim]
    splitlens = np.zeros(nsplit).astype(int)
    splitlens[:] = int(dimlen/nsplit)
    splitlens[:dimlen%nsplit] += 1

    start = 0
    end = 0
    for splitlen in splitlens:
        end += splitlen
        if start == 0:
            xslice = slice(start,end)
        else:
            xslice = slice(start-extra_cols,end)
        indexer = {split_dim:xslice}
        dataset_out = dataset_in.isel(indexers=indexer)
        endm1=end-1
        if start == 0 or extra_cols == 0:
            dataset_out.to_netcdf(f"{outfile}_{start:04}_{endm1:04}.nc")
        else:
            dataset_out.to_netcdf(f"{outfile}_{start:04}_{endm1:04}_X{extra_cols}.nc")
        start += splitlen

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", action="store",dest="infile",
                    help="input file")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",
                    help="output file stem")
    parser.add_argument("-N", "--nsplit", action="store",dest="nsplit",type=int,
                    help="number of target files to split into (default 10)")
    parser.add_argument("-X", "--extra_cols", action="store",dest="extra_cols",type=int,
                    help="number of extra columns to add to each output file")
    parser.add_argument("-D", "--split_dim", action="store",dest="split_dim",
                    help="dimension to split along - default x")
 
    args = parser.parse_args()

    presplit(infile=args.infile, outfile=args.outfile, nsplit=args.nsplit, 
             extra_cols=args.extra_cols, split_dim=args.split_dim)

