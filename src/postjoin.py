#! /usr/bin/env python
'''
Rejoin multiple files that were processed in parallel by Spice. 
Possibly with redundant rows/columns that need to be stripped off first.

@author: Dave Storkey
@date: Dec 2021
'''

import xarray as xr
import numpy as np
import re

def postjoin(infiles=None, outfile=None, extra_cols=None, join_dim=None):

    extra_cols_flag=(extra_cols is None)

    if join_dim is None:
        join_dim = 'x'

    datasets_in = []
    for infile in infiles:
        dataset_in = xr.open_dataset(infile)
        if extra_cols_flag:
            # if we haven't explicitly set extra_cols check filename for "Xn.nc" at the end.
            if re.findall(r"cn\.\dX",infile[::-1]):
               extra_cols = int(infile[-4])
            else:
                extra_cols = 0
        if extra_cols > 0:
            xslice = slice(extra_cols,None,None)
            indexer = {join_dim:xslice}
            datasets_in.append(dataset_in.isel(indexers=indexer))
        else:
            datasets_in.append(dataset_in)

    dataset_out = xr.concat(datasets_in,dim=join_dim)
    dataset_out.to_netcdf(outfile)


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", action="store",dest="infiles",nargs="+",
                    help="input files")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",
                    help="output file")
    parser.add_argument("-X", "--extra_cols", action="store",dest="extra_cols",type=int,
                    help="number of extra columns to add to each output file")
    parser.add_argument("-D", "--join_dim", action="store",dest="join_dim",
                    help="dimension to concatenate along - default x")
 
    args = parser.parse_args()

    postjoin(infiles=args.infiles, outfile=args.outfile, 
             extra_cols=args.extra_cols, join_dim=args.join_dim)

