#! /usr/bin/env python

'''
Routine to compare the land sea mask between two
bathymetries. Takes land as zero and sea as nonzero
so works for tmask as well as Bathymetry fields. 

@author: Dave Storkey
@date: Nov 2021
'''

import xarray as xr
import numpy as np

def lsm_diff(infiles=None,invar=None,outfile=None):

    if len(infiles) != 2:
        raise Exception("Error: must specify two input files.")

    if invar is None:
        invar="Bathymetry"

    bathy=[]
    for infile in infiles:
        bathy_in = xr.open_dataset(infile)
        bathy.append(bathy_in.get(invar))
        # Replace masked (land) values with 0 and sea values with 1.
        bathy[-1] = xr.where(np.isnan(bathy[-1].values),0.0,bathy[-1])
        bathy[-1] = xr.where(bathy[-1] != 0.0, 1.0, bathy[-1])
#    lats = bathy[-1].lat.values
#    lons = bathy[-1].lon.values
    
    lsm_cf = xr.where(bathy[1] == bathy[0],bathy[0],5.0)
    lsm_cf_count = np.count_nonzero(lsm_cf.values == 5.0)
    print('Number of non-matching locations : ',lsm_cf_count)

    out_dataset = xr.Dataset(
        {
            "landsea_cf" : lsm_cf
        }
    )
    
    out_dataset.to_netcdf(outfile)

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", action="store",dest="infiles",nargs="+",
                         help="names of two input files")
    parser.add_argument("-v", "--invar", action="store",dest="invar",
                         help="name of input variable")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",
                         help="name of output file")

    args = parser.parse_args()
 
    lsm_diff(infiles=args.infiles,invar=args.invar,outfile=args.outfile)
