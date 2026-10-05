#! /usr/bin/env python

'''
Read in closea masks and use them to copy the bathymetry 
for inland seas from one bathymetry file to another. 

Oct 2026 : Update for new treatment of closea mask fields (from NEMO 4.2). DS.

@author: Dave Storkey
@date: Nov 2021
'''

import xarray as xr
import numpy as np

def closea_copy_bathy(bathy_source=None, bathy_target=None, bathy_out=None, 
                      closea_file=None, closea_mask=None, mask_csundef=None,
                      mask_csglo=None, mask_csemp=None, mask_csrnf=None, fill=None):

    with xr.open_dataset(closea_file) as closea_data:
        mask_field={}
        mask_indices={}
        for mask, mask_name in zip([closea_mask, mask_csundef, mask_csglo, mask_csemp, mask_csrnf],
                                   ['closea_mask', 'mask_csundef', 'mask_csglo', 'mask_csemp', 'mask_csrnf']):
            if mask is not None:
                mask_field[mask_name] = getattr(closea_data,mask_name).squeeze()
                # nan_to_num converts NaNs to zeroes
                mask_indices_all = np.unique(np.nan_to_num(mask_field[mask_name].values.astype(int)))
                if mask[0] == -1:
                    # index.item() converts from numpy.int64 to native python "int" type
                    mask_indices[mask_name] = [index.item() for index in mask_indices_all]
                else:
                    if set(mask) <= set(mask_indices_all) :
                        mask_indices[mask_name] = mask
                    else:
                        raise Exception("Error : list of indices does not match field for "+mask_name)
                if 0 in mask_indices[mask_name]:
                    # in the old set up closea_mask=0 was the global ocean minus the closed seas.
                    mask_indices[mask_name].remove(0)
            else:
                mask_field[mask_name] = None
                mask_indices[mask_name] = None
                    
    if not fill:
        with xr.open_dataset(bathy_source) as source_data:
            bathy_source = source_data.Bathymetry

    with xr.open_dataset(bathy_target) as target_data:
        coords={}
        for coordname in ['nav_lat','nav_lon']:
            try:
                coords[coordname] = getattr(target_data,coordname)
            except(AttributeError):
                pass
        bathy_target = target_data.Bathymetry

    for mask_name in ['closea_mask', 'mask_csundef', 'mask_csglo', 'mask_csemp', 'mask_csrnf']:
        index_list=mask_indices[mask_name]
        if index_list is not None:
            print('Processing field : ',mask_name)
            for lake_index in index_list:
                if fill:
                    print('Filling lake index : ',lake_index)
                    print('Number of lake points : ',np.count_nonzero(mask_field[mask_name].astype(int) == lake_index) )
                    bathy_target.values = np.where(mask_field[mask_name].astype(int) == lake_index, 0.0, bathy_target.values)
                else:
                    print('Copying lake index : ',lake_index)
                    print('Number of lake points : ',np.count_nonzero(mask_field[mask_name].astype(int) == lake_index) )
                    bathy_target.values = np.where(mask_field[mask_name].astype(int) == lake_index, bathy_source.values, bathy_target.values)

    outdata = bathy_target.to_dataset()    
    if len(coords.keys()) > 0:
        for key in coords.keys():
            outdata[key] = coords[key]
    outdata.to_netcdf(bathy_out)
                
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-C", "--closea_file", action="store",dest="closea_file",
                         help="name of file with closea_mask in it")
    parser.add_argument("-S", "--bathy_source", action="store",dest="bathy_source",
                         help="source bathymetry to copy from")
    parser.add_argument("-T", "--bathy_target", action="store",dest="bathy_target",
                         help="source bathymetry to copy to")
    parser.add_argument("-o", "--outfile", action="store",dest="bathy_out",
                         help="name of output file")
    parser.add_argument("-F", "--fill", action="store_true",dest="fill",
                         help="fill the specified lakes rather than copying them over from the source bathymetry")
    parser.add_argument("--closea_mask", action="store",dest="closea_mask",type=int,nargs="+",
                         help="List of integer labels for closea_mask to be copied. (Set to -1 to copy all).")
    parser.add_argument("--mask_csundef", action="store",dest="mask_csundef",type=int,nargs="+",
                         help="List of integer labels for mask_csundef to be copied. (Set to -1 to copy all).")
    parser.add_argument("--mask_csglo", action="store",dest="mask_csglo",type=int,nargs="+",
                         help="List of integer labels for mask_csglo to be copied. (Set to -1 to copy all).")
    parser.add_argument("--mask_csemp", action="store",dest="mask_csemp",type=int,nargs="+",
                         help="List of integer labels for mask_csemp to be copied. (Set to -1 to copy all).")
    parser.add_argument("--mask_csrnf", action="store",dest="mask_csrnf",type=int,nargs="+",
                         help="List of integer labels for mask_csrnf to be copied. (Set to -1 to copy all).")

    args = parser.parse_args()

    closea_copy_bathy(closea_file=args.closea_file,bathy_source=args.bathy_source,bathy_target=args.bathy_target,
                      bathy_out=args.bathy_out, fill=args.fill, closea_mask=args.closea_mask, mask_csundef=args.mask_csundef,
                      mask_csglo=args.mask_csglo, mask_csemp=args.mask_csemp, mask_csrnf=args.mask_csrnf )
