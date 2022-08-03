#! /usr/bin/env python
"""
    Script to extend fields on standard ORCA025 grid for extended ORCA025 grid.
    Fill in extra area with zeroes or FillValues.
"""

__author__ = 'Dave Storkey (dave.storkey@metoffice.gov.uk)'
__date__ = '$Date: June 2014 $'

import sys
import numpy as np
import netCDF4
from ncutil import copy_ncfile as copy

def extended_orca_input_fields(input_file,coords_file,output_file,fields=None,extrap=None):

    dataset_in = netCDF4.Dataset(input_file, mode='r')
    dims_in = dataset_in.dimensions
    vars_in = dataset_in.variables
    
    dataset_out = copy(coords_file,output_file,fields=['nav_lat','nav_lon'])

    ysize_in = len(dims_in['y'])
    ysize_out = len(dataset_out.dimensions['y'])   

    if fields is None:
        # in this case try to extend all 2D, 3D or 4D variables in file   
        fields = vars_in.keys()

    for dim in dims_in.keys():
        if dim not in dataset_out.dimensions.keys():
            if dims_in[dim].isunlimited():
                dataset_out.createDimension(dim,size=None)
            else:
                dataset_out.createDimension(dim,size=len(dims_in[dim]))

    for var in vars_in.keys():
        dataset_out.createVariable(var,vars_in[var].dtype,dimensions=vars_in[var].dimensions)
        var_out = dataset_out.variables[var]
        # Need to define the attributes before we fill in the data. The netcdf libraries
        # don't like it if you try to define the _FillValue *after* filling the data arrays.
        atts = vars_in[var].__dict__
        var_out.setncatts(atts)    
        if len(var_out.shape) == 1:
           var_out[:] = vars_in[var][:]
        elif len(var_out.shape) == 2 and var in fields:
           var_out[-ysize_in:,:] = vars_in[var][:,:]
           var_out[0:-ysize_in,:] = 0.0
           if extrap is not None:
               for j in range(-ysize_in-extrap,-ysize_in):
                   var_out[j,:] = var_out[-ysize_in+1,:]
        elif len(var_out.shape) == 3 and var in fields:
           var_out[:,-ysize_in:,:] = vars_in[var][:,:,:]
           var_out[:,0:-ysize_in,:] = 0.0
           if extrap is not None:
               for j in range(-ysize_in-extrap,-ysize_in):
                   var_out[:,j,:] = var_out[:,-ysize_in+1,:]
        elif len(var_out.shape) == 4 and var in fields:
           var_out[:,:,-ysize_in:,:] = vars_in[var][:,:,:,:]
           var_out[:,:,0:-ysize_in,:] = 0.0
           if extrap is not None:
               for j in range(-ysize_in-extrap,-ysize_in):
                   var_out[:,:,j,:] = var_out[:,:,-ysize_in+1,:]

    dataset_in.close()
    dataset_out.close()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="name of input file")
    parser.add_argument("coords_file", help="name of coordinates file of target grid")
    parser.add_argument("output_file", help="name of output file")
    parser.add_argument("-f", "--fields", action="store",nargs='?',dest="fields",
                    default=None,help="fields to extend")
    parser.add_argument("-x", "--extrap", action="store",dest="extrap",type=int,
                    default=None,help="number of rows to fill with extrapolated values")

    args = parser.parse_args()

    extended_orca_input_fields(args.input_file, args.coords_file, args.output_file, 
                               fields=args.fields, extrap=args.extrap)
                        
