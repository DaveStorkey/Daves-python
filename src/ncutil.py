#! /usr/bin/env python
"""
    Module to copy a netcdf file with the option to leave the output
    file open for editing. 
"""

__author__ = 'Dave Storkey (dave.storkey@metoffice.gov.uk)'
__date__ = '$Date: 2013/30/7 $'

import sys
import numpy as np
import netCDF4

def copy_ncfile(infile,outfile,close='false',format='NETCDF4',fields=None):

    file_in = netCDF4.Dataset(infile, mode='r')
    file_out = netCDF4.Dataset(outfile, mode='w',format=format)

    print ' '
    print 'Dimensions : '
    for dim in file_in.dimensions.keys():
        dim_in = file_in.dimensions[dim]
        print dim, len(dim_in), dim_in.isunlimited()
        if dim_in.isunlimited():
            file_out.createDimension(dim,size=None)
        else:
            file_out.createDimension(dim,size=len(dim_in))
        
    print ' '
    print 'Variables : '
    for var in file_in.variables.keys():
        var_in = file_in.variables[var]
        print var, var_in.dimensions
        if fields is None or var in fields:
            file_out.createVariable(var,var_in.dtype,dimensions=var_in.dimensions)
            var_out = file_out.variables[var]
            atts = var_in.__dict__
            var_out.setncatts(atts)
            if var_in.ndim == 1:
                var_out[:] = var_in[:]
            elif var_in.ndim == 2:
                var_out[:,:] = var_in[:,:]
            elif var_in.ndim == 3:
                var_out[:,:,:] = var_in[:,:,:]
            elif var_in.ndim == 4:
                var_out[:,:,:,:] = var_in[:,:,:,:]

    if close == 'true':
        file_out.close()
    else:
        return file_out

if __name__ == "__main__":
    copy_ncfile(infile=sys.argv[1],outfile=sys.argv[2],close='true')
