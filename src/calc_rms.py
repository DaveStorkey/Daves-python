#! /usr/bin/env python
# 
# Routine to calculate the RMS difference between two 3D fields.
#
# DS. Nov 2025
#

import netCDF4 as nc
import numpy as np
import numpy.ma as ma

def rms(x):
    rms_x = np.sqrt(np.mean(np.square(x)))
    return rms_x

def calc_rms(filenames=None, varnames=None, outfile=None, maskfilename=None, maskname=None,
             invert_mask=None):
    
    if len(filenames) != 2:
        raise Exception("Error : must specify two filenames.")

    if type(varnames) is not list:
        varnames=[varnames]

    fields=[[],[]]
    for ii, filename in enumerate(filenames):
        with nc.Dataset(filename,'r') as infile:
            for var in varnames:
                fields[ii].append(infile.variables[var][:])

    if maskfilename is not None:
        with nc.Dataset(maskfilename,'r') as maskfile:
            mask = maskfile.variables[maskname][:]
        if invert_mask:
            if type(mask[0]) is np.bool_:
                mask[:] = ~mask[:]
            else:
                mask[:] = 1 - mask[:]

        for fieldlist in fields:
            for field in fieldlist:
                field.mask = mask
    
    # if only one varname supplied this will flatten the 2D array to a 1D array
    # otherwise we end up with a single concatenated 2D array. Either works for
    # calculating the RMS.     
    field1 = ma.concatenate((fields[0]))
    field2 = ma.concatenate((fields[1]))

    field_diff = field2 - field1
    rms_x = ma.sqrt(ma.mean(field_diff*field_diff))

    print("RMS is ",rms_x)
    return(rms_x)

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--filenames", action="store",dest="filenames",nargs=2,
                    help="names of two input file(s)")
    parser.add_argument("-v", "--varnames", action="store",dest="varnames",nargs="+",
                    help="name(s) of variable(s) to plot")
    parser.add_argument("-M", "--maskfilename", action="store",dest="maskfilename",
                    help="name of file containing mask field")
    parser.add_argument("-m", "--maskname", action="store",dest="maskname",
                    help="name of mask field")
    parser.add_argument("-X", "--invert_mask", action="store_true",dest="invert_mask",
                    help="invert the mask field before applying")
 
    args = parser.parse_args()

    calc_rms(filenames=args.filenames, varnames=args.varnames, invert_mask=args.invert_mask,
             maskfilename=args.maskfilename, maskname=args.maskname)
