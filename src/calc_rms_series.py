#! /usr/bin/env python

'''
Script to calculate RMS of arbitrary variables between:

             1. successive files in a timeseries of files;
             2. pairs of files from two timeseries of files;
             3. files in a timeseries of files and an "endpoint" file. 

This is modelled on calc_res_norm.py in my AA_spin repos, but it doesn't
use Samar's machinery for handling the NEMO fields - just uses standard
masking and numpy.ma functionality.

@author: Dave Storkey
@date: Dec 2025
'''

import netCDF4 as nc
import numpy as np
import numpy.ma as ma

def get_fields(infile=None, varnames=None, mask=None ):

    fields_out=[]
    with nc.Dataset(infile,'r') as indata:
        for tname in ["t", "time", "time_counter"]:
            if tname in indata.dimensions.keys():
                tdim=True
                break
        else:
            tdim=False
        for var in varnames:
            if tdim:
                # take the first record if there's a time dimension
                fields_out.append(indata.variables[var][0][:])
            else:
                fields_out.append(indata.variables[var][:])
            if mask is not None:
                fields_out[-1].mask = mask

    # if only one varname supplied ma.concatenate will flatten the 2D array
    # to a 1D array otherwise we end up with a single concatenated 2D array.
    # Either works for calculating the RMS.     
    return ma.concatenate((fields_out))


def calc_rms_series(files_in=None, files_in2=None, varnames=None, maskfilename=None, maskname=None,
                    invert_mask=None, end_file_in=None, file_out_stem=None):

    if files_in is None:
        raise Exception("Error : must specify at least two input files.")

    if files_in2 is None:
        files_in2 = [None]*len(files_in)
    elif len(files_in2) != len(files_in):
        raise Exception("Error : second list of input files must be same length as primary list of intput files.")
        
    if varnames is None:
        raise Exception("Error : must specify at least one variable (varnames)")
        
    mask=None
    if maskfilename is not None:
        if maskname is None:
            maskname="tmask"
        with nc.Dataset(maskfilename,'r') as maskfile:
            mask = maskfile.variables[maskname][:]
        if invert_mask:
            if type(mask[0]) is np.bool_:
                mask[:] = ~mask[:]
            else:
                mask[:] = 1 - mask[:]

    if end_file_in is not None:
        endfield = get_fields(infile=end_file_in, varnames=varnames, mask=mask)
                
    if file_out_stem is None:
        file_out_stem="RMS_diffs"
        
    rms_seq=[]
    rms_pairwise=[]
    rms_wrt_endpoint=[]
    field1_prev=None
    for file1, file2 in zip(files_in, files_in2):
        print("Working on file "+file1)
        field1 = get_fields(infile=file1, varnames=varnames, mask=mask)
        if field1_prev is not None:
            field_diff = field1 - field1_prev
            rms_seq.append( ma.sqrt(ma.mean(field_diff*field_diff)) )
        field1_prev = field1
        if file2 is not None:
            field2 = get_fields(infile=file2, varnames=varnames, mask=mask)
            field_diff = field2 - field1
            rms_pairwise.append( ma.sqrt(ma.mean(field_diff*field_diff)) )
        if end_file_in is not None:
            field_diff = field1 - endfield
            rms_wrt_endpoint.append( ma.sqrt(ma.mean(field_diff*field_diff)) )

    with open(file_out_stem+"_seq.dat","w") as f:
        for rms_out in rms_seq:
            print("rms_out : ",rms_out)
            f.write(str(rms_out)+"\n")
                
    if files_in2[0] is not None:
        with open(file_out_stem+"_pairwise.dat","w") as f:
            for rms_out in rms_pairwise:
                f.write(str(rms_out)+"\n")
                
    if end_file_in is not None:
        with open(file_out_stem+"_wrt_endpoint.dat","w") as f:
            for rms_out in rms_wrt_endpoint:
                f.write(str(rms_out)+"\n")
                                    
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--files_in", action="store",dest="files_in",nargs="+",
                         help="at least two input files")
    parser.add_argument("-j", "--files_in2", action="store",dest="files_in2",nargs="+",
                         help="optional second list of input files")
    parser.add_argument("-v", "--varnames", action="store",dest="varnames",nargs="+",
                         help="name of field(s) to use")
    parser.add_argument("-M", "--maskfilename", action="store",dest="maskfilename",
                    help="name of file containing mask field")
    parser.add_argument("-m", "--maskname", action="store",dest="maskname",
                    help="name of mask field")
    parser.add_argument("-X", "--invert_mask", action="store_true",dest="invert_mask",
                    help="invert the mask field before applying")
    parser.add_argument("-o", "--file_out", action="store",dest="file_out_stem",
                         help="filename stem of output file")
    parser.add_argument("-e", "--end_file_in", action="store",dest="end_file_in",
                         help="input end file")

    args = parser.parse_args()

    calc_rms_series(files_in=args.files_in,files_in2=args.files_in2,varnames=args.varnames,
                  file_out_stem=args.file_out_stem, end_file_in=args.end_file_in,
                  maskfilename=args.maskfilename, maskname=args.maskname, invert_mask=args.invert_mask)
