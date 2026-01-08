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

def get_fields(infile=None, varnames=None, masks=None ):

    if masks is None:
        masks=[None]

    fields_out=[]
    with nc.Dataset(infile,'r') as indata:
        for tname in ["t", "time", "time_counter"]:
            if tname in indata.dimensions.keys():
                tdim=True
                break
        else:
            tdim=False
        for mask in masks:
            for varname in varnames:
                fields_in=[]
                vars_to_read=varname.split("+")
                for var in vars_to_read:
                    if tdim:
                        # take the first record if there's a time dimension
                        fields_in.append(indata.variables[var][0][:])
                    else:
                        fields_in.append(indata.variables[var][:])
                    if mask is not None:
                        fields_in[-1].mask = mask
                if len(fields_in) > 1:
                    fields_out.append( ma.concatenate((fields_in)) )
                else:
                    fields_out.append( fields_in[0] )
    return fields_out


def calc_rms_series(files_in=None, files_in2=None, varnames=None, maskfilename=None, masknames=None,
                    invert_mask=None, end_file_in=None, file_out_stem=None, append=None):

    if files_in is None:
        raise Exception("Error : must specify at least two input files.")

    if files_in2 is None:
        files_in2 = [None]*len(files_in)
    elif len(files_in2) != len(files_in):
        raise Exception("Error : second list of input files must be same length as primary list of intput files.")
        
    if varnames is None:
        raise Exception("Error : must specify at least one variable (varnames)")
    else:
        nvar=len(varnames)
    
    if maskfilename is not None:
        if masknames is None:
            masknames=["tmask"]
        elif type(masknames) is not list:
            masknames=[masknames]
        masks=[]
        with nc.Dataset(maskfilename,'r') as maskfile:
            for maskname in masknames:
                masks.append(maskfile.variables[maskname][:])
                if invert_mask:
                    if type(masks[-1]) is np.bool_:
                        masks[-1][:] = ~masks[-1][:]
                    else:
                        masks[-1][:] = 1 - masks[-1][:]
    else:
        masknames=["global"]
        masks=[None]
                
    if end_file_in is not None:
        endfields = get_fields(infile=end_file_in, varnames=varnames, masks=masks)
                
    if file_out_stem is None:
        file_out_stem="RMS_diffs"
        
    rms_seq=[]
    rms_pairwise=[]
    rms_wrt_endpoint=[]
    fields1_prev=None
    for file1, file2 in zip(files_in, files_in2):
        print("Working on file "+file1)
        fields1 = get_fields(infile=file1, varnames=varnames, masks=masks)
        if fields1_prev is not None:
            fields_diff = [field1-field1_prev for field1,field1_prev in zip(fields1,fields1_prev)]
            rms_seq.append( [ma.sqrt(ma.mean(field_diff*field_diff)) for field_diff in fields_diff] )
        fields1_prev = fields1
        if file2 is not None:
            fields2 = get_fields(infile=file2, varnames=varnames, masks=masks)
            fields_diff = [field2-field1 for field1,field2 in zip(fields1,fields2)]
            rms_pairwise.append( [ma.sqrt(ma.mean(field_diff*field_diff)) for field_diff in fields_diff] )
        if end_file_in is not None:
            fields_diff = [field1-endfield for field1,endfield in zip(fields1,endfields)]
            rms_wrt_endpoint.append( [ma.sqrt(ma.mean(field_diff*field_diff)) for field_diff in fields_diff] )

    if append:
        mode="a"
    else:
        mode="w"

    for ii, maskname in enumerate(masknames):
        range_to_write=slice(ii*nvar,(ii+1)*nvar)
        with open(file_out_stem+"_seq_"+maskname+".dat",mode) as f:
            if mode == "w":
                f.write(",".join([varname for varname in varnames])+"\n")
            for rms_out in rms_seq:
                f.write(",".join([str(rms_write) for rms_write in rms_out[range_to_write]])+"\n")
                
    if files_in2[0] is not None:
        for ii, maskname in enumerate(masknames):
            range_to_write=slice(ii*nvar,(ii+1)*nvar)
            with open(file_out_stem+"_pairwise_"+maskname+".dat",mode) as f:
                if mode == "w":
                    f.write(",".join([varname for varname in varnames])+"\n")
                for rms_out in rms_pairwise:
                    f.write(",".join([str(rms_write) for rms_write in rms_out[range_to_write]])+"\n")
                
    if end_file_in is not None:
        for ii, maskname in enumerate(masknames):
            range_to_write=slice(ii*nvar,(ii+1)*nvar)
            with open(file_out_stem+"_wrt_endpoint_"+maskname+".dat",mode) as f:
                if mode == "w":
                    # for the RMS w.r.t. endpoint write the endpoint filename to the .dat file for reference.
                    f.write(end_file_in+"\n")
                    f.write(",".join([varname for varname in varnames])+"\n")
                for rms_out in rms_wrt_endpoint:
                    f.write(",".join([str(rms_write) for rms_write in rms_out[range_to_write]])+"\n")
                                    
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
    parser.add_argument("-m", "--masknames", action="store",dest="masknames",nargs="+",
                    help="name(s) of mask field")
    parser.add_argument("-X", "--invert_mask", action="store_true",dest="invert_mask",
                    help="invert the mask field before applying")
    parser.add_argument("-o", "--file_out", action="store",dest="file_out_stem",
                         help="filename stem of output file")
    parser.add_argument("-A", "--append", action="store_true",dest="append",
                    help="append data to existing files")
    parser.add_argument("-e", "--end_file_in", action="store",dest="end_file_in",
                         help="input end file")

    args = parser.parse_args()

    calc_rms_series(files_in=args.files_in,files_in2=args.files_in2,varnames=args.varnames,
                    file_out_stem=args.file_out_stem, end_file_in=args.end_file_in, append=args.append,
                    maskfilename=args.maskfilename, masknames=args.masknames, invert_mask=args.invert_mask)
