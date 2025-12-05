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
from calc_rms import calc_rms

def calc_rms_filelist(files_in=None, files_in2=None, varnames=None, maskfilename=None, maskname=None,
                      invert_mask=None, file_out=None, file_out2=None, end_file_in=None, end_file_out=None):

    if files_in is None:
        raise Exception("Error : must specify at least two input files.")

    if files_in2 is not None:
        if len(files_in2) != len(files_in):
            raise Exception("Error : second list of input files must be same length as primary list of intput files.")

    ## norms between successive files in the first list
    outmode="w"
    for file1,file2 in [(files_in[i],files_in[i+1]) for i in range(len(files_in)-1)]:
        rms_str = str( calc_rms(filenames=[file1,file2],varnames=varnames,
                                maskfilename=maskfilename,maskname=maskname,invert_mask=invert_mask) )
        if file_out is None:
            print("RMS : "+rms_str)
        else:
            with open(file_out,outmode) as f:
                f.write(rms_str+"\n")
            outmode="a"
                
    ## norms between files in the first and second lists
    outmode="w"
    if files_in2 is not None:
        for file1,file2 in zip(files_in[:],files_in2[:]):
            rms_str = str( calc_rms(filenames=[file1,file2],varnames=varnames,
                                    maskfilename=maskfilename,maskname=maskname,invert_mask=invert_mask) )
            if file_out2 is None:
                print("pairwise RMS : "+rms_str)
            else:
                with open(file_out2,outmode) as f:
                    f.write(rms_str+"\n")
                outmode="a"

    ## norms between files in the first list and an end point
    outmode="w"
    if end_file_in is not None:
        for file1 in files_in:
            rms_str = str( calc_rms(filenames=[file1,end_file_in],varnames=varnames,
                                    maskfilename=maskfilename,maskname=maskname,invert_mask=invert_mask) )
            if end_file_out is None:
                print("RMS w.r.t. endpoint : "+rms_str)
            else:
                with open(end_file_out,outmode) as f:
                    f.write(rms_str+"\n")
                outmode="a"
                    
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
    parser.add_argument("-o", "--file_out", action="store",dest="file_out",
                         help="optional output file")
    parser.add_argument("-p", "--file_out2", action="store",dest="file_out2",
                         help="optional output file for list1-list2 norms")
    parser.add_argument("-e", "--end_file_in", action="store",dest="end_file_in",
                         help="input end file")
    parser.add_argument("-E", "--end_file_out", action="store",dest="end_file_out",
                         help="output end file")

    args = parser.parse_args()

    calc_rms_filelist(files_in=args.files_in,files_in2=args.files_in2,varnames=args.varnames,
                  file_out=args.file_out, file_out2=args.file_out2, end_file_in=args.end_file_in,
                  end_file_out=args.end_file_out, maskfilename=args.maskfilename, maskname=args.maskname,
                  invert_mask=args.invert_mask)
