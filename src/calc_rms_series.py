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
import csv

def decomment(csvfile):
    for row in csvfile:
        raw = row.split('#')[0].strip()
        if raw: yield raw

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
                    invert_mask=None, end_files_in=None, points_file=None, file_out_stem=None, append=None):

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
                
    if end_files_in is not None:
        if not isinstance(end_files_in,list):
            end_files_in=[end_files_in]
        endfields_list=[]
        for end_file in end_files_in:
            endfields_list.append(get_fields(infile=end_file, varnames=varnames, masks=masks))

    if points_file is not None:
        with open(points_file, 'r') as f:
            # turn the iterator into a list to make it easy to iterate over it more than once.
            lines=list(csv.reader(decomment(f), delimiter=','))
            # NB. for advanced indexing of numpy arrays the list of indices
            #     has to be a tuple, *not* a list or a numpy array.
            points_tuple=tuple([[] for _ in range(len(lines[0]))])
            for line in lines:
                print('line : ',line)
                for ii, idx in enumerate(line):
                    print('ii, idx : ',ii,idx)
                    points_tuple[ii].append(int(idx))
        print('points_tuple: ',points_tuple)
        if len(points_tuple[0]) == 0:
            raise Exception('Could not read points file '+points_file)
        
    if file_out_stem is None:
        file_out_stem="RMS_diffs"
        
    dates=[]
    rms_seq=[]
    point_values={}
    for varname in varnames:
        point_values[varname]=[]
    rms_pairwise=[]
    if end_files_in is not None:
        rms_wrt_endpoints=[[] for _ in range(len(end_files_in))]
    fields1_prev=None
    # NB. We always assume that the first file in the input list is a "prev" field which only gets
    #     used to calculate sequential residuals, *not* the pairwise RMS or the RMS w.r.t. endpoint.
    #     This facilitates the usual mode of operation where a lot of files are restored from MASS
    #     in chunks and this script called iteratively for each chunk from calc_rms_massget.sh.
    #     To keep things simple we assume the same number of files in the two file lists, the first
    #     file in the second file list being ignored.
    fields1_prev = get_fields(infile=files_in[0], varnames=varnames, masks=masks)
    for file1, file2 in zip(files_in[1:], files_in2[1:]):
        print("Working on file "+file1)
        # assuming restart file of form RUNID_DATE_...
        date1 = file1.split("_")[1]
        dates.append(date1)
        fields1 = get_fields(infile=file1, varnames=varnames, masks=masks)
        if points_file is not None:
            # just pick out the global mask - assume first in the list
            for varname,field1 in zip(varnames,fields1[:len(varnames)]):
                point_values[varname].append(field1[points_tuple])
        if fields1_prev is not None:
            fields_diff = [field1-field1_prev for field1,field1_prev in zip(fields1,fields1_prev)]
            rms_seq.append( [ma.sqrt(ma.mean(field_diff*field_diff)) for field_diff in fields_diff] )
        fields1_prev = fields1
        if file2 is not None:
            date2 = file2.split("_")[1]
            if date2 != date1:
                raise Exception("Error : dates in two filelists do not match.")
            fields2 = get_fields(infile=file2, varnames=varnames, masks=masks)
            fields_diff = [field2-field1 for field1,field2 in zip(fields1,fields2)]
            rms_pairwise.append( [ma.sqrt(ma.mean(field_diff*field_diff)) for field_diff in fields_diff] )
        if end_files_in is not None:
            for ii, endfields in enumerate(endfields_list):
                fields_diff = [field1-endfield for field1,endfield in zip(fields1,endfields)]
                rms_wrt_endpoints[ii].append( [ma.sqrt(ma.mean(field_diff*field_diff)) for field_diff in fields_diff] )
            
    if append:
        mode="a"
    else:
        mode="w"

    for ii, maskname in enumerate(masknames):
        range_to_write=slice(ii*nvar,(ii+1)*nvar)
        with open(file_out_stem+"_seq_"+maskname+".dat",mode) as f:
            if mode == "w":
                f.write(",".join([varname for varname in varnames])+"\n")
            for date1, rms_out in zip(dates, rms_seq):
                f.write(str(date1)+":"+",".join([str(rms_write) for rms_write in rms_out[range_to_write]])+"\n")

    if points_file is not None:            
        for varname in varnames:
            with open(file_out_stem+"_"+varname+"_points.dat",mode) as f:
                # point_value is a list of values for the points specified
                # for a particular variable at a particular time.
                for date1, point_value in zip(dates,point_values[varname]):
                    f.write(str(date1)+":"+",".join([str(value_out) for value_out in point_value])+"\n")
        
    if files_in2[0] is not None:
        for ii, maskname in enumerate(masknames):
            range_to_write=slice(ii*nvar,(ii+1)*nvar)
            with open(file_out_stem+"_pairwise_"+maskname+".dat",mode) as f:
                if mode == "w":
                    f.write(",".join([varname for varname in varnames])+"\n")
                for date1, rms_out in zip(dates, rms_pairwise):
                    f.write(str(date1)+":"+",".join([str(rms_write) for rms_write in rms_out[range_to_write]])+"\n")
                
    if end_files_in is not None:
        for end_file_in, rms_wrt_endpoint in zip(end_files_in, rms_wrt_endpoints):
            for ii, maskname in enumerate(masknames):
                range_to_write=slice(ii*nvar,(ii+1)*nvar)
                with open(file_out_stem+"_wrt_"+end_file_in.replace(".nc","")+"_"+maskname+".dat",mode) as f:
                    if mode == "w":
                        # for the RMS w.r.t. endpoint write the endpoint filename to the .dat file for reference.
                        # f.write(end_file_in+"\n")
                        f.write(",".join([varname for varname in varnames])+"\n")
                    for date1, rms_out in zip(dates, rms_wrt_endpoint):
                        f.write(str(date1)+":"+",".join([str(rms_write) for rms_write in rms_out[range_to_write]])+"\n")
                                    
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
    parser.add_argument("-e", "--end_files_in", action="store",dest="end_files_in",nargs="*",
                         help="input endpoint files")
    parser.add_argument("-p", "--points_files", action="store",dest="points_file",
                         help="file with list of points to be sampled")

    args = parser.parse_args()

    calc_rms_series(files_in=args.files_in,files_in2=args.files_in2,varnames=args.varnames,
                    file_out_stem=args.file_out_stem, end_files_in=args.end_files_in, append=args.append,
                    maskfilename=args.maskfilename, masknames=args.masknames, invert_mask=args.invert_mask,
                    points_file=args.points_file)
