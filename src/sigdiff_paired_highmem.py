#! /usr/bin/env python

'''
Routine to calculate the "signficant" part of a time-mean difference
field given the time series of the differences or the two time series
of the fields to be differenced. Insignificant differences are set
to zero in the output field. 

Significance can be measured using a Student's t-test or by comparison
to the standard deviation of one of the fields differenced. 

@author: Dave Storkey
@date: Jan 2017
'''

import netCDF4 as nc
import numpy as np
import numpy.ma as ma
from scipy import stats

def sigdiff(filenames1_in=None,filenames2_in=None,differences_in=None,filename_out=None,
            var_list=None,siglevel=None,nstddev=None):

    #### A few consistency checks ####

    if nstddev is None and siglevel is None:
        raise Exception('Must specify one of nstddev and siglevel')

    if nstddev is not None and siglevel is not None:
        raise Exception("Can't specify both of nstddev and siglevel")

    if nstddev is not None and filenames1_in is None:
        raise Exception("Require infiles1 to be set if nstddev is set.")

    #### Open input files #### 

    need_to_make_diffs = True

    if filenames1_in is not None:

        if len(filenames1_in) == 1:
            # only one file:
            infile1 = nc.Dataset(filenames1_in[0])
        else:
            infile1 = nc.MFDataset(filenames1_in)

        if var_list is None:
            var_list = infile1.variables.keys()

        template_file = infile1

    if filenames2_in is not None:

        if len(filenames2_in) == 1:
            # only one file:
            infile2 = nc.Dataset(filenames2_in[0])
        else:
            infile2 = nc.MFDataset(filenames2_in)

    if differences_in is not None:
        
        need_to_make_diffs = False

        if len(differences_in) == 1:
            # only one file:
            diff_file = nc.Dataset(differences_in[0])
        else:
            diff_file = nc.MFDataset(differences_in)

        if var_list is None:
            var_list = diff_file.variables.keys()

        template_file = diff_file

    #### Open output file and copy relevant dimensions  ####
    #### and coordinate variables from input file       ####

    outfile = nc.Dataset(filename_out,'w')

    dims = template_file.dimensions
    for dim in dims.keys():
        if dim.find('time') == -1 and dim != 't':
            outfile.createDimension(dim, len(dims[dim]))

    for latname in 'lat','latitude','nav_lat':
        try:
            lat = template_file.variables[latname]
        except KeyError:
            pass
        else:
            outfile.createVariable(latname,datatype='f',dimensions=lat.dimensions)
            outfile.variables[latname][:] = lat[:]
            if latname in var_list:
                var_list.remove(latname)
            break

    for lonname in 'lon','longitude','nav_lon':
        try:
            lon = template_file.variables[lonname]
        except KeyError:
            pass
        else:
            outfile.createVariable(lonname,datatype='f',dimensions=lon.dimensions)
            outfile.variables[lonname][:] = lon[:]
            if lonname in var_list:
                var_list.remove(lonname)
            break

    for depthname in 'depth','deptht','depthu','depthv','depthw':
        try:
            depth = template_file.variables[depthname]
        except KeyError:
            pass
        else:
            outfile.createVariable(depthname,datatype='f',dimensions=depth.dimensions)
            outfile.variables[depthname][:] = depth[:]
            if depthname in var_list:
                var_list.remove(depthname)
            break

    #### Loop over variable names in the list and calculate significant differences ####

    for varname in var_list:

        if need_to_make_diffs:
            var1 = infile1.variables[varname][:]
            var2 = infile2.variables[varname][:]
            var = var2 - var1
        else:
            if filenames1_in is not None:
                var1 = infile1.variables[varname][:]
            var = diff_file.variables[varname][:]

        var_template = template_file.variables[varname]
        dim_count=-1
        dim_to_mean_number=-1
        dims_to_write=[]
        for dim in var_template.dimensions:
            dim_count = dim_count+1
            if dim.find('time') != -1 or dim == 't':
                dim_to_mean = dims[dim]
                dim_to_mean_number = dim_count
            else:
                if dim.find('depth') != -1 or dim == 'z':
                    dimdepth = dim
                if dim.find('y') != -1:
                    dimy = dim
                dims_to_write = dims_to_write + [dim]

        print 'dim_to_mean : ',dim_to_mean
        print 'dim_to_mean_number : ',dim_to_mean_number
        print 'dim_to_mean length : ',len(dim_to_mean)

        # Can apply meaning to the boolean mask - will treat True as 1 and False as 0 
        # - then we convert back to boolean.
        # This allows us to have any dimension as the time dimension.
        # Hopefully the mask doesn't change with time :)
        output_mask = np.bool_(np.mean(var[:].mask,axis=dim_to_mean_number))
       
        if nstddev is not None:
            var_mean_diff = ma.average(var[:],axis=dim_to_mean_number)
            var_stddev = ma.std(var1[:],axis=dim_to_mean_number)
            var_sigdiff = ma.where(abs(var_mean_diff) < nstddev*var_stddev, np.zeros(var_mean_diff.shape), var_mean_diff)
            if "_FillValue" in var_template.ncattrs():
                outfile.createVariable(varname+'sigdiff',datatype='f',dimensions=dims_to_write,fill_value=var_template._FillValue)
            else:
                outfile.createVariable(varname+'sigdiff',datatype='f',dimensions=dims_to_write)
            var_sigdiff[output_mask]=var[:].fill_value
            outfile.variables[varname+'sigdiff'][:] = var_sigdiff[:]

        elif siglevel is not None:
            var_mean_diff = ma.average(var,axis=dim_to_mean_number)
            var_tstat = ma.zeros(var[0].shape)
            var_pvalue = ma.zeros(var[0].shape)
            (var_tstat,var_pvalue) = stats.ttest_1samp(var,0.0,axis=dim_to_mean_number)
            var_sigdiff = ma.where(var_pvalue > siglevel, np.zeros(var_mean_diff.shape), var_mean_diff)
            if "_FillValue" in var_template.ncattrs():
                outfile.createVariable(varname+'_tstat',datatype='f',dimensions=dims_to_write,fill_value=var_template._FillValue)
                outfile.createVariable(varname+'_pvalue',datatype='f',dimensions=dims_to_write,fill_value=var_template._FillValue)
                outfile.createVariable(varname+'_meandiff',datatype='f',dimensions=dims_to_write,fill_value=var_template._FillValue)
                outfile.createVariable(varname+'_sigdiff',datatype='f',dimensions=dims_to_write,fill_value=var_template._FillValue)
            else:
                outfile.createVariable(varname+'_tstat',datatype='f',dimensions=dims_to_write)
                outfile.createVariable(varname+'_pvalue',datatype='f',dimensions=dims_to_write)
                outfile.createVariable(varname+'_meandiff',datatype='f',dimensions=dims_to_write)
                outfile.createVariable(varname+'_sigdiff',datatype='f',dimensions=dims_to_write)
            var_tstat[output_mask]=var[:].fill_value
            outfile.variables[varname+'_tstat'][:] = var_tstat[:]
            var_pvalue[output_mask]=var[:].fill_value
            outfile.variables[varname+'_pvalue'][:] = var_pvalue[:]
            var_mean_diff[output_mask]=var[:].fill_value
            outfile.variables[varname+'_meandiff'][:] = var_mean_diff[:]
            var_sigdiff[output_mask]=var[:].fill_value
            outfile.variables[varname+'_sigdiff'][:] = var_sigdiff[:]

    try:
        infile1.close()
    except(NameError):
        pass
    try:
        infile2.close()
    except(NameError):
        pass
    try:
        diff_file.close()
    except(NameError):
        pass
    outfile.close()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles1", action="store",dest="filenames1_in", nargs="+", default=None,
                    help="name(s) of input file(s)")
    parser.add_argument("-j", "--infiles2", action="store",dest="filenames2_in", nargs="+", default=None,
                    help="name(s) of input file(s)")
    parser.add_argument("-d", "--diff_files", action="store",dest="differences_in", nargs="+", default=None,
                    help="name(s) of input file(s)")
    parser.add_argument("-o", "--outfile", action="store",dest="filename_out",default=None,
                    help="plot to output file - filetype by extension. If unset plot to GUI.")
    parser.add_argument("-v", "--var_list", action="store",dest="var_list",nargs='+',default=None,
                    help="variables to time-mean")
    parser.add_argument("-n", "--nstddev", action="store",dest="nstddev",default=None,type=float,
                    help="number of std dev cut off")
    parser.add_argument("-s", "--siglevel", action="store",dest="siglevel",default=None,type=float,
                    help="significance level cut off")

    args = parser.parse_args()

    sigdiff( filenames1_in=args.filenames1_in, filenames2_in=args.filenames2_in, differences_in=args.differences_in, 
             filename_out=args.filename_out, var_list=args.var_list, nstddev=args.nstddev, siglevel=args.siglevel ) 
