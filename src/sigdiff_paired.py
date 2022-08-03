#! /usr/bin/env python

'''
Routine to calculate the "signficant" part of a time-mean difference
field given the time series of the differences or the two time series
of the fields to be differenced. Insignificant differences are set
to zero in the output field. 

This version calculates a *paired* Student's t-test in a memory-efficient
way.

@author: Dave Storkey
@date: Feb 2017
'''

import netCDF4 as nc
import numpy as np
import numpy.ma as ma
from scipy import stats

def sigdiff(filenames1_in=None,filenames2_in=None,filename_out=None,
            varname_list=None,siglevel=None):

    #### Default sig level ####

    if siglevel is None:
        print "Defaulting to significance level of 5%"
        siglevel = 0.05

    #### Open input files #### 

    need_to_make_diffs = True

    if filenames1_in is not None:

        if len(filenames1_in) == 1:
            # only one file:
            infile1 = nc.Dataset(filenames1_in[0])
        else:
            infile1 = nc.MFDataset(filenames1_in)

        if varname_list is None:
            varname_list = infile1.variables.keys()

        template_file = infile1

    if filenames2_in is not None:

        if len(filenames2_in) == 1:
            # only one file:
            infile2 = nc.Dataset(filenames2_in[0])
        else:
            infile2 = nc.MFDataset(filenames2_in)

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
            if latname in varname_list:
                varname_list.remove(latname)
            break

    for lonname in 'lon','longitude','nav_lon':
        try:
            lon = template_file.variables[lonname]
        except KeyError:
            pass
        else:
            outfile.createVariable(lonname,datatype='f',dimensions=lon.dimensions)
            outfile.variables[lonname][:] = lon[:]
            if lonname in varname_list:
                varname_list.remove(lonname)
            break

    for depthname in 'depth','deptht','depthu','depthv','depthw':
        try:
            depth = template_file.variables[depthname]
        except KeyError:
            pass
        else:
            outfile.createVariable(depthname,datatype='f',dimensions=depth.dimensions)
            outfile.variables[depthname][:] = depth[:]
            if depthname in varname_list:
                varname_list.remove(depthname)
            break

    #### Loop over variable names in the list and calculate significant differences ####

    for varname in varname_list:

        var_template = template_file.variables[varname]
        dim_count=-1
        dim_to_mean=-1
        dims_to_write=[]
        dims_var = var_template.dimensions
        for dim in dims_var:
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

        ntimes = len(dim_to_mean)

        print 'dim_to_mean : ',dim_to_mean
        print 'dim_to_mean_number : ',dim_to_mean_number
        print 'dim_to_mean length : ',ntimes
        
        # 1. Calculate time mean of the difference field incrementally (to save memory).
        var_diff_mean = ma.zeros(var_template[0].shape)
        for t in range(ntimes):

            var1 = infile1.variables[varname][t]
            var2 = infile2.variables[varname][t]
            var_diff = var2 - var1
            # accumulate
            var_diff_mean[:] = var_diff_mean[:] + var_diff[:]

        var_diff_mean[:] = var_diff_mean[:]/ntimes

        # 2. Calculate variance incrementally (to save memory).
        var_diff_variance = ma.zeros(var_template[0].shape)
        for t in range(ntimes):

            # wasteful of cpu to read in var1 and var2 again and calculate difference for this 
            # time level but we are aiming for memory efficiency here. 
            var1 = infile1.variables[varname][t]
            var2 = infile2.variables[varname][t]
            var_diff = var2 - var1
            # accumulate
            var_diff_variance[:] = var_diff_variance[:] + (var_diff[:]-var_diff_mean[:])*(var_diff[:]-var_diff_mean[:])

        var_diff_variance[:] = var_diff_variance[:]/ntimes

        var_diff_stddev = ma.sqrt(var_diff_variance)
        var_tstat = ma.abs(var_diff_mean)*ma.sqrt(ntimes)/var_diff_stddev
        # stats.t.sf doesn't handle masking (and there is no equivalent mstats function)
        # so we create an unmasked array and then mask it.
        var_pval_unmasked = stats.t.sf(ma.abs(var_tstat), ntimes-1)*2  # two-sided pvalue 
        var_pval = ma.copy(var_tstat)
        var_pval.harden_mask()  # to prevent masked points getting unmasked in the next line
        var_pval[:] = var_pval_unmasked[:]
        var_sigdiff = ma.where(var_pval > siglevel, np.zeros(var_diff_mean.shape), var_diff_mean)

        if "_FillValue" in var_template.ncattrs():
            outfile.createVariable(varname+'_diff_stddev',datatype='f',dimensions=dims_to_write,fill_value=var_template._FillValue)
            outfile.createVariable(varname+'_tstat',datatype='f',dimensions=dims_to_write,fill_value=var_template._FillValue)
            outfile.createVariable(varname+'_pval',datatype='f',dimensions=dims_to_write,fill_value=var_template._FillValue)
            outfile.createVariable(varname+'_diffmean',datatype='f',dimensions=dims_to_write,fill_value=var_template._FillValue)
            outfile.createVariable(varname+'_sigdiff',datatype='f',dimensions=dims_to_write,fill_value=var_template._FillValue)
        else:
            outfile.createVariable(varname+'_diff_stddev',datatype='f',dimensions=dims_to_write)
            outfile.createVariable(varname+'_tstat',datatype='f',dimensions=dims_to_write)
            outfile.createVariable(varname+'_pval',datatype='f',dimensions=dims_to_write)
            outfile.createVariable(varname+'_diffmean',datatype='f',dimensions=dims_to_write)
            outfile.createVariable(varname+'_sigdiff',datatype='f',dimensions=dims_to_write)

        var_diff_stddev[var_diff_stddev.mask]=var_diff_stddev.fill_value
        outfile.variables[varname+'_diff_stddev'][:] = var_diff_stddev[:]
        var_tstat[var_tstat.mask]=var_tstat.fill_value
        outfile.variables[varname+'_tstat'][:] = var_tstat[:]
        var_pval[var_pval.mask]=var_pval.fill_value
        outfile.variables[varname+'_pval'][:] = var_pval[:]
        var_diff_mean[var_diff_mean.mask]=var_diff_mean.fill_value
        outfile.variables[varname+'_diffmean'][:] = var_diff_mean[:]
        var_sigdiff[var_sigdiff.mask]=var_sigdiff.fill_value
        outfile.variables[varname+'_sigdiff'][:] = var_sigdiff[:]

    try:
        infile1.close()
    except(NameError):
        pass
    try:
        infile2.close()
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
    parser.add_argument("-o", "--outfile", action="store",dest="filename_out",default=None,
                    help="plot to output file - filetype by extension. If unset plot to GUI.")
    parser.add_argument("-v", "--varname_list", action="store",dest="varname_list",nargs='+',default=None,
                    help="variables to time-mean")
    parser.add_argument("-s", "--siglevel", action="store",dest="siglevel",default=None,type=float,
                    help="significance level cut off")

    args = parser.parse_args()

    sigdiff( filenames1_in=args.filenames1_in, filenames2_in=args.filenames2_in, 
             filename_out=args.filename_out, varname_list=args.varname_list, siglevel=args.siglevel ) 
