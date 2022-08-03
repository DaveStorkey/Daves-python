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

def sigdiff(filenames1_in=None,filenames2_in=None,filename_out=None,
            varname_list=None,siglevel=None,nstddev=None):

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
                dim_to_mean = dim
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
        
        # 1. Calculate time mean of each population incrementally (to save memory).
        var_mean_list = [ma.zeros(var_template[0].shape),ma.zeros(var_template[0].shape)]
        for t in range(ntimes):

            var_list=[]
            var_list.append(infile1.variables[varname][t])
            var_list.append(infile2.variables[varname][t])
            
            for (var,var_mean) in zip(var_list,var_mean_list):            
                # accumulate
                var_mean[:] = var_mean[:] + var[:]

        for var_mean in var_mean_list:
            var_mean[:] = var_mean[:]/ntimes

        # 2. Calculate variance incrementally (to save memory).
        var_variance_list=[ma.zeros(var_template[0].shape),ma.zeros(var_template[0].shape)]
        for t in range(ntimes):

            var_list=[]
            var_list.append(infile1.variables[varname][t])
            var_list.append(infile2.variables[varname][t])

            for (var,var_mean,var_variance) in zip(var_list,var_mean_list,var_variance_list):            
                # accumulate
                var_variance[:] = var_variance[:] + (var[:]-var_mean[:])*(var[:]-var_mean[:])

        for var_variance in var_variance_list:
            var_variance[:] = var_variance[:]/ntimes

        var_tstat = ma.zeros(var_template[0].shape)
        var_pvalue = ma.zeros(var_template[0].shape)
        # Only available at scipy v 0.16.0 D'oh!!
        #(var_tstat,var_pvalue) = stats.ttest_ind_from_stats(var_mean_list[0], var_stddev_list[0], ntimes, 
        #                                                    var_mean_list[1], var_stddev_list[1], ntimes )

        var_mean_diff = var_mean_list[1]-var_mean_list[0]
        var_pooled_stddev = ma.sqrt( 0.5*(var_variance_list[0]+var_variance_list[1]) )
        var_tstat = ma.abs(var_mean_diff)/var_pooled_stddev
        # 1.96 is the significance threshold for a p-value of 0.05 and an infinite number of degrees of freedom. 
        var_sigdiff = ma.where(var_tstat > 1.96, var_mean_diff, np.zeros(var_mean.shape))

        if "_FillValue" in var_template.ncattrs():
            outfile.createVariable(varname+'_mean1',datatype='f',dimensions=dims_to_write,fill_value=var_template._FillValue)
            outfile.createVariable(varname+'_stddev1',datatype='f',dimensions=dims_to_write,fill_value=var_template._FillValue)
            outfile.createVariable(varname+'_mean2',datatype='f',dimensions=dims_to_write,fill_value=var_template._FillValue)
            outfile.createVariable(varname+'_stddev2',datatype='f',dimensions=dims_to_write,fill_value=var_template._FillValue)
            outfile.createVariable(varname+'_tstat',datatype='f',dimensions=dims_to_write,fill_value=var_template._FillValue)
            outfile.createVariable(varname+'_meandiff',datatype='f',dimensions=dims_to_write,fill_value=var_template._FillValue)
            outfile.createVariable(varname+'_sigdiff',datatype='f',dimensions=dims_to_write,fill_value=var_template._FillValue)
        else:
            outfile.createVariable(varname+'_mean1',datatype='f',dimensions=dims_to_write)
            outfile.createVariable(varname+'_stddev1',datatype='f',dimensions=dims_to_write)
            outfile.createVariable(varname+'_mean2',datatype='f',dimensions=dims_to_write)
            outfile.createVariable(varname+'_stddev2',datatype='f',dimensions=dims_to_write)
            outfile.createVariable(varname+'_tstat',datatype='f',dimensions=dims_to_write)
            outfile.createVariable(varname+'_meandiff',datatype='f',dimensions=dims_to_write)
            outfile.createVariable(varname+'_sigdiff',datatype='f',dimensions=dims_to_write)

        var_mean_list[0][var.mask]=var.fill_value
        outfile.variables[varname+'_mean1'][:] = var_mean_list[0][:]
        var_mean_list[1][var.mask]=var.fill_value
        outfile.variables[varname+'_mean2'][:] = var_mean_list[1][:]
        var_variance_list[0][var.mask]=var.fill_value
        outfile.variables[varname+'_stddev1'][:] = ma.sqrt(var_variance_list[0][:])
        var_variance_list[1][var.mask]=var.fill_value
        outfile.variables[varname+'_stddev2'][:] = ma.sqrt(var_variance_list[1][:])
        var_tstat[var.mask]=var.fill_value
        outfile.variables[varname+'_tstat'][:] = var_tstat[:]
        var_mean_diff[var.mask]=var.fill_value
        outfile.variables[varname+'_meandiff'][:] = var_mean_diff[:]
        var_sigdiff[var.mask]=var.fill_value
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
    parser.add_argument("-n", "--nstddev", action="store",dest="nstddev",default=None,type=float,
                    help="number of std dev cut off")
    parser.add_argument("-s", "--siglevel", action="store",dest="siglevel",default=None,type=float,
                    help="significance level cut off")

    args = parser.parse_args()

    sigdiff( filenames1_in=args.filenames1_in, filenames2_in=args.filenames2_in, 
             filename_out=args.filename_out, varname_list=args.varname_list, nstddev=args.nstddev, siglevel=args.siglevel ) 
