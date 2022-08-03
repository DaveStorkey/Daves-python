#! /usr/bin/env python

'''
Routine to calculate the time-variance of a field
excluding the seasonal variability. 

Method:

Assume one time level per file. Calculate the variance over the files
for a particular time of year. Then combine all the variances by averaging. 
This is equivalent to removing the seasonal signal from the raw data 
first and then calculating the variance. 

@author: Dave Storkey
@date: March 2017
'''

import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import re
import os.path

def nonseasonal_variance(filenames_in=None,filename_out=None,varname_list=None):

    #### Open output file and copy relevant dimensions  ####
    #### and coordinate variables from input file       ####

    outfile = nc.Dataset(filename_out,'w')
    template_file = nc.Dataset(filenames_in[0])
    dims = template_file.dimensions
    for dim in dims.keys():
        if dim.find('time') == -1 and dim != 't':
            outfile.createDimension(dim, len(dims[dim]))
        else:
            dim_time = dim

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

    ## Calculate the full variance and the nonseasonal variance.

    # assume input files look like "jobname_1m_yyyymmdd_...." etc.
    # converting to a set and back removes duplication
    years = list(set([filename.split("_")[2][0:4] for filename in filenames_in]))
    years.sort()
    nyears = len(years)
    ntimes = len(filenames_in)
    if ntimes%nyears != 0:
        nyears=nyears-1
    ntimes_per_year = ntimes/nyears
    print 'number of years : ',nyears
    print 'number of files (total) : ',ntimes
    print 'number of files per year : ',ntimes_per_year
    filenames_firstyear = filenames_in[0:ntimes_per_year]

    for varname in varname_list:

        # Calculate the full time-series mean
        files_in = nc.MFDataset(filenames_in)
        for itime in range(ntimes):
            var_in = files_in.variables[varname][itime]
            # accumulate
            if itime == 0:
                var_mean = var_in[:]
            else:
                var_mean[:] = var_mean[:] + var_in[:]
        var_mean[:] = var_mean[:]/ntimes
        
        # Calculate the full and nonseasonal variance
        nfiles_last = -1
        for filename in filenames_firstyear:
            year0 = filename.split("_")[2][0:4]
            years_to_loop = list(years)
            years_to_loop.remove(year0)
            filename_list=[filename]
            for year in years_to_loop:
                 if os.path.isfile(filename.replace(year0,year)):
                     filename_list.append(filename.replace(year0,year))
            print 'len(filename_list) : ',len(filename_list)
            print 'filename_list : ',filename_list
            if nfiles_last != -1 and len(filename_list) != nfiles_last:
                raise Exception("Error: inconsistent numbers of files to calculate stats")            
            else:
                nfiles_last = len(filename_list)
            files_in = nc.MFDataset(filename_list)

            # Calculate time mean of the field incrementally (to save memory).
            for iyear in range(nyears):
                var_in = files_in.variables[varname][iyear]
                # accumulate
                if iyear == 0:
                    var_submean = var_in[:]
                else:
                    var_submean[:] = var_submean[:] + var_in[:]
            var_submean[:] = var_submean[:]/nyears
            
            # Calculate variances of the field incrementally (to save memory).
            for iyear in range(nyears):
                var_in = files_in.variables[varname][iyear]
                # accumulate
                if filename == filenames_firstyear[0] and iyear == 0:
                    var_variance = (var_in[:]-var_mean[:]) * (var_in[:]-var_mean[:])
                else:
                    var_variance[:] = var_variance[:] + (var_in[:]-var_mean[:]) * (var_in[:]-var_mean[:])
                if iyear == 0:
                    var_subvariance = (var_in[:]-var_submean[:]) * (var_in[:]-var_submean[:])
                else:
                    var_subvariance[:] = var_subvariance[:] + (var_in[:]-var_submean[:]) * (var_in[:]-var_submean[:])

            if filename == filenames_firstyear[0]:
                var_nsvariance = var_subvariance[:]
            else:
                var_nsvariance[:] = var_nsvariance[:] + var_subvariance[:]
             
            files_in.close()

        var_variance[:] = var_variance[:] / len(filenames_in)
        var_nsvariance[:] = var_nsvariance[:] / len(filenames_in)

        var_template = template_file.variables[varname]
        dims_to_write = list(var_template.dimensions)
        dims_to_write.remove(dim_time)
        if "_FillValue" in var_template.ncattrs():
            outfile.createVariable(varname+'_stddev',datatype='f',dimensions=dims_to_write,fill_value=var_template._FillValue)
            outfile.createVariable(varname+'_ns_stddev',datatype='f',dimensions=dims_to_write,fill_value=var_template._FillValue)
        else:
            outfile.createVariable(varname+'_stddev',datatype='f',dimensions=dims_to_write)
            outfile.createVariable(varname+'_ns_stddev',datatype='f',dimensions=dims_to_write)

        var_stddev = ma.sqrt(var_variance)
        var_stddev[var_stddev.mask] = var_stddev.fill_value
        outfile.variables[varname+'_stddev'][:] = var_stddev
        var_nsstddev = ma.sqrt(var_nsvariance)
        var_nsstddev[var_nsstddev.mask] = var_nsstddev.fill_value
        outfile.variables[varname+'_ns_stddev'][:] = var_nsstddev

    outfile.close()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", action="store",dest="filenames_in", nargs="+", default=None,
                    help="name(s) of input file(s)")
    parser.add_argument("-o", "--outfile", action="store",dest="filename_out",default=None,
                    help="plot to output file - filetype by extension. If unset plot to GUI.")
    parser.add_argument("-v", "--varname_list", action="store",dest="varname_list",nargs='+',default=None,
                    help="variables to time-mean")

    args = parser.parse_args()

    nonseasonal_variance( filenames_in=args.filenames_in, filename_out=args.filename_out, varname_list=args.varname_list ) 
