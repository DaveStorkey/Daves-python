#! /usr/bin/env python

'''
Use Iris to mean a series of PP/netcdf files over a single dimension
(default "time") and write to a PP/netcdf file. 

@author: Dave Storkey
@date: Dec 2023
'''

import iris
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", action="store", dest="infiles", nargs="+", 
                         help="names of input files")
    parser.add_argument("-o", "--outfile", action="store", dest="outfile",
                         help="name of output file (format by extension)")
    parser.add_argument("-d", "--dim", action="store",dest="dim",
                         help="name of dimension to mean (default time)")
    args = parser.parse_args()

    files_in = args.infiles
    print("files_in : ",files_in)
    cubes = iris.load(files_in)

    dim_to_mean = args.dim
    if dim_to_mean is None:
        dim_to_mean = "time"
    file_out = args.outfile
    if file_out is None:
        file_out = files_in[0].replace(".pp","_"+dim_to_mean+"_mean.pp")
        file_out = files_in[0].replace(".nc","_"+dim_to_mean+"_mean.nc")

    cubes_mean=[cube.collapsed(dim_to_mean,iris.analysis.MEAN) for cube in cubes]

    iris.save(cubes_mean, file_out)


