#! /usr/bin/env python

'''
Use Iris to convert from PP to netcdf. 

@author: Dave Storkey
@date: Dec 2023
'''

import iris
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="name of input file")
    args = parser.parse_args()

    file_in = args.infile
    cubes = iris.load(file_in)

    file_out = file_in.replace('.pp','.nc')
    iris.save(cubes, file_out)


