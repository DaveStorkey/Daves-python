#! /usr/bin/env python

'''
Take a file with no meta data (eg. a NEMO restart file) and use the meta data
from another file to wrap specified fields from the original file. 

@author: Dave Storkey
@date: Oct 2025
'''

import iris
import iris.analysis
import iris.util

def read_cube(filename,fieldname):
    '''
    Read a variable from a netcdf file as an Iris cube.
    Try to match name to standard_name, long_name or var_name
    Remove duplicate time dimension if necessary.
    '''

    constraints = [ iris.NameConstraint(standard_name=fieldname),
                    iris.NameConstraint(long_name=fieldname),
                    iris.NameConstraint(var_name=fieldname) ]

    for constraint in constraints:
        try:
            cube = iris.load_cube(filename, constraint)
        except iris.exceptions.ConstraintMismatchError:
            pass
        else:
            break
    else:
        raise Exception("Could not find field ",fieldname," in file ",filename)

    for coord in cube.coords(dim_coords=True):
        if coord.var_name == "time_counter" and "time_centered" in [coord.var_name for coord in cube.coords()]:
            cube.remove_coord(coord)
            iris.util.promote_aux_coord_to_dim_coord(cube,"time")
            break

    return cube

def metawrap(infile=None,invars=None,metafile=None,metavars=None,outfile=None,maskfields=None,thick=None):

    if infile is None:
        raise Exception("Error: must specify input file")

    if invars is None:
        raise Exception("Error: must specify at least one input variable")

    if metafile is None:
        raise Exception("Error: must specify file containing vars with meta data")

    if metavars is None:
        raise Exception("Error: must specify variables with meta data")
    elif len(metavars) != len(invars):
        raise Exception("Error: number of meta variables must match number of input variables")

    if outfile is None:
        outfile=".".join(infile.split(".")[:-1])+"_wrapped."+infile.split(".")[-1]

    cubes_in=[]
    for var in invars:
        cubes_in.append(read_cube(infile,var))
    # save the input fill_value attribute so that we can force the output fields
    # to use the same value.
    fill_value=cubes_in[0].data.fill_value    

    cubes_meta=[]
    for var in metavars:
        cubes_meta.append(read_cube(metafile,var))

    for cube_in, cube_meta in zip(cubes_in, cubes_meta):
        if maskfields:
            mask = cube_meta.data.mask.copy()
        cube_meta.data[:] = cube_in.data[:]
        if maskfields:
            cube_meta.data.mask[:] = mask[:]
        
    if thick:
        cubes_meta.append(read_cube(metafile,"cell_thickness"))
        
    iris.save(cubes_meta, outfile, fill_value=fill_value)

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", action="store", dest="infile", 
                         help="name of input file")
    parser.add_argument("-v", "--invars", action="store", dest="invars", nargs="+", 
                         help="names of input variables")
    parser.add_argument("-M", "--metafile", action="store", dest="metafile",
                         help="name of file with meta data")
    parser.add_argument("-m", "--wgtsvars", action="store", dest="metavars", nargs="+",
                         help="names of variables with meta data")
    parser.add_argument("-o", "--outfile", action="store", dest="outfile",
                         help="name of output file (format by extension)")
    parser.add_argument("-K", "--mask", action="store_true", dest="maskfields",
                         help="apply mask from meta field to input field")
    parser.add_argument("-T", "--thick", action="store_true", dest="thick",
                         help="write thickness field from meta file to output file")
    args = parser.parse_args()

    metawrap(infile=args.infile,invars=args.invars,outfile=args.outfile,
             metafile=args.metafile,metavars=args.metavars,maskfields=args.maskfields,
             thick=args.thick)



