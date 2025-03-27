#! /usr/bin/env python

'''
Wrapper to make it easy to call cube.collapsed from the command line. 

@author: Dave Storkey
@date: March 2025
'''

import iris
import iris.analysis

def read_cube(filename,fieldname):
    '''
    Read a variable from a netcdf file as an Iris cube.
    Try matching varname to the standard_name first and if that fails
    try matching to the variable name in the file. 
    '''

    try:
        cube = iris.load_cube(filename, fieldname)
    except iris.exceptions.ConstraintMismatchError:
        try: 
            cube = iris.load_cube(filename,
                      iris.Constraint(cube_func=lambda nemoname: nemoname.var_name == fieldname))
        except iris.exceptions.ConstraintMismatchError:
            raise Exception('Could not find '+fieldname+' in file '+filename)

    return cube

def get_weights(wgtsfile,wgtsname,cube):
    '''
    Try to read in a weights field from the specified file, or if the 
    file is specified as "measures", try to read the weights field as
    a CellMeasure of the supplied cube.
    '''

    if wgtsfile == "measures":
        for cell_measure in cube.cell_measures():
            if cell_measure.standard_name == wgtsname:
                wgts = cell_measure.core_data()
                break
            else:
                pass
        else:
            raise Exception("Could not find "+wgtsname+" in cell measures of "+cube.var_name)
    else:
        wgts = read_cube(wgtsfile,wgtsname)
    
    return wgts

def reduce_fields(infile=None,invars=None,coords=None,wgtsfiles=None,wgtsnames=None,
                  aggr=None,outfile=None,east=None,west=None,south=None,north=None,
                  top=None,bottom=None):

    aggregators = { "mean"     :  iris.analysis.MEAN ,
                    "min"      :  iris.analysis.MIN  ,
                    "max"      :  iris.analysis.MAX    }

    if infile is None:
        raise Exception("Error: must specify input file")

    if aggr is None:
        aggr="mean"

    if outfile is None:
        outfile=infile.split(".")[:-1]+"_reduced."+infile.split(".")[-1]

    if invars is None:
        cubes = iris.load(infile)
    else:
        cubes = [read_cube(infile,varname) for varname in invars]
        
    if coords is None:
        coords = "time"

    if wgtsnames is not None:
        if not isinstance(wgtsnames,list):
            wgtsnames=[wgtsnames]
        if wgtsfiles is None:
            print("No wgtsfile specified. Looking for weights in input file.")
            wgtsfiles = [infile]
        elif not isinstance(wgtsfiles,list):
            wgtsfiles=[wgtsfiles]
        if len(wgtsfiles) == 1:
            wgtsfiles = wgtsfiles*len(wgtsnames)
        wgtsfiles = [infile if wf == "self" else wf for wf in wgtsfiles]
        if len(wgtsfiles) != len(wgtsnames):
            raise Exception("Must specify one weights file or the same number as the number of weights fields")
        if len(wgtsfiles) > 1 and wgtsfiles[0] == "measures":
            # multiply won't work with a CellMeasure object as the first argument
            wgtsfiles.append(wgtsfiles.pop(0))
            wgtsnames.append(wgtsnames.pop(0))
        wgts_list = [get_weights(wgtsfile,wgtsname,cubes[0]) for (wgtsfile,wgtsname) in zip(wgtsfiles,wgtsnames)]
        wgts=wgts_list[0]        
        if len(wgtsnames) > 1:
            for wgts_to_multiply in wgts_list[1:]:
                wgts = iris.analysis.maths.multiply(wgts, wgts_to_multiply, in_place=True)
    else:
        wgts = None
                
    cubes_reduced=[cube.collapsed(coords, aggregators[aggr], weights=wgts) for cube in cubes]

    iris.save(cubes_reduced, outfile)


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", action="store", dest="infile", 
                         help="names of input file")
    parser.add_argument("-v", "--vars", action="store", dest="invars", nargs="+", 
                         help="names of input variables")
    parser.add_argument("-G", "--wgtsfiles", action="store", dest="wgtsfiles", nargs="+",
                         help="names of weights file or 'self' if input file or 'measures' if a cell measure")
    parser.add_argument("-g", "--wgtsnames", action="store", dest="wgtsnames", nargs="+",
                         help="names of weighting variable")
    parser.add_argument("-o", "--outfile", action="store", dest="outfile",
                         help="name of output file (format by extension)")
    parser.add_argument("-c", "--coords", action="store",dest="coords",nargs="+",
                         help="name of coordinates to reduce over (default time)")
    parser.add_argument("-A", "--aggr", action="store",dest="aggr",
                         help="name of aggregator: mean, max, min")
    parser.add_argument("-W", "--west", action="store",dest="west",
                         help="western limit of area to reduce")
    parser.add_argument("-E", "--east", action="store",dest="east",
                         help="eastern limit of area to reduce")
    parser.add_argument("-S", "--south", action="store",dest="south",
                         help="southern limit of area to reduce")
    parser.add_argument("-N", "--north", action="store",dest="north",
                         help="northern limit of area to reduce")
    parser.add_argument("-T", "--top", action="store",dest="top",
                         help="top limit of volume to reduce")
    parser.add_argument("-B", "--bottom", action="store",dest="bottom",
                         help="bottom limit of volume to reduce")
    args = parser.parse_args()

    reduce_fields(infile=args.infile,invars=args.invars,outfile=args.outfile,
                  wgtsfiles=args.wgtsfiles,wgtsnames=args.wgtsnames,coords=args.coords,aggr=args.aggr,
                  south=args.south,north=args.north,west=args.west,east=args.east,
                  top=args.top,bottom=args.bottom)



