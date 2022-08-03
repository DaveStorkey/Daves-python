#! /usr/bin/env python

'''
Integrate field (x optional scalar) over specified model levels.

@author: Dave Storkey
@date: March 2022
'''

import xarray as xr
import numpy as np

def vert_integrate(file_in=None, var_in=None, file_out=None, var_out=None, levels=None, factor=None, avg=None):

    if not file_in or not var_in or not file_out:
        raise Exception("Error: must specify file_in, var_in and file_out.")

    if not isinstance(var_in,list):
        var_in=[var_in]

    if not var_out:
        var_out=[]
        if avg:
            suffix="_depth_avg"
        else:
            suffix="_depth_int"
        for var in var_in:
            var_out.append(var+suffix)
    else:
        if len(var_out) != len(var_in):
            raise Exception("Error: number of output variables must match number of input variables")

    if not factor:
        factor=1.0

    for invar,outvar in zip(var_in,var_out):
        with xr.open_dataset(file_in) as indata:
            try:
                field = getattr(indata,invar)
            except AttributeError:
                raise Exception("Can't find "+invar+" in input file.")
            for thkvar in ["thkcello","e3t","e3tn"]:
                try:
                    thick = getattr(indata,thkvar)
                except AttributeError:
                    pass
                else:
                    break
            else:
                raise Exception("Can't find thickness field in input file.")

        # vertically integrate
        field_int3D = field * thick * factor
        if avg:
            if levels:
                thick_int = thick.isel(deptht=slice(levels[0],levels[1])).sum(dim="deptht")
            else:
                thick_int = thick.sum(dim="deptht")
        else:
            thick_int = 1.0
        if levels:
            field_int = field_int3D.isel(deptht=slice(levels[0],levels[1])).sum(dim="deptht") / thick_int
        else:
            field_int = field_int3D.sum(dim="deptht") / thick_int

        if invar == var_in[0]:
            outdata = field_int.to_dataset(name=outvar)
        else:
            outdata[outvar] = field_int

    outdata.to_netcdf(file_out)
                
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--file_in", action="store",dest="file_in",
                         help="input file")
    parser.add_argument("-v", "--var_in", action="store",dest="var_in",nargs="+",
                         help="name of input variable")
    parser.add_argument("-o", "--file_out", action="store",dest="file_out",
                         help="output file")
    parser.add_argument("-w", "--var_out", action="store",dest="var_out",nargs="+",
                         help="name of output variable")
    parser.add_argument("-L", "--levels", action="store",dest="levels",type=int,nargs=2,
                         help="start and end levels for integration")
    parser.add_argument("-X", "--factor", action="store",dest="factor",type=float,
                         help="scalar multiplicative factor")
    parser.add_argument("-A", "--avg", action="store_true",dest="avg",
                         help="calculate depth mean rather than depth integral")

    args = parser.parse_args()

    vert_integrate(file_in=args.file_in,var_in=args.var_in,file_out=args.file_out,var_out=args.var_out,
                   levels=args.levels,factor=args.factor,avg=args.avg)

