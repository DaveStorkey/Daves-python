#! /usr/bin/env python

'''
Volume integrate a 3D T-grid field or area integrate a 2D T-grid field.
Outputs the weighted integral and the average (mean). 

Possibility to specify bounds on lon/lat and depth or supply a 2D or 3D
mask to specify the volume/area of integration. 

NB. Currently assumes only 1 time level of data in the input file.

@author: Dave Storkey
@date: May 2022
'''

import xarray as xr
import numpy as np

def geo_integrate(grid=None, file_in=None, var_in=None, file_out=None, 
                  maskfile=None, maskvar=None, levels=None, factor=None,
                  gt=None, lt=None, north=None, south=None, east=None, west=None):

    if grid is None or file_in is None or var_in is None or file_out is None:
        raise Exception("Error: must specify grid, file_in, var_in and file_out.")

    if not isinstance(var_in,list):
        var_in=[var_in]

    print("var_in : ",var_in[:])

    var_int=[]
    var_avg=[]
    for var in var_in:
        var_int.append(var+"_int")
        var_avg.append(var+"_avg")

    if factor is None:
        factor=1.0

    with xr.open_dataset(grid) as griddata:
        try:
            tmask = griddata.tmask.squeeze()
        except AttributeError:
            raise Exception("Can't find tmask in grid file.")
        sfactor=[]
        for sfname in ['e1t','e2t']: 
            try:
                sfactor.append(getattr(griddata,sfname).squeeze())
            except AttributeError:
                raise Exception("Can't find "+sfname+" in grid file.")
        if north is not None or south is not None or east is not None or west is not None:
            try:
                lons = griddata.glamt.squeeze()
            except AttributeError:
                raise Exception("Can't find longitude field glamt in grid file.")
            try:
                lats = griddata.gphit.squeeze()
            except AttributeError:
                raise Exception("Can't find latitude field gphit in grid file.")
   
    fields=[]
    fieldshapes=[]
    space_metrics=[]
    with xr.open_dataset(file_in) as indata:
        thick = None
        for invar in var_in:
            try:
                fields.append(getattr(indata,invar).squeeze())
            except AttributeError:
                raise Exception("Can't find "+invar+" in input file.")
            if len(fields[-1].shape) == 2:
                fieldshapes.append("2D")
                space_metrics.append("area")
            elif len(fields[-1].shape) == 3:
                fieldshapes.append("3D")
                space_metrics.append("volume")
                # only need thickness field for the 3D case:
                if thick is None:
                    for thkvar in ["thkcello","e3t","e3tn"]:
                        try:
                            thick = getattr(indata,thkvar).squeeze()
                        except AttributeError:
                            pass
                        else:
                            break
                    else:
                        raise Exception("Can't find thickness field in input file.")
            else:
                raise Exception("Error: can only handle 2D or 3D fields (possibly with degenerate time dim).")
       
    if maskfile is not None:
        if maskvar is None:
            maskvar = "mask"
        with xr.open_dataset(maskfile) as maskdata:
            try:
                mask = getattr(maskdata,maskvar)
            except(AttributeError):
                raise Exception("Can't find "+maskvar+" in "+maskfile)

    outdata = None
    for field,fieldshape,space_metric,intvar,avgvar in zip(fields,fieldshapes,space_metrics,var_int,var_avg):
        print("calculations for "+field.name+" "+intvar+" "+avgvar)
        # create mask to mask out parts of the field we want to exclude from integration.
        fldmask=None
        print('gt,lt : ',gt,lt)
        if gt is not None or lt is not None:
            fldmask = field.copy()
            fldmask.values[:] = 1.0
            if gt is not None:
                fldmask = xr.where(field > gt, fldmask, 0.0)
            if lt is not None:
                fldmask = xr.where(field < lt, fldmask, 0.0)
        if north is not None or south is not None or east is not None or west is not None:
            if fldmask is None:
                fldmask = field.copy()
                fldmask.values[:] = 1.0
            if north is not None:
                fldmask = xr.where(lats < north, fldmask, 0.0)
            if south is not None:
                fldmask = xr.where(lats > south, fldmask, 0.0)
            if east is not None:
                fldmask = xr.where(lons > east, fldmask, 0.0)
            if west is not None:
                fldmask = xr.where(lats < west, fldmask, 0.0)

        # Need tmask to have same coordinate names as the other arrays for the multiplication to work
        # and rename doesn't seem to work. 
        tmask_copy = field.copy()
        if fieldshape == "2D":
            tmask_copy.values = tmask[0].values.squeeze()
        else:
            tmask_copy.values = tmask.values
        if maskfile is not None:
            tmask_copy.values = tmask_copy.values * mask.values
        if fldmask is not None:
            tmask_copy.values = tmask_copy.values * fldmask.values

        maskdata = tmask_copy.to_dataset(name='tmask')
        maskdata.to_netcdf('geo_mask.nc')

        if fieldshape == "2D":
            weights = (sfactor[0].transpose() * sfactor[1].transpose() \
                        * tmask_copy.transpose()).transpose()
        else:
            weights = (sfactor[0].transpose() * sfactor[1].transpose() \
                        * thick.transpose() * tmask_copy.transpose()).transpose()

        # vertically sum for 3D fields
        if fieldshape == "3D":
            if levels:
                weights_2D = weights.isel(deptht=slice(levels[0],levels[1])).sum(dim="deptht")
            else:
                weights_2D = weights.sum(dim="deptht")
        else:
            weights_2D = weights

        # horizontally sum
        weights_int = weights_2D.sum(dim=["x","y"])
        if outdata is None:
            outdata = weights_int.to_dataset(name=space_metric)
        elif space_metric not in outdata:
            outdata[space_metric] = weights_int

        # multiply by weights and optional factor
        field_weighted = field * weights * factor
        # vertically sum for 3D fields
        if fieldshape == "3D":
            if levels:
                field_2D = field_weighted.isel(deptht=slice(levels[0],levels[1])).sum(dim="deptht")
            else:
                field_2D = field_weighted.sum(dim="deptht")
        else:
            field_2D = field_weighted
        # horizontally sum
        field_int = field_2D.sum(dim=["x","y"])
        field_avg = field_int / weights_int
        print("writing "+intvar+" and "+avgvar)
        outdata[intvar] = field_int
        outdata[avgvar] = field_avg

    outdata.to_netcdf(file_out)
                
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--grid", action="store",dest="grid",
                         help="grid (mesh) file (for horizontal scale factors and tmask)")
    parser.add_argument("-M", "--maskfile", action="store",dest="maskfile",
                         help="name of mask file")
    parser.add_argument("-m", "--maskvar", action="store",dest="maskvar",
                         help="mask variable name")
    parser.add_argument("-i", "--file_in", action="store",dest="file_in",
                         help="input file")
    parser.add_argument("-v", "--var_in", action="store",dest="var_in",nargs="+",
                         help="name of input variable")
    parser.add_argument("-o", "--file_out", action="store",dest="file_out",
                         help="output file")
    parser.add_argument("-N", "--north", action="store",dest="north",type=float,
                         help="northern latitude limit of integration")
    parser.add_argument("-S", "--south", action="store",dest="south",type=float,
                         help="southern latitude limit of integration")
    parser.add_argument("-E", "--east", action="store",dest="east",type=float,
                         help="eastern longitude limit of integration")
    parser.add_argument("-W", "--west", action="store",dest="west",type=float,
                         help="western longitude limit of integration")
    parser.add_argument("-L", "--levels", action="store",dest="levels",type=int,nargs=2,
                         help="start and end levels for integration")
    parser.add_argument("-X", "--factor", action="store",dest="factor",type=float,
                         help="scalar multiplicative factor")
    parser.add_argument("--gt", action="store",dest="gt",type=float,
                         help="only sum areas where the field takes values greater than this.")
    parser.add_argument("--lt", action="store",dest="lt",type=float,
                         help="only sum areas where the field takes values less than this.")

    args = parser.parse_args()

    geo_integrate(grid=args.grid,maskfile=args.maskfile,maskvar=args.maskvar,
                   file_in=args.file_in,var_in=args.var_in,file_out=args.file_out,
                   levels=args.levels,factor=args.factor,gt=args.gt,lt=args.lt,
                   north=args.north,south=args.south,east=args.east,west=args.west)

