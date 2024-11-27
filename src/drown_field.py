#! /usr/bin/env python

'''
Script to apply "drowning", ie. extrapolating ocean values over
land points horizontally to avoid problems with eg. climatologies
when you change the land-sea mask.

@author: Dave Storkey
@date: July 2022
'''

import xarray as xr
import numpy as np
import re

def drown_field(file_in=None, file_out=None, varnames=None, vert_fill=None,
                niter=None, min_points=None, new_mask_file=None, new_mask_name=None):

    if niter is None:
        niter = 10

    if min_points is None:
        min_points = 1

    if varnames is None:
        raise Exception("Error: must specify at least one field to drown.")
    elif not isinstance(varnames,list):
        varnames=[varnames]

    fields=[]
    with xr.open_dataset(file_in) as indata:
        for timedimname in ['t','time','time_counter','time_centered']:
            try:
                timedim = getattr(indata,timedimname)
            except(AttributeError):
                pass
            else:
                break
        else: 
            timedim = None
        try:
            nav_lat = indata.nav_lat
        except(AttributeError):
            nav_lat = None
        try:
            nav_lon = indata.nav_lon
        except(AttributeError):
            nav_lon = None

        for varname in varnames:
            fields.append(getattr(indata,varname))

    new_mask = None
    if new_mask_file is not None:
        if new_mask_name is None:
            new_mask_name="tmask"
        with xr.open_dataset(new_mask_file) as new_mask_in:
            new_mask = getattr(new_mask_in,new_mask_name)

    for field in fields:
        drowned_field = field
        if vert_fill:
            # just fill in from points immediately above
            nn = 0
            while nn < niter:
                nn += 1
                field_zeros = xr.where( np.isnan(drowned_field), 0.0, drowned_field)
                field_up = np.roll(field_zeros.values,shift=-1,axis=-3)
                drowned_field.values = xr.where( np.isnan(drowned_field.values),
                                           field_up, drowned_field )
        else:
            # do horizontal averaging of neighbouring non-NaN points
            nn = 0
            while nn < niter:
                nn += 1
                field_zeros = xr.where( np.isnan(drowned_field), 0.0, drowned_field)
                field_n = np.roll(field_zeros.values,shift=-1,axis=-2)
                field_s = np.roll(field_zeros.values,shift=+1,axis=-2)
                field_e = np.roll(field_zeros.values,shift=-1,axis=-1)
                field_w = np.roll(field_zeros.values,shift=+1,axis=-1)
                field_ne = np.roll(field_zeros.values,shift=(-1,-1),axis=(-2,-1))
                field_nw = np.roll(field_zeros.values,shift=(-1,+1),axis=(-2,-1))
                field_se = np.roll(field_zeros.values,shift=(+1,-1),axis=(-2,-1))
                field_sw = np.roll(field_zeros.values,shift=(+1,+1),axis=(-2,-1))
                weights = np.sign(field_n)*np.sign(field_n) + \
                          np.sign(field_s)*np.sign(field_s) + \
                          np.sign(field_e)*np.sign(field_e) + \
                          np.sign(field_w)*np.sign(field_w) + \
                          np.sign(field_ne)*np.sign(field_ne) + \
                          np.sign(field_nw)*np.sign(field_nw) + \
                          np.sign(field_se)*np.sign(field_se) + \
                          np.sign(field_sw)*np.sign(field_sw) 
                drowned_field.values = xr.where( np.isnan(drowned_field.values) & (weights > min_points-1),
                            ( field_n  + field_s  + field_e  + field_w +
                              field_nw + field_ne + field_sw + field_se ) / weights, drowned_field)

        if field.name == fields[0].name:
            outdata = drowned_field.to_dataset()    
        else:
            outdata[field.name] = drowned_field

    if new_mask is not None:
        n_fail = np.count_nonzero( (np.isnan(drowned_field[0].values) & new_mask == 1) )
        if n_fail > 0:
            print(n_fail," NaN points at sea point in new mask.")
            print(np.where( (np.isnan(drowned_field[0].values) & new_mask == 1) ))
        else:
            print("No more NaNs at sea points in new mask.")

    if nav_lat is not None and nav_lon is not None:
        outdata.update({'nav_lat':nav_lon ,
                        'nav_lat':nav_lon })

    print('timedimname: ',timedimname)
    if timedim is not None:
        outdata.to_netcdf(file_out,unlimited_dims=timedimname)
    else:
        outdata.to_netcdf(file_out)

                
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--file_in", action="store",dest="file_in",
                         help="input file")
    parser.add_argument("-o", "--file_out", action="store",dest="file_out",
                         help="output file")
    parser.add_argument("-v", "--varnames", action="store",dest="varnames",nargs="+",
                         help="name of field(s) to drown")
    parser.add_argument("-V", "--vert_fill", action="store_true",dest="vert_fill",
                         help="do vertical filling instead of horizontal filling")
    parser.add_argument("-n", "--niter", action="store",dest="niter",type=int,
                         help="number of iterations - default 10")
    parser.add_argument("-m", "--min_points", action="store",dest="min_points",type=int,
                         help="minimum number of nonland neighbours for filling to occur - default 1")
    parser.add_argument("-K", "--new_mask_file", action="store",dest="new_mask_file",
                         help="name of file with new mask to check drowned field against")
    parser.add_argument("-k", "--new_mask_name", action="store",dest="new_mask_name",
                         help="new mask field to check drowned field against")

    args = parser.parse_args()

    drown_field(file_in=args.file_in,file_out=args.file_out,varnames=args.varnames,vert_fill=args.vert_fill,
                niter=args.niter,min_points=args.min_points, new_mask_file=args.new_mask_file,
                new_mask_name=args.new_mask_name )

