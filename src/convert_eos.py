#! /usr/bin/env python

'''
Routine to convert between eos80 and teos10 thermodynamic variables. 
Based on routines by Catherine Guiavarc'h

@author: Dave Storkey
@date: Sep 2022
'''

import xarray as xr
import numpy as np
import gsw as gsw

def convert_eos(infile=None,outfile=None,eos80=None):

    if eos80:
       tem_in_name = ["thetao_con"]
       sal_in_name = ["so_abs"]
       tem_out_name = "thetao_pot"
       sal_out_name = "so_pra"
    else:
       tem_in_name = ["thetao_pot","thetao","votemper"]
       sal_in_name = ["so_pra","so","vosaline"]
       tem_out_name = "thetao_con"
       sal_out_name = "so_abs"

    with xr.open_dataset(infile) as indata:
        for timedimname in ['t','time','time_counter','time_centered']:
            if timedimname in indata.dims.keys():
                break
            else:
                pass
        else: 
            timedim = None
        for lonname in ["longitude", "lon", "nav_lon"]:
            try:
                lon = getattr(indata,lonname)
            except(AttributeError):
                pass
            else:
                break
        else:
            raise Exception("Error: could not find longitude coordinate in file.")
        for latname in ["latitude", "lat", "nav_lat"]:
            try:
                lat = getattr(indata,latname)
            except(AttributeError):
                pass
            else:
                break
        else:
            raise Exception("Error: could not find latitude coordinate in file.")
        for depname in ["depth","deptht"]:
            try:
                depth = getattr(indata,depname)
            except(AttributeError):
                pass
            else:
                break
        else:
            raise Exception("Error: could not find depth coordinate in input file.")
        for temname in tem_in_name:
            try:
                tem_in = getattr(indata,temname)
            except(AttributeError):
                pass
            else:
                break
        else:
            raise Exception("Error: could not find temperature field in input file.")
        for salname in sal_in_name:
            try:
                sal_in = getattr(indata,salname)
            except(AttributeError):
                pass
            else:
                break
        else:
            raise Exception("Error: could not find salinity field in input file.")

    tem_in, lonb = xr.broadcast(tem_in, lon)
    tem_in, latb = xr.broadcast(tem_in, lat)
    depthb, tem_in = xr.broadcast(depth, tem_in)
    print("tem_in.shape :",tem_in.shape)
    print("sal_in.shape :",sal_in.shape)
    print("lonb.shape :",lonb.shape)
    print("latb.shape :",latb.shape)
    print("depthb.shape :",depthb.shape)
    pressure = gsw.p_from_z(-depthb.values, latb.values)
    print("pressure.shape :",pressure.shape)

    sal_out = sal_in.copy()
    sal_out.name = sal_out_name
    tem_out = tem_in.copy()
    tem_out.name = tem_out_name

    if eos80:
        # convert to potential temperature and practical salinity from
        # Conservative Temperature and Absolute Salinity.
        sal_out.values = gsw.SP_from_SA(sal_in.values, pressure, lonb.values, latb.values)
        tem_out.values = gsw.pt_from_CT(sal_in.values,tem_in.values)
        tem_out.attrs["standard_name"]="sea_water_potential_temperature"
        tem_out.attrs["long_name"]="Potential Temperature"
        sal_out.attrs["standard_name"]="sea_water_salinity"
        sal_out.attrs["long_name"]="Practical Salinity"
    else:
        # convert to Conservative Temperature and Absolute Salinity from
        # potential temperature and practical salinity.
        sal_out.values = gsw.SA_from_SP(sal_in.values, pressure, lonb.values, latb.values)
        tem_out.values = gsw.CT_from_pt(sal_in.values,tem_in.values)
        tem_out.attrs["standard_name"]="sea_water_conservative_temperature"
        tem_out.attrs["long_name"]="Conservative Temperature"
        sal_out.attrs["standard_name"]="sea_water_absolute_salinity"
        sal_out.attrs["long_name"]="Absolute Salinity"

#    #Create a mask
#       tmsk = np.ones(shape=SA.shape)
#       tmsk[SA>100.] = 0.
#    # Masking land
#       SP[tmsk==0] = 1.e+20 #np.nan
#       PT[tmsk==0] = 1.e+20 #np.nan

    outdata = tem_out.to_dataset()
    outdata[sal_out_name] = sal_out
    outdata[lonname] = lon
    outdata[latname] = lat
    outdata[depname] = depth

    if timedim is not None:
        outdata.to_netcdf(outfile,unlimited_dims=timedimname)
    else:
        outdata.to_netcdf(outfile)

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", action="store",dest="infile",
                         help="input file")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",
                         help="output file")
    parser.add_argument("-E", "--eos80", action="store_true",dest="eos80",
                         help="convert to EOS80 variables from TEOS-10 variables")

    args = parser.parse_args()

    convert_eos(infile=args.infile,outfile=args.outfile,eos80=args.eos80)
