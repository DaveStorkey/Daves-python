#! /usr/bin/env python

'''
Cut out a domain from a netcdf file
based on min/max lats/lons.

(Originally called cutout_GEBCO.py).

@author: Dave Storkey
@date: Oct 2021
'''

import xarray as xr
import numpy as np

def cutout_domain(infile=None,invar=None,outfile=None,
                 south=None,north=None,west=None,east=None):

    if invar is None:
        invar="elevation"

    bathy_in = xr.open_dataset(infile)
    elev = bathy_in.get(invar)
    lats = elev.lat.values
    lons = elev.lon.values

    if south is None:
        south = np.min(lats)
    if north is None:
        north = np.max(lats)
    if west is None:
        west = np.min(lons)
    if east is None:
        east = np.max(lons)

    print("west, east before check : ",west,east)

    # Restrict west and east to the range [-360,360]
    # and make sure that east > west within this range.
    # NB. python modulus function (%) doesn't behave the
    # same as mathematical modulus for negative numbers!
    west = ((west*np.sign(west))%360.0)*np.sign(west)
    east = ((east*np.sign(east))%360.0)*np.sign(east)
    if east < west:
        if west < 0.0:
            east+=360.0
        else:
            west-=360.0
    print("west, east after check : ",west,east)
    # Fix longitudes to appropriate range depending on the values of [west,east]
    if west < -180.0:
        # use range [-360,0]
        lons[:] = lons[:] - np.maximum(0.0,np.sign(lons[:]))*360.0
        print("Using range [-360,0] : ",np.min(lons),np.max(lons))
    elif east > 180.0:
        # use range [0,360]
        lons[:] = lons[:] - np.minimum(0.0,np.sign(lons[:]))*360.0
        print("Using range [0,+360] : ",np.min(lons),np.max(lons))
    else:
        # use range [-180,+180]
        lons[:] = lons[:] - np.maximum(0.0,np.sign(lons[:]-180.0))*360.0
        lons[:] = lons[:] - np.minimum(0.0,np.sign(lons[:]+180.0))*360.0
        print("Using range [-180,+180] : ",np.min(lons),np.max(lons))
    
    latsel = lats[np.where((lats >= south) & (lats <= north))]
    lonsel = lons[np.where((lons >= west ) & (lons <= east ))] 

    elev_cutout = elev.sel(lat=latsel,lon=lonsel)
    elev_cutout.values = np.maximum(0.0, -1.0 * elev_cutout.values)

    out_dataset = xr.Dataset(
        {
            "elevation" : elev_cutout.load()
        },
        coords={
            "lon"      : (('lon'),lonsel,{'standard_name':'longitude','units':'degrees_east','_FillValue':False}),
            "lat"      : (('lat'),latsel,{'standard_name':'latitude','units':'degrees_north','_FillValue':False}),
        },
    )
    
    out_dataset.to_netcdf(outfile)

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", action="store",dest="infile",
                         help="name of input file")
    parser.add_argument("-v", "--invar", action="store",dest="invar",
                         help="name of input variable")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",
                         help="name of output file")
    parser.add_argument("-W", "--west", action="store",dest="west",type=float,default=None,
                    help="western limit of area to plot")
    parser.add_argument("-E", "--east", action="store",dest="east",type=float,default=None,
                    help="eastern limit of area to plot")
    parser.add_argument("-S", "--south", action="store",dest="south",type=float,default=None,
                    help="southern limit of area to plot")
    parser.add_argument("-N", "--north", action="store",dest="north",type=float,default=None,
                    help="northern limit of area to plot")

    args = parser.parse_args()
 
    cutout_domain(infile=args.infile,invar=args.invar,outfile=args.outfile,
                 south=args.south,north=args.north,west=args.west,east=args.east)
