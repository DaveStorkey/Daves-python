#! /usr/bin/env python
"""
    Script to extract a subsampled, subsea bathymetry under
    the ice shelves from BEDMAP2
"""

__author__ = 'Dave Storkey (dave.storkey@metoffice.gov.uk)'
__date__ = '$Date: December 2013 $'

import numpy as np
import netCDF4

def make_sub_iceshelf_bathy(bedmap,outfile,subsampling=4):

    bedmap_file = netCDF4.Dataset(bedmap, mode='r')

    bathy = bedmap_file.variables['Bathymetry']  
    isf_draft = bedmap_file.variables['isf_draft']

    lon = bedmap_file.variables['lon']
    lat = bedmap_file.variables['lat']

# subsampled arrays in memory:

    bathy_subs = bathy[0:-1:subsampling,0:-1:subsampling]
    draft_subs = isf_draft[0:-1:subsampling,0:-1:subsampling]

    lon_subs = lon[0:-1:subsampling,0:-1:subsampling]
    lat_subs = lat[0:-1:subsampling,0:-1:subsampling]

    print 'bathy_subs.shape : ',bathy_subs.shape
    print 'draft_subs.shape : ',draft_subs.shape

    print 'bathy_subs.min/max : ',bathy_subs.min(),bathy_subs.max()
    print 'draft_subs.min/max : ',draft_subs.min(),draft_subs.max()

    bathy_out = np.where( ( bathy_subs < 0.0 ) & ( abs(draft_subs) < abs(bathy_subs) ),bathy_subs,0.0)

    output_file = netCDF4.Dataset(outfile, mode='w')
    output_file.createDimension('y',bathy_out.shape[0])
    output_file.createDimension('x',bathy_out.shape[1])
    output_file.createVariable('Bathymetry',bathy.dtype,dimensions=('y','x'),fill_value=bathy.missing_value)
    output_file.createVariable('longitude',bathy.dtype,dimensions=('y','x'))
    output_file.createVariable('latitude',bathy.dtype,dimensions=('y','x'))

    bathyvar_out = output_file.variables['Bathymetry']
#    bathyvar_out._FillValue = bathy.missing_value
    bathyvar_out[:,:] = bathy_out[:,:]

    lon_out = output_file.variables['longitude']
    lon_out[:,:] = lon_subs[:,:]

    lat_out = output_file.variables['latitude']
    lat_out[:,:] = lat_subs[:,:]

    bedmap_file.close()
    output_file.close()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("bedmap", help="name of BEDMAP file")
    parser.add_argument("outfile", help="name of output file")
    parser.add_argument("-s", "--subsample", action="store",type=int,dest="subsampling",
                    default=4,help="subsampling frequency (in x and y) for BEDMAP data")

    args = parser.parse_args()

    make_sub_iceshelf_bathy(args.bedmap,args.outfile,subsampling=args.subsampling)

