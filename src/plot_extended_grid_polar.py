#! /usr/bin/env python
"""
    Script to do a polar stereographic plot of 
    extended ORCA grid coordinates.
"""

__author__ = 'Dave Storkey (dave.storkey@metoffice.gov.uk)'
__date__ = '$Date: December 2013 $'

import sys
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.basemap import Basemap
import netCDF4
from mpl_toolkits.mplot3d import Axes3D

def plot_coords_grid(coordsfile,map,spacing=None,colors='black'):

    coords_id = netCDF4.Dataset(coordsfile, mode='r')
    lam = np.squeeze(np.copy(coords_id.variables['glamt']))
    phi = np.squeeze(np.copy(coords_id.variables['gphit']))

    if spacing is None:
        spacing = 10

    i = np.arange(lam.shape[1])
    ncont_i = len(i)/spacing
    i2d = np.ones(lam.shape)
    i2d = i * i2d  # broadcasting

    j = np.arange(lam.shape[0])
    if spacing > 3:
        # hack to make sure I have a contour at j=0:
        j[1] = -spacing
    ncont_j = len(j)/spacing
    j2d = np.ones(lam.shape).transpose()
    j2d = j * j2d # broadcasting
    j2d = j2d.transpose()
 
    [lon,lat] = map(lam,phi)
    print 'ncont_i: ',ncont_i
    map.contour(lon, lat, i2d, ncont_i, colors=colors)
    map.contour(lon, lat, j2d, ncont_j, colors=colors)


def plot_polar_coords(coordsfile1,bathyfile=None,bathyfile2=None,coordsfile2=None,spacing=10):

    a = 1.0

    map = Basemap(projection='spstere',boundinglat=-45,lon_0=180,resolution='h')
    map.drawcoastlines()
#    map.drawmeridians(np.arange(-360, 360, 30))
#    map.drawparallels(np.arange(-90, 90, 30))
    cmap = cm.jet
    cmap.set_bad('white')

    if bathyfile is not None:
        bathy_id = netCDF4.Dataset(bathyfile, mode='r')
        for lonname in ['glamt', 'longitude', 'lon']:
            try:
                lambat = np.squeeze(np.copy(bathy_id.variables[lonname]))
            except KeyError:
                pass
            else:
                break
        else:
            raise KeyError
        for latname in ['gphit', 'latitude', 'lat']:
            try:
                phibat = np.squeeze(np.copy(bathy_id.variables[latname]))
            except KeyError:
                pass
            else:
                break
        else:
            raise KeyError
        bathy = np.squeeze(np.copy(bathy_id.variables['Bathymetry']))
        bathy_masked = ma.masked_where(bathy == 0.0, bathy)
        values = 1.0 - bathy_masked / 5000.0

        [lon,lat] = map(lambat,phibat)
        print lon.min(),lon.max()
        print lat.min(),lat.max()
        map.contourf(lon, lat, values, 30)

    if bathyfile2 is not None:
        bathy2_id = netCDF4.Dataset(bathyfile2, mode='r')
        for lonname in ['glamt', 'longitude', 'lon']:
            try:
                lambat2 = np.squeeze(np.copy(bathy2_id.variables[lonname]))
            except KeyError:
                pass
            else:
                break
        else:
            raise KeyError
        for latname in ['gphit', 'latitude', 'lat']:
            try:
                phibat2 = np.squeeze(np.copy(bathy2_id.variables[latname]))
            except KeyError:
                pass
            else:
                break
        else:
            raise KeyError
        bathy2var = bathy2_id.variables['Bathymetry']
        bathy2 = ma.squeeze(ma.copy(bathy2var))
        bathy2_masked = ma.masked_where(( bathy2 == 0.0 ) | (bathy2 == bathy2var._FillValue), bathy2)
        values2 = bathy2_masked
       
        [lon,lat] = map(lambat2,phibat2)
        print lon.min(),lon.max()
        print lat.min(),lat.max()
        print 'values2 min/max :',values2.min(),values2.max()
        map.contourf(lon, lat, values2, 30)

    plot_coords_grid(coordsfile1,map,spacing=spacing,colors='black')

    if coordsfile2 is not None:
        plot_coords_grid(coordsfile2,map,spacing=spacing,colors='blue')

    plt.show()


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("coordsfile", help="name of coorddinates file")
    parser.add_argument("-b", "--bathyfile", action="store",dest="bathyfile",
                    default=None,help="bathymetry file (optional)")
    parser.add_argument("-B", "--bathyfile2", action="store",dest="bathyfile2",
                    default=None,help="second bathymetry file (optional)")
    parser.add_argument("-c", "--coordsfile2", action="store",dest="coordsfile2",
                    default=None,help="second coordinates file (optional)")
    parser.add_argument("-s", "--spacing", action="store",type=int,dest="spacing",
                    default=10,help="spacing for coordinate line plotting")

    args = parser.parse_args()

    plot_polar_coords(args.coordsfile,bathyfile=args.bathyfile,bathyfile2=args.bathyfile2,
                      coordsfile2=args.coordsfile2,spacing=args.spacing)

