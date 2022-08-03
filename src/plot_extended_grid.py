#! /usr/bin/env python
"""
    Script to do a 3D plot of extended ORCA grid coordinates.
"""

__author__ = 'Dave Storkey (dave.storkey@metoffice.gov.uk)'
__date__ = '$Date: December 2013 $'

import sys
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib import cm
import netCDF4
from mpl_toolkits.mplot3d import Axes3D

def plot_3D_coords(extcoordsfile,bathyfile=None):

    a = 1.0
    deg2rad = np.pi/180.0

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='polar')
    ax.set_axis_off()
    cmap = cm.jet
    cmap.set_bad('white')

    if bathyfile is not None:
        bathy_id = netCDF4.Dataset(bathyfile, mode='r')
        lambat = np.squeeze(np.copy(bathy_id.variables['glamt']))
        phibat = np.squeeze(np.copy(bathy_id.variables['gphit']))
        bathy = np.squeeze(np.copy(bathy_id.variables['Bathymetry']))
        bathy_masked = ma.masked_where(bathy == 0.0, bathy)

        print 'minimum phibat :',phibat.min()

        ind_south = np.where(phibat < 0.0)
        j_equator = ind_south[0].max() + 1
        print 'j_equator : ',j_equator

        lambat_south = lambat[0:j_equator-1,:]
        phibat_south = phibat[0:j_equator-1,:]
        bathy_south = bathy_masked[0:j_equator-1,:]

        Xbat = a * np.cos(deg2rad*phibat_south) * np.cos(deg2rad*lambat_south)
        Ybat = a * np.cos(deg2rad*phibat_south) * np.sin(deg2rad*lambat_south)
        Zbat = a * np.sin(deg2rad*phibat_south)
        Nbat = 1.0 - bathy_south / 5000.0

        surf = ax.plot_surface(Xbat, Ybat, Zbat, rstride=4, cstride=4, facecolors=cmap(Nbat), 
                               alpha=1, linewidth=0.4, shade=False)

    extcoords_id = netCDF4.Dataset(extcoordsfile, mode='r')
    lam = np.squeeze(np.copy(extcoords_id.variables['glamt']))
    phi = np.squeeze(np.copy(extcoords_id.variables['gphit']))

    ind_south = np.where(phi < 0.0)
    j_equator = ind_south[0].max() + 1
    print 'j_equator : ',j_equator

    lam_south = lam[0:j_equator-1,:]
    phi_south = phi[0:j_equator-1,:]

    X = a * np.cos(deg2rad*phi_south) * np.cos(deg2rad*lam_south)
    Y = a * np.cos(deg2rad*phi_south) * np.sin(deg2rad*lam_south)
    Z = a * np.sin(deg2rad*phi_south)

    ax.plot_wireframe(X, Y, Z,  rstride=15, cstride=15, color='black')

    ax.view_init(-90, 0) 
    plt.show()


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("coordsfile", help="name of coorddinates file")
    parser.add_argument("-b", "--bathyfile", action="store",dest="bathyfile",
                    default=None,help="bathymetry file (optional)")

    args = parser.parse_args()

    plot_3D_coords(args.coordsfile,bathyfile=args.bathyfile)

