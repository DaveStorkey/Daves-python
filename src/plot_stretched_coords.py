#! /usr/bin/env python
"""
    Script to do a 3D plot of ORCA grid coordinates.
"""

__author__ = 'Dave Storkey (dave.storkey@metoffice.gov.uk)'
__date__ = '$Date: 2013/12/2 $'

import sys
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib import cm
import netCDF4
from mpl_toolkits.mplot3d import Axes3D

def plot_3D_coords(coordsfile,bathyfile,viewpoint='north',lat_cut=None,scale_factor=None):

    deg2rad = np.pi/180.0

    coords_id = netCDF4.Dataset(coordsfile, mode='r')
    bathy_id = netCDF4.Dataset(bathyfile, mode='r')
    lam = np.squeeze(np.copy(coords_id.variables['glamt']))
    phi = np.squeeze(np.copy(coords_id.variables['gphit']))
    bathy = np.squeeze(np.copy(bathy_id.variables['Bathymetry']))
    bathy_masked = ma.masked_where(bathy == 0.0, bathy)

    a = 1.0

    print 'minimum phi :',phi.min()

    ind_south = np.where(phi < 0.0)
    j_equator = ind_south[0].max() + 1
    print 'j_equator : ',j_equator

    if lat_cut is not None:
        ind_cut = np.where(phi >= lat_cut)
        j_cut = ind_cut[0].min()
        print 'j_cut : ',j_cut

    lam_south = lam[0:j_equator-1,:]
    phi_south = phi[0:j_equator-1,:]
    bathy_south = bathy_masked[0:j_equator-1,:]

    lam_north = lam[j_equator-1:-1,:]
    phi_north = phi[j_equator-1:-1,:]
    bathy_north = bathy_masked[j_equator-1:-1,:]

    lam_cut = lam[j_cut-1:-1,:]
    phi_cut = phi[j_cut-1:-1,:]

    if scale_factor is not None:
        phi_cut = 90.0 - ( (90.0-phi_cut)*scale_factor )

    if viewpoint == 'north':
        X = a * np.cos(deg2rad*phi_north) * np.cos(deg2rad*lam_north)
        Y = a * np.cos(deg2rad*phi_north) * np.sin(deg2rad*lam_north)
        Z = a * np.sin(deg2rad*phi_north)
        N = 1.0 - bathy_north / 5000.0
    elif viewpoint == 'south':
        X = a * np.cos(deg2rad*phi_south) * np.cos(deg2rad*lam_south)
        Y = a * np.cos(deg2rad*phi_south) * np.sin(deg2rad*lam_south)
        Z = a * np.sin(deg2rad*phi_south)
        N = 1.0 - bathy_south / 5000.0

    if lat_cut is not None:
        X_cut = a * np.cos(deg2rad*phi_cut) * np.cos(deg2rad*lam_cut)
        Y_cut = a * np.cos(deg2rad*phi_cut) * np.sin(deg2rad*lam_cut)
        Z_cut = a * np.sin(deg2rad*phi_cut)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d',frame_on=False)
    ax.set_axis_off()
    cmap = cm.jet
    cmap.set_bad('white')

    surf = ax.plot_surface(X, Y, Z, rstride=4, cstride=4, facecolors=cmap(N), alpha=1, linewidth=0.4, shade=False)
    if lat_cut is not None:
        ax.plot_wireframe(X_cut, Y_cut, Z_cut,  rstride=15, cstride=15, color='black')
    else:
        ax.plot_wireframe(X, Y, Z,  rstride=15, cstride=15, color='black')

    if viewpoint == 'north' or viewpoint == '20N':
        ax.view_init(90, 0) 
    elif viewpoint == 'south':
        ax.view_init(-90, 0) 
    plt.show()


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("coordsfile", help="name of coorddinates file")
    parser.add_argument("bathyfile", help="name of bathymetry file")
    parser.add_argument("-N", "--north", action="store_const",dest="viewpoint",const="north",
                    default="north",help="set initial viewpoint to north")
    parser.add_argument("-S", "--south", action="store_const",dest="viewpoint",const="south",
                    default="north",help="set initial viewpoint to south")
    parser.add_argument("-c", "--lat_cut", action="store",type=float,dest="lat_cut",
                    default=None,help="set cut off latitude")
    parser.add_argument("-s", "--scale_factor", action="store",type=float,dest="scale_factor",
                    default=None,help="scaling for latitudes")

    args = parser.parse_args()

    plot_3D_coords(args.coordsfile,args.bathyfile,viewpoint=args.viewpoint,lat_cut=args.lat_cut,
                   scale_factor=args.scale_factor)        

