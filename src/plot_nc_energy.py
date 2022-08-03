#! /usr/bin/env python
"""
    Script to calculate and plot timeseries of perturbation 
    potential energy and kinetic energy for the barotropic
    part of the solution from a model SSH and 3D velocity
    field. 

    Refer Gill (1982) 'Atmosphere-Ocean Dynamics', equation (5.7.4)

"""

__author__ = 'Dave Storkey (dave.storkey@metoffice.gov.uk)'
__date__ = '$Date: 2013/24/6 $'

import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import netCDF4

def plot_nc_energy(datafiles,field):

    # water depth
    H = 0.4
    # gravitational acceleration
    g = 9.8

    dataidT = netCDF4.Dataset(datafiles[0], mode='r')
    dataidU = netCDF4.Dataset(datafiles[1], mode='r')
    dataidV = netCDF4.Dataset(datafiles[2], mode='r')
    ssh = np.copy(dataidT.variables['sossheig'])
    u = np.copy(dataidU.variables['vozocrtx'])
    v = np.copy(dataidV.variables['vomecrty'])

    print 'ssh.shape : ',ssh.shape    
    pe = np.zeros(ssh.shape[0])

    # depth-mean currents (ASSUMES UNIFORM VERTICAL GRID SPACING)
    # don't include bottom layer
    u_bar = u[:,0:-1,:,:].mean(1)
    v_bar = v[:,0:-1,:,:].mean(1)

    print 'u.shape : ',u.shape    
    print 'u_bar.shape : ',u_bar.shape    
    jpi = u.shape[-1]
    jpj = u.shape[-2]
    ke = np.zeros(u.shape[0])

    # interpolate to set of T points excluding outermost row/column
    # which is land anyway in NEMO. 
    ssh_int = ssh[:,1:jpj-1,1:jpi-1]
    u_bar_int = 0.5 * ( u_bar[:,1:jpj-1,0:jpi-2] + u_bar[:,1:jpj-1,1:jpi-1] )  
    v_bar_int = 0.5 * ( v_bar[:,0:jpj-2,1:jpi-1] + v_bar[:,1:jpj-1,1:jpi-1] )  

    for time in range(ssh.shape[0]):
        pe[time] = g * np.sum( ssh_int[time,:,:]*ssh_int[time,:,:] )

    for time in range(u.shape[0]):
        ke[time] = H * np.sum( u_bar_int[time,:,:]*u_bar_int[time,:,:] 
            + v_bar_int[time,:,:]*v_bar_int[time,:,:] )

    if field == 'all':
        plt.plot(pe)
        plt.plot(ke)
        plt.plot(pe+ke)
        plt.legend(['PE','KE','total'],loc=1)
    elif field == 'peke':
        plt.plot(pe)
        plt.plot(ke)
        plt.legend(['PE','KE'],loc=3)
    elif field == 'total':
        plt.plot(pe+ke)
        plt.legend(['PE + KE'],loc=3)

    plt.xlim([0,600])
    plt.ylim([0,9])
    plt.xlabel('time (days)')
    plt.ylabel('perturbation energy (J/kg)')
    plt.show()

if __name__ == "__main__":
    plot_nc_energy(datafiles=[sys.argv[1],sys.argv[2],sys.argv[3]],field=sys.argv[4])
