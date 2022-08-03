#! /usr/bin/env python
"""
    Script to calculate and plot timeseries of:
       - min/max of a field
       - variance of a field
"""

__author__ = 'Dave Storkey (dave.storkey@metoffice.gov.uk)'
__date__ = '$Date: 2013/24/6 $'

import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import netCDF4

def plot_nc_minmaxvar(datafile,field,func):

    dataid = netCDF4.Dataset(datafile, mode='r')
    fld = np.copy(dataid.variables[field])

    print 'fld.shape : ',fld.shape    
    functoplot = np.zeros(fld.shape[0])

    if func == 'min':
        functoplot = fld.min(1).min(1)
        plt.ylabel('SSH amplitude (m)')
    elif func == 'max':
        functoplot = fld.max(1).max(1)
        plt.ylabel('SSH amplitude (m)')
    elif func == 'minmax':
        functoplot = fld.min(1).min(1)
        functoplot2 = fld.max(1).max(1)
        plt.ylabel('SSH amplitude (m)')
    elif func == 'var':
        for time in range(fld.shape[0]):
            functoplot[time] = sp.var(fld[time,:,:])
        plt.ylabel('SSH variance (m)')

    plt.xlabel('time (days)')
    plt.plot(functoplot)
    if func == 'minmax':
        plt.plot(functoplot2)
        plt.legend( ('SSH minimum','SSH maximum'), loc=2 )
    plt.xlim( [0,450] )
    plt.ylim( [-0.2,+0.2] )
    plt.show()

if __name__ == "__main__":
    plot_nc_minmaxvar(datafile=sys.argv[1],field=sys.argv[2],func=sys.argv[3])
