#! /usr/bin/env python
'''
Plot timeseries of spatially-meaned KE.
(Expecting output from CDFTOOLS routines
cdfeke and cdfmean). 

DS. Jan 2022
'''

import argparse
import netCDF4 as nc
import matplotlib.pyplot as plt

def plot_ke(infiles=None,labels=None,time_units=None,outfile=None):

    if time_units is None:
        time_units="months"

    for infile in infiles:
        with nc.Dataset(infile,'r') as data1:
            mke=data1.variables['mean_3D_vomke_global'][:]
            eke=data1.variables['mean_3D_voeke_global'][:]

        plt.plot(mke[:,0,0])
        plt.plot(eke[:,0,0])

    if labels is not None:
        plt.legend(labels)
    plt.gca().set_xlabel('time ('+time_units+')')
    plt.gca().set_ylabel('energy (J/m3)')

    if outfile is None:
        plt.show()
    else:
        plt.gcf().savefig(outfile,dpi=200)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles",action="store",dest="infiles", nargs="+",
                    help="name(s) of input file(s).")
    parser.add_argument("-L", "--labels",action="store",dest="labels", nargs="+",
                    help="labels for legend.")
    parser.add_argument("-T", "--time_units", action="store",dest="time_units",
                    help="time units (default months)")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",
                    help="plot to output file - filetype by extension. If unset plot to GUI.")

    args = parser.parse_args()

    plot_ke(infiles=args.infiles,labels=args.labels,outfile=args.outfile,
            time_units=args.time_units)
