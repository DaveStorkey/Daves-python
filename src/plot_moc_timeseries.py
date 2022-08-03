#! /usr/bin/env python
'''
Routine to plot MOC timeseries at a particular latitude from a set of model runs.
Takes as input the timeseries of "moc.nc" files for each run produced by Ocean Assess when
it calculates the MOC. 

TO DO: remove hardcoding for Atlantic MOC. Option to plot any of MOC fields. 

Created June 2017

@author: Dave Storkey
'''

import glob
import netCDF4 as nc
import numpy as np
import model
from nemo_moc_tools import moc_max
import matplotlib
import matplotlib.pyplot as plt

def plot_moc_timeseries(dirnames=None, titles=None, lat=None, outfile=None, time_offset=None,
                        legend_loc=None, noshow=False, xtitle=None, ytitle=None):

    if len(dirnames) != len(titles):
        raise Exception('Error: must have the same number of input directories and titles.')

    if legend_loc is None:
        legend_loc = 1

    moc_ts = {}
    for (dirname,title) in zip(dirnames,titles):
        moc_ts[title]=[]
        filelist = glob.glob(dirname+'/*moc.nc')
        filelist.sort()
        for filename in filelist:
            dataset = nc.Dataset(filename,'r')
            mocfield = dataset.variables['zomsfatl'][:]
            for latname in ['nav_lat','nav_lat_grid_V']:
                try:
                    latfield = dataset.variables[latname][:]
                except KeyError:
                    pass
                else:
                    break
            else:
                raise Exception('Could not find latitude field in file '+filename)
            depthfield = dataset.variables['depthw'][:]
            moc_ts[title].append(moc_max(mocfield,lat=latfield,depth=depthfield,lat_range=lat)[0])

    len_time_axis = 0
    for title in titles:
        len_time_axis = max(len_time_axis,len(moc_ts[title]))        

    time_axis = np.array(range(len_time_axis))
    if time_offset is not None:
        time_axis[:] = time_axis[:] + time_offset  

    for title in titles:
        print title+':',moc_ts[title]
        plt.plot(time_axis[0:len(moc_ts[title])],moc_ts[title])     

    plt.legend(titles,loc=legend_loc,fontsize='medium')
    if xtitle is None:
        xtitle='year of spin up'
    plt.gca().set_xlabel(xtitle)
    if ytitle is None:
        ytitle='maximum overturning (Sv) at latitude '+str(lat)
    plt.gca().set_ylabel(ytitle)

    if outfile is not None:
        matplotlib.rcParams['font.size'] =8
        plt.savefig(outfile,dpi=200)
    else:
        if not noshow:
            plt.show()


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dir", action="store",dest="dirnames", nargs="+",
                        help="name of input directories (one for each model run)")
    parser.add_argument("-t", "--titles", action="store",dest="titles", nargs="+",
                        help="titles (one for each model run)")
    parser.add_argument("-x", "--xtitle", action="store",dest="xtitle", nargs="+",
                        help="title for x axis")
    parser.add_argument("-y", "--ytitle", action="store",dest="ytitle", nargs="+",
                        help="title for y axis")
    parser.add_argument("-l", "--lat", action="store",dest="lat",type=float,
                        help="latitude of MOC to plot")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",
                        help="output file, format by extension, if unset plot to GUI")
    parser.add_argument("-n", "--noshow", action="store_true",dest="noshow",
                        help="don't run plt.show()")
    parser.add_argument("-T", "--time_offset", action="store",dest="time_offset",type=float,default=None,
                        help="offset to add to time axis")
    parser.add_argument("-L", "--legend_loc", action="store",dest="legend_loc",type=int,default=None,
                        help="location for legend (see matplotlib doc)")

    args = parser.parse_args()

    plot_moc_timeseries(dirnames=args.dirnames, titles=args.titles, lat=args.lat, outfile=args.outfile,
                        time_offset=args.time_offset, legend_loc=args.legend_loc, xtitle=args.xtitle, 
                        ytitle=args.ytitle, noshow=args.noshow )
