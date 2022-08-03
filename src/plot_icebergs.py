#! /usr/bin/env python
'''
Routine to plot positions of icebergs from icebergs trajectory file
or icebergs restart file.

Currently just plots the first position of each iceberg in the file
with the choice of northern or southern stereopolar projection. 

Created Jan 2018

Dec 2019 : Update to python 3. DS.

@author: Dave Storkey
'''

import netCDF4 as nc
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import textwrap

def plot_icebergs(filename,outfile=None,proj=None,marker=None,size=None,rec=None,title=None):

    if rec is None:
        rec=0

    with nc.Dataset(filename,'r') as file_in:
        lons = file_in.variables['lon'][:]
        lats = file_in.variables['lat'][:]
        ntraj,ntime = lons.shape
        print ('number of trajectories: ',ntraj)
        print ('number of times: ',ntime)

    points_to_plot = zip(lons[:,rec],lats[:,rec])

    plt.figure(1)
    if proj == 'northps':
        map = Basemap(projection='npstere',boundinglat=50.0,lon_0=0.0,resolution='h')
    elif proj == 'southps':
        map = Basemap(projection='spstere',boundinglat=-50.0,lon_0=180.0,resolution='h')
    else:
        raise Exception('Unrecognised map projection.')

    map.drawlsmask()
    for point in points_to_plot:
        x, y = map(point[0], point[1])
        map.plot(x,y,marker='o',color='Blue',markersize=2)

    if title is not None:
        title=textwrap.fill(title,70)
        plt.gcf().suptitle(title, fontsize=12, y=0.95)    

    if outfile is not None:
        matplotlib.rcParams['font.size'] =8
        plt.savefig(outfile,dpi=200)
    else:
        plt.show()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="name of input file")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",default=None,
                    help="plot to output file - filetype by extension. If unset plot to GUI.")
    parser.add_argument("-m", "--marker", action="store",dest="marker",default='o',
                    help="symbol to plot for each iceberg")
    parser.add_argument("-p", "--proj", action="store",dest="proj",default='northps',
                    help="projection: northps, southps")
    parser.add_argument("-r", "--rec", action="store",dest="rec",type=int,default=0,
                    help="record number to plot (defaults to zero)")
    parser.add_argument("-s", "--size", action="store",dest="size",type=int,default=2,
                    help="size of symbol to plot for each iceberg")
    parser.add_argument("-t", "--title", action="store",dest="title",default=None,
                    help="title for plot")

    args = parser.parse_args()

    plot_icebergs(args.filename,proj=args.proj,outfile=args.outfile,marker=args.marker,size=args.size,
                  rec=args.rec,title=args.title)        
