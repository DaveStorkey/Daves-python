#! /usr/bin/env python
"""
Script to calculate and plot timeseries of a volume-mean of a field. 
"""

__author__ = 'Dave Storkey (dave.storkey@metoffice.gov.uk)'
__date__ = '$Date: July 2017 $'

import socket
import matplotlib
if 'spice' in socket.gethostname():
    # Note this disables plt.show()
    matplotlib.use('Agg')
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import iris
import netCDF4 as nc
import general_tools as gt

def volmean_timeseries(files=None,fieldname=None,meshfile=None,region=None,firstlevel=None,lastlevel=None,mask=None):
    # Based on some of Tim Graham's code.
    # TO DO:  1. check that field is masked (perhaps separate routine?)
    #         2. add option to import a mask
    #         3. select volume based on depths in metres rather than levels
    files.sort()
    dummy_cube = gt.read_cube(files[0],fieldname)
    lon=dummy_cube.coord('longitude')
    lat=dummy_cube.coord('latitude')
    region_mask=gt.mask_gen(lon.points,lat.points,region)

    dummy_cube.data.mask=gt.combine_masks(dummy_cube.data.mask,region_mask)
    dummy_cube,cube_index=gt.trim_cube(dummy_cube,dims=[-2,-1])
    if firstlevel is not None and lastlevel is not None:
        # if levels are not specified this is already initialised 
        # to slice(None) which selects all levels.
        cube_index[-3]=np.arange(lastlevel-firstlevel+1)+firstlevel
    dummy_cube = gt.read_cube(files[0],fieldname)
    
    ncid=nc.Dataset(meshfile)
    for e3tvar in ['e3t','e3t_0']:
        try:
            vol=ncid.variables[e3tvar]
        except KeyError:
            pass
        else:
            break
    else:
        raise Exception('Could not find e3t field in mesh mask file.')
    e1t=ncid.variables['e1t']
    e2t=ncid.variables['e2t']
    if vol.ndim == 4:
        cube_index_with_time = [0]+cube_index[-3:]
        cube_index_no_depth = [0]+cube_index[-2:] 
        e1t=e1t[tuple(cube_index_no_depth)]
        e2t=e2t[tuple(cube_index_no_depth)]
        vol=vol[tuple(cube_index_with_time)]
    else:
        e1t=e1t[tuple(cube_index[-2:])]
        e2t=e2t[tuple(cube_index[-2:])]
        vol=vol[tuple(cube_index[-3:])]
        if dummy_cube.ndim == 4:
            vol=iris.util.broadcast_weights(vol,dummy_cube.data,[1,2,3])
    
    vol=vol[:]*e1t[:]*e2t[:]
    
    ncid.close()
    
    volmean_ts=[gt.read_cube(fname,fieldname)[tuple(cube_index)].collapsed(['depth','latitude','longitude'],
                               iris.analysis.MEAN,
                               weights=vol) 
                               for fname in files]
    
    volmean_ts=[cube.data for cube in volmean_ts]
     
    return volmean_ts

def plot_volmean_timeseries(infiles=None,fieldname=None,meshfile=None,mask=None,firstlevel=None,lastlevel=None,
                            north=None,south=None,east=None,west=None,outfile=None,title=None,
                            xmin=None,xmax=None,ymin=None,ymax=None,noshow=False):

    region = [north,east,south,west]
    volmean_ts = volmean_timeseries(files=infiles,fieldname=fieldname,meshfile=meshfile,region=region,
                                    firstlevel=firstlevel,lastlevel=lastlevel,mask=mask)   

    plt.plot(volmean_ts)
    if outfile is not None:
        matplotlib.rcParams['font.size'] =8
        plt.savefig(outfile,dpi=200)
    elif not noshow:
        plt.show()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", action="store",dest="infiles",nargs='+',default=None,
                    help="name(s) of input file(s).")
    parser.add_argument("-f", "--fieldname", action="store",dest="fieldname",default=None,
                    help="name of field - standard_name or variable name.")
    parser.add_argument("-m", "--meshfile", action="store",dest="meshfile",default=None,
                    help="name of meshfile.")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",default=None,
                    help="plot to output file - filetype by extension. If unset plot to GUI.")
    parser.add_argument("-l", "--firstlevel", action="store",dest="firstlevel",type=int,default=None,
                    help="first level in volume required (counting from zero)")
    parser.add_argument("-L", "--lastlevel", action="store",dest="lastlevel",type=int,default=None,
                    help="last level in volume required (counting from zero)")
    parser.add_argument("-x", "--xmin", action="store",dest="xmin",type=float,default=None,
                    help="start of x-axis")
    parser.add_argument("-X", "--xmax", action="store",dest="xmax",type=float,default=None,
                    help="end of x-axis")
    parser.add_argument("-y", "--ymin", action="store",dest="ymin",type=float,default=None,
                    help="start of y-axis")
    parser.add_argument("-Y", "--ymax", action="store",dest="ymax",type=float,default=None,
                    help="end of y-axis")
    parser.add_argument("-t", "--title", action="store",dest="title",default=None,
                    help="title for plot")
    parser.add_argument("-W", "--west", action="store",dest="west",type=float,default=None,
                    help="western limit of area to plot")
    parser.add_argument("-E", "--east", action="store",dest="east",type=float,default=None,
                    help="eastern limit of area to plot")
    parser.add_argument("-S", "--south", action="store",dest="south",type=float,default=None,
                    help="southern limit of area to plot")
    parser.add_argument("-N", "--north", action="store",dest="north",type=float,default=None,
                    help="northern limit of area to plot")
    parser.add_argument("-n", "--noshow", action="store_true",dest="noshow",default=False,
                    help="don't show GUI window even if outfile unset")

    args = parser.parse_args()

    plot_volmean_timeseries(infiles=args.infiles,fieldname=args.fieldname, meshfile=args.meshfile, 
                            north=args.north, south=args.south, west=args.west, east=args.east,
                            firstlevel=args.firstlevel,lastlevel=args.lastlevel, outfile=args.outfile, 
                            title=args.title, xmin=args.xmin, xmax=args.xmax, 
                            ymin=args.ymin, ymax=args.ymax, noshow=args.noshow)
