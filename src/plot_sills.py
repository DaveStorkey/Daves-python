#! /usr/bin/env python

'''
Wrapper script for plot_nemo.py to plot the bathymetry in the
region of the sills/channels identified in "Sills of the Global
Ocean" by Simon Thompson (1995).

Read in details of each location from a simple database.

@author: Dave Storkey
@date: Apr 2024
'''

import os
import distutils
import xarray as xr
import numpy as np
import plot_nemo as pn
from cutout_domain import cutout_domain

def plot_sills(database=None, filenames=None, vars=None, titles=None, cutout=False):

    if database is None:
        raise Exception("Error: must specify database file.")

    if filenames is None:
        raise Exception("Error: must specify at least one input file.")
    elif not isinstance(filenames,list):
        filenames=[filenames]

    if vars is None:
        vars=len(filenames)*["Bathymetry"]
    elif not isinstance(vars,list):
        vars=[vars]
    if len(vars) == 1 and len(filenames) > 1:
        vars=len(filenames)*vars

    if titles is None:
        titles=""
    if not isinstance(titles,list):
        titles=[titles]
    if len(titles) == 1 and len(filenames) > 1:
        titles=len(filenames)*titles

    if cutout is None:
        cutout=False
    if not isinstance(cutout,list):
        cutout=[cutout]
    if len(cutout) == 1 and len(filenames) > 1:
        cutout=len(filenames)*cutout

    # parse database
    s=[]
    try:
        with open(database) as db:
            for cline in db:
                if cline[0]=="#" or "|" not in cline:
                    # allow comments in database
                    continue
                att=cline.split('|')
                s.append({})
                s[-1]["name"] = att[0].strip()
                s[-1]["depth"] = float(att[1].strip())
                s[-1]["width"] = att[2].strip()
                s[-1]["lat"] = float(att[3].strip())
                s[-1]["lon"] = float(att[4].strip())
                s[-1]["ref"] = att[5].strip()
    except FileNotFoundError:
        raise Exception("Error: file "+database+" not found.")

    print(len(s)," entries found in database "+database)

    # loop over entries in database and plot 
    for sill in s:
        print("")
        print("Plotting "+sill["name"])
        print("")
        draw_points = [sill["lon"]-0.1,sill["lat"]-0.1,sill["lon"]+0.1,sill["lat"]+0.1,
                       sill["lon"]-0.1,sill["lat"]+0.1,sill["lon"]+0.1,sill["lat"]-0.1]
        depmax = sill["depth"]*2.0
        compressed_name = sill["name"].replace(" ","")

        for filename,var,title,precut in zip(filenames,vars,titles,cutout):
            filestem=filename.replace(".nc","")

            south = sill["lat"]-2.5
            north = sill["lat"]+2.5
            west = sill["lon"]-2.5
            east = sill["lon"]+2.5
            if filename == filenames[0]:
                print("south, north, west, east : ",south,north,west,east)

            print("")
            print("precut : ",precut)
            print("")
            if precut:
                cutout_domain(infile=filename,invar=var,outfile=filestem+"_"+compressed_name+".nc",
                             south=south,north=north,west=west,east=east)
                filename=filestem+"_"+compressed_name+".nc"
                south = None
                north = None
                west  = None
                east  = None

            (cslines, cscolor, csarrows) = pn.plot_nemo(filenames=filename,sca_names=var,
                                           plot_types="b",cmap="tab20b_r",
                                           mnfld=0.0,mxfld=depmax,west=west,east=east, 
                                           south=south,north=north,vertbar=True,nlevs=21,
                                           outfile=filestem+"_"+compressed_name+".png",
                                           title=title+": "+sill["name"],facecolor="white",
                                           draw_points=draw_points)

        if precut:
            try:
                os.remove(filename)
            except FileNotFoundError:
                pass

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--database", action="store",dest="database",
        help="name of database file with details of locations to plot.")
    parser.add_argument("-i", "--filenames", action="store",dest="filenames",nargs="+",
        help="name(s) of database file with details of locations to plot.")
    parser.add_argument("-v", "--vars", action="store",dest="vars",nargs="+",
        help="name(s) of variables to plot")
    parser.add_argument("-t", "--titles", action="store",dest="titles",nargs="+",
        help="descriptions of the data from each input file eg. GEBCO2021 or ORCA025")
    parser.add_argument("-C", "--cutout", action="store",dest="cutout",nargs="+",
        type=lambda x:bool(distutils.util.strtobool(x)),
        help="if True cutout the data prior to plotting (for large data sets)")
    args = parser.parse_args()
 
    plot_sills(database=args.database,filenames=args.filenames,vars=args.vars,
               titles=args.titles,cutout=args.cutout)
