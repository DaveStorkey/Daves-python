#! /usr/bin/env python

'''
Utility routine to read and extract data from NEMO icebergs trajectory files.

Currently just extracts the iceberg counts, globally and for north and south
hemispheres.  

Created Nov 2016.

@author: Dave Storkey
'''

import netCDF4 as nc
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import plot_numbers as pn

def icebergs(infiles=None,mean=None):

    ntraj_global_read=[]
    ntraj_north_read=[]
    ntraj_south_read=[]
    for filename in infiles:
        file = nc.Dataset(filename,"r")
        ntraj_global_read.append(len(file.dimensions["ntraj"]))
        # read initial latitude for each iceberg
        lat = file.variables["lat"][:,0]
        # count number of icebergs in each hemisphere
        # assume that icebergs don't move between hemispheres :)
        ntraj_north_read.append(len(np.extract(lat>0,lat)))
        ntraj_south_read.append(len(np.extract(lat<0,lat)))

    filestem = infiles[0].split("_")[0]
    for (ntraj_read,label) in zip([ntraj_global_read,ntraj_north_read,ntraj_south_read],["global","north","south"]):

        if mean is not None:
            # lose any end values so len(array) is divided by mean exactly:
            nmeans = len(ntraj_read)/mean
            lentrim = nmeans*mean       
            ntraj_reshape = np.array(ntraj_read[0:lentrim]).reshape(nmeans,mean) 
            print "ntraj_reshape.shape: ",ntraj_reshape.shape
            ntraj = np.average(ntraj_reshape,axis=-1)
        else:
            ntraj = np.array(ntraj_read)

        with open(filestem+"_iceberg_count_"+label+".dat","w") as f:
            for ntraj_val in ntraj:
                f.write(str(ntraj_val)+"\n")
        
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", action="store", dest="infiles", nargs="+", default=None,
                    help="name of input file")
    parser.add_argument("-m", "--mean", action="store",dest="mean",type=int,default=None,
                    help="mean data in groups of this number")

    args = parser.parse_args()

    icebergs(infiles=args.infiles,mean=args.mean)
