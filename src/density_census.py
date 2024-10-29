#! /usr/bin/env python
'''
Routine to perform a density census given an 
input file with 3D density, area and thickness
variables.  

@author: Dave Storkey
@date: October 2024
'''

import argparse
import iris
import numpy as np
import numpy.ma as ma
import general_tools as gt

def density_census(infile=None,dens_name=None,threshold=None,outfile=None):

    if dens_name is None:
        dens_name="vosigmainsitu"

    thick = gt.read_cube(infile,"thkcello")
    dens = gt.read_cube(infile,dens_name)
    for cell_measure in thick.cell_measures():
        if cell_measure.measure == "area":
            area = cell_measure.data[:]
            break
        else:
            pass
    else:
        raise Exception("Could not find area cell_measure in thickness field")

    volume = thick.copy()
    volume.var_name = "volume"
    # cell_volume isn't a CF standard name (ho hum) so just unset it...
    volume.standard_name = None
    volume.long_name = "cell volume"
    volume.units = "m3"
    volume.data[:] = volume.data[:] * area[:]

    print("threshold : ",threshold)
    print("min/max densities : ",dens.data.min(),dens.data.max())
    if threshold is not None:
        thick.data = ma.masked_where(dens.data[:] < threshold, thick.data[:])
        volume.data = ma.masked_where(dens.data[:] < threshold, volume.data[:])
        dens.data = ma.masked_where(dens.data[:] < threshold, dens.data[:])
        depth = thick.collapsed("depth",iris.analysis.SUM)
        depth.var_name = "sigdepth"
        volsum = volume.collapsed(["depth","latitude","longitude"],iris.analysis.SUM)
        volsum.var_name = "volsum"

    if outfile is not None:
        iris.save([volume,dens,depth,volsum],outfile)

#    print(volume)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file",action="store",dest="infile",
        help="name of input file")
    parser.add_argument("-d", "--dens_name",action="store",dest="dens_name",
        help="name of density variable")
    parser.add_argument("-T", "--threshold",action="store",dest="threshold",type=float,
        help="density threshold")
    parser.add_argument("-o", "--output_file",action="store",dest="outfile",
        help="name of output file")

    args = parser.parse_args()

    density_census(infile=args.infile,dens_name=args.dens_name,threshold=args.threshold,outfile=args.outfile)

