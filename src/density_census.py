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

def density_census(infile=None,dens_name=None,threshold=None,box=None,outfile=None):

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
#        dens.data   = ma.masked_where(dens.data[:] < threshold, dens.data[:])
        thick.data  = ma.masked_where(dens.data[:] < threshold, thick.data[:])
        volume.data = ma.masked_where(dens.data[:] < threshold, volume.data[:])
        if box is not None:
            # Note haven't been very clever with longitudes here...
            ones = np.ones(volume.shape)
            lats = ones * volume.coord("latitude").points
            lons = ones * volume.coord("longitude").points
            print("box : ",box)
            if box[0] != "None":
                # west limit
                thick.data = ma.masked_where(lons < float(box[0]), thick.data[:]) 
                volume.data = ma.masked_where(lons < float(box[0]), volume.data[:]) 
            if box[1] != "None":
                # east limit
                thick.data = ma.masked_where(lons > float(box[1]), thick.data[:]) 
                volume.data = ma.masked_where(lons > float(box[1]), volume.data[:]) 
            if box[2] != "None":
                # south limit
                thick.data = ma.masked_where(lats < float(box[2]), thick.data[:]) 
                volume.data = ma.masked_where(lats < float(box[2]), volume.data[:]) 
            if box[3] != "None":
                # north limit
                thick.data = ma.masked_where(lats > float(box[3]), thick.data[:]) 
                volume.data = ma.masked_where(lats > float(box[3]), volume.data[:]) 
        thickbottom = thick.collapsed("depth",iris.analysis.SUM)
        thickbottom.var_name = "thickAABW"
        thickbottom.long_name = "thickness of AABW"
        volsum = volume.collapsed(["depth","latitude","longitude"],iris.analysis.SUM)
        # convert to km^3
        volsum.data[:] = volsum.data[:] * 1.0e-9
        volsum.var_name = "volsum"

        # to get zonal mean first slice the "thick" cube to get a cube with the right
        # dimensions, then do a numpy weighted sum over the i-dimension of the original
        # "thick" field and overwrite the data of the sliced cube. For zonal sums in the
        # southern hemisphere we can just sum along x grid lines. 
        volslice = thickbottom[:,:,0]
        volslice.var_name = "volAABW_zonalint"
        volslice.data[:] = ma.sum(thickbottom.data[:]*area[:],axis=-1)
        delta_lat = volslice.coord('latitude').bounds.max(axis=-1) - \
                    volslice.coord('latitude').bounds.min(axis=-1)
        volslice.data[:] = volslice.data[:]/delta_lat[:]

    if outfile is not None:
        iris.save([volume,dens,thickbottom,volslice,volsum],outfile)

#    print(volume)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file",action="store",dest="infile",
        help="name of input file")
    parser.add_argument("-d", "--dens_name",action="store",dest="dens_name",
        help="name of density variable")
    parser.add_argument("-B", "--box",action="store",dest="box",nargs=4,
        help="area to integrate volume over: W, E, S, N. Enter None for no limit.")
    parser.add_argument("-T", "--threshold",action="store",dest="threshold",type=float,
        help="density threshold")
    parser.add_argument("-o", "--output_file",action="store",dest="outfile",
        help="name of output file")

    args = parser.parse_args()

    density_census(infile=args.infile,dens_name=args.dens_name,threshold=args.threshold,box=args.box,
                   outfile=args.outfile)

