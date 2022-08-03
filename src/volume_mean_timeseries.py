#! /usr/bin/env python
'''
Module to calculate volume mean timeseries from 
a timeseries of 3D fields. 

NB. If coordinates file (or mesh mask file) is not specified then
it assumes uniform grid spacing in the horizontal
and vertical

@author: Dave Storkey
'''

import iris
import netCDF4 as nc

def volume_mean_timeseries(infile=None,coordsfile=None,outfile=None,fieldname=None):

    incube = iris.load_cube(infile,fieldname)
    if coordsfile is None:
        print 'WARNING: no coordinates file specified. Assuming uniform cell volumes.'
        incube_volmean = incube.collapsed(('longitude','latitude','Vertical T levels'),iris.analysis.MEAN)
    else:
        with nc.Dataset(coordsfile,'r') as coords:
            e1t = coords.variables['e1t'][:] 
            e2t = coords.variables['e2t'][:] 
        with nc.Dataset(infile,'r') as indata:
            e3t = indata.variables['e3t'][:]
        volume = e1t*e2t*e3t
        incube_volmean = incube.collapsed(('longitude','latitude','Vertical T levels'),iris.analysis.MEAN,
                                           weights=volume)

    iris.save(incube_volmean,outfile)

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", action="store",dest="infile",default=None,
                    help="input filename")
    parser.add_argument("-c", "--coordsfile", action="store",dest="coordsfile",default=None,
                    help="input filename")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",default=None,
                    help="output filename")
    parser.add_argument("-f", "--fieldname", action="store",dest="fieldname",default=None,
                    help="field name")
 
    args = parser.parse_args()

    volume_mean_timeseries(infile=args.infile, coordsfile=args.coordsfile, outfile=args.outfile, fieldname=args.fieldname)


