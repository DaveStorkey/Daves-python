#! /usr/bin/env python
'''
Module to do area or volume integration for NEMO fields

To do:
   1. Put in volume integration (and manage choice of area or volume integration).
   2. Auto-select meshfile depending on size of input array (with option to manually specify meshmask). 
   3. Remove hardwiring for grid_T (maybe move to using Iris cubes). 

@author: Dave Storkey
'''

import netCDF4 as nc
import numpy as np
import numpy.ma as ma

def geo_integrate(meshfile=None,datafile=None,fieldname=None,
                  area=None, volume=None,record=None,
                  east=None, west=None, south=None, north=None,
                  top=None, bottom=None):

    if record is None:
        record=0

    meshfile_in = nc.Dataset(meshfile,'r')

    e1 = meshfile_in.variables['e1t'][:]
    e2 = meshfile_in.variables['e2t'][:]
#    e3 = meshfile_in.variables['e3t'][:]
    tmaskutil = meshfile_in.variables['tmaskutil'][:]

    datafile_in = nc.Dataset(datafile,'r')

    for lonvar in 'nav_lon', 'nav_lon_grid_T':
        try:
            lons = datafile_in.variables[lonvar][:]
        except KeyError:
            pass
        else:
            break
    else:
        raise Exception('Could not find longitude field in file.')

    for latvar in 'nav_lat', 'nav_lat_grid_T':
        try:
            lats = datafile_in.variables[latvar][:]
        except KeyError:
            pass
        else:
            break
    else:
        raise Exception('Could not find latitude field in file.')

#    depths = datafile_in.variables['deptht'][:]
    select_mask = np.where(lons < -180.0,1,0)
    lons = lons + 360.0*select_mask
    select_mask = np.where(lons > 180.0,1,0)
    lons = lons - 360.0*select_mask

    field_in_file = datafile_in.variables[fieldname]
    # make a copy in case we have to modify the mask later
    try:
        field = ma.copy(datafile_in.variables[fieldname][record,:,:])
    except(ValueError):
        # no record dimension
        field = ma.copy(datafile_in.variables[fieldname][:,:])
      
    (area_integral,area_average) = area_integrate(field=field,e1=e1,e2=e2,mask=tmaskutil,
               lats=lats,lons=lons,east=east,west=west,south=south,north=north)

    if "units" in field_in_file.ncattrs():
        print "area_integral : ",area_integral,field_in_file.units,"*m2"
        print "area_average : ",area_average,field_in_file.units
    else:
        print "area_integral : ",area_integral," units?"
        print "area_average : ",area_average," units?"

    meshfile_in.close()
    datafile_in.close()        


def area_integrate(field=None,e1=None,e2=None,mask=None,lats=None,lons=None,
                   east=None, west=None, south=None, north=None, plainsum=None):

    if east is not None:
        if east < -180.0:
            east=east+360.0
        if east > 180.0:
            east=east-360.0
        field = ma.masked_where(lons > east, field)
        
    if west is not None:
        if west < -180.0:
            west=west+360.0
        if west > 180.0:
            west=west-360.0
        field = ma.masked_where(lons < west, field)

    if south is not None:
        field = ma.masked_where(lats < south, field)
        
    if north is not None:
        field = ma.masked_where(lats > north, field)
        
#    print 'ma.count(field) : ',ma.count(field)
#    print 'ma.count_masked(field) : ',ma.count_masked(field)

    if mask is None:
        mask=np.ones(e1.shape)

    if plainsum:
        area_integral = ma.sum(field[:]*mask[:])
        area_average = ma.sum(field[:]*mask[:])/ma.sum(e1[:]*e2[:])
    else:
        area_integral = ma.sum(field[:]*e1[:]*e2[:]*mask[:])
        area_average = ma.sum(field[:]*e1[:]*e2[:]*mask[:])/ma.sum(e1[:]*e2[:])

    return(area_integral, area_average)


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--meshfile", action="store",dest="meshfile",default=None,
                    help="meshmask file with e1 and e2 scale factors")
    parser.add_argument("-d", "--datafile", action="store",dest="datafile",default=None,
                    help="data file")
    parser.add_argument("-f", "--fieldname", action="store",dest="fieldname",default=None,
                    help="field name")
    parser.add_argument("-W", "--west", action="store",dest="west",type=float,default=None,
                    help="western limit of area to integrate (deg lon)")
    parser.add_argument("-E", "--east", action="store",dest="east",type=float,default=None,
                    help="eastern limit of area to integrate (deg lon)")
    parser.add_argument("-S", "--south", action="store",dest="south",type=float,default=None,
                    help="southern limit of area to integrate (deg lat)")
    parser.add_argument("-N", "--north", action="store",dest="north",type=float,default=None,
                    help="northern limit of area to integrate (deg lat)")
    parser.add_argument("-B", "--bottom", action="store",dest="bottom",type=float,default=None,
                    help="lower limit of volume to integrate (metres below sea surface)")
    parser.add_argument("-T", "--top", action="store",dest="top",type=float,default=None,
                    help="upper limit of volume to integrate (metres below sea surface)")
 
    args = parser.parse_args()

    geo_integrate(meshfile=args.meshfile, datafile=args.datafile, fieldname=args.fieldname,
                  west=args.west,east=args.east,south=args.south,north=args.north,
                  bottom=args.bottom,top=args.top)

