#! /usr/bin/env python
'''
Routine to integrate freshwater fluxes around Antarctica.

@author: Dave Storkey
'''

import netCDF4 as nc
import numpy.ma as ma
import geo_integrate as geo

def antarctic_fw_flux(filenames=None,meshfilename=None,fieldlist=None,record=None,
                  east=None, west=None, south=None, north=None):

    meshfile = nc.Dataset(meshfilename,'r')

    e1 = meshfile.variables['e1t'][:]
    e2 = meshfile.variables['e2t'][:]
    e3 = meshfile.variables['e3t'][:]
    lats = meshfile.variables['gphit'][:]
    lons = meshfile.variables['glamt'][:]

    if record is None:
        record=0

    infiles=[]
    for filename in filenames:
        infiles.append(nc.Dataset(filename,'r'))

    for fieldname in fieldlist:
        field_file_in = None
        for infile in infiles:
            try:
                field_file_in = infile.variables[fieldname]
            except KeyError:
                pass
            else:
                infile_to_use = infile
                break
        if field_file_in is None:        
            raise Exception("Error: Can't find "+fieldname+" in any input file")

        field_in = ma.copy(infile_to_use.variables[fieldname][record,:,:])

# Hack for diagnostic CICE fields which are missing the southernmost row that NEMO has:
        if field_in.shape[0] == lats.shape[0] - 1:
            field = ma.zeros(lats.shape)
            field[1:,1:-1] = field_in[:,:]
        else:
            field = field_in

        plainsum = None
        Gt_per_year_conversion_factor = 1.0e-12 * 3.1536e+07
        if fieldname == "fresh_ai":
            plainsum = False           
            # fresh_ai is in units of cm/day!    
            Gt_per_year_conversion_factor = 1.0e-02 * 1000.0 * 1.0e-12 * 365.0

        (field_ttl,field_avg) = geo.area_integrate(field=field,e1=e1,e2=e2,
               lats=lats,lons=lons,east=east,west=west,south=south,north=north,plainsum=plainsum) 

#        print "Gt_per_year_conversion_factor is ",Gt_per_year_conversion_factor
        field_ttl_Gt_per_yr = field_ttl * Gt_per_year_conversion_factor 
#        print fieldname," area integral : ",field_ttl,field_file_in.units,"*m2"
        print fieldname," area integral : ",field_ttl_Gt_per_yr," Gt/yr"


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--filenames", action="store", dest="filename", default=None, nargs="+",
                    help="name(s) of input file(s)")
    parser.add_argument("-m", "--meshfilename", action="store",dest="meshfilename",default=None,
                    help="meshmask file with e1 and e2 scale factors")
    parser.add_argument("-f", "--fieldlist", action="store",dest="fieldlist",default=None,nargs="+",
                    help="meshmask file with e1 and e2 scale factors")
    parser.add_argument("-r", "--record", action="store",dest="record",type=int,default=None,
                    help="record number to read in file (unlimited dimension)")
    parser.add_argument("-W", "--west", action="store",dest="west",type=float,default=None,
                    help="western limit of area to integrate (deg lon)")
    parser.add_argument("-E", "--east", action="store",dest="east",type=float,default=None,
                    help="eastern limit of area to integrate (deg lon)")
    parser.add_argument("-S", "--south", action="store",dest="south",type=float,default=None,
                    help="southern limit of area to integrate (deg lat)")
    parser.add_argument("-N", "--north", action="store",dest="north",type=float,default=None,
                    help="northern limit of area to integrate (deg lat)")

    args = parser.parse_args()

    antarctic_fw_flux(args.filename,meshfilename=args.meshfilename,fieldlist=args.fieldlist,
                  west=args.west,east=args.east,south=args.south,north=args.north)
