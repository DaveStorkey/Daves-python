#! /usr/bin/env python
'''
Module to calculate a timeseries of
domain-integrated kinetic energy from
model U and V fields. Option to calculate
for a subdomain.

NB. Assumes uniform horizontal areas and
vertical thicknesses if it can't find the
relevant variables in the file.

Note that for now it doesn't do much to
check that the timeseries of U field and V
fields that you give it are valid for the
same set of times. Just checks that you
have the same number of time levels for each.

@author: Dave Storkey
'''

import iris
import numpy.ma as ma
import general_tools as gt
import area_integral_iris as ai

def ke_integrate(ufiles=None,ufieldname=None,vfiles=None,vfieldname=None,
                  east=None, west=None, south=None, north=None,
                  top=None, bottom=None):

    ke_timeseries_u=[]
    ke_timeseries_v=[]
    
    for ufile in ufiles:
        uvel = gt.read_cube(ufile,ufieldname)
        try:
            cell_area = gt.read_cube(ufile,"cell_area")
        except:
            cell_area = None
        try:
            cell_thickness = gt.read_cube(ufile,"cell_thickness")
        except:
            cell_thickness = None

    

    return(field_area_integral,field_area_average)

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-u", "--ufiles", action="store",dest="ufiles",default=None,nargs="+",
                    help="data file for zonal velocity")
    parser.add_argument("-U", "--ufieldname", action="store",dest="ufieldname",default=None,
                    help="field name for zonal velocity")
    parser.add_argument("-v", "--vfiles", action="store",dest="vfiles",default=None,nargs="+",
                    help="data file for meridional velocity")
    parser.add_argument("-V", "--vfieldname", action="store",dest="vfieldname",default=None,
                    help="field name for meridional velocity")
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

    ke_integrate(ufiles=args.ufiles, ufieldname=args.ufieldname, vfiles=args.vfiles, vfieldname=args.vfieldname, 
                  west=args.west,east=args.east,south=args.south,north=args.north,bottom=args.bottom,top=args.top)

