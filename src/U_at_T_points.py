#! /usr/bin/env python

'''
Routine to . 

@author: Dave Storkey
@date: August 2018
'''

import iris
import numpy as np
import general_tools as gt

def rotate(ucube,vcube):

    if ucube.shape[-1] == 362 and ucube.shape[-2] == 332:
        model='eORCA1'
    elif ucube.shape[-1] == 1442 and ucube.shape[-2] == 1207:
        model='eORCA025'
    elif ucube.shape[-1] == 4322 and ucube.shape[-2] == 3606:
        model='eORCA12'
    else:
        raise Exception("Error: I can't guess which model it is to get the angles for vector rotation")

    sin_t = gt.read_cube('/project/nemo/interp/nemo_remap/angle_'+model+'.nc','gsint')
    cos_t = gt.read_cube('/project/nemo/interp/nemo_remap/angle_'+model+'.nc','gcost')

#    if ucube.shape != cos_t.shape:
#        sin_t = sin_t[1:-1,1:-1]
#        cos_t = cos_t[1:-1,1:-1]

    ucube.data = (ucube.data * cos_t.data) - (vcube.data * sin_t.data)
    vcube.data = (vcube.data * cos_t.data) + (ucube.data * sin_t.data)
    
    return ucube, vcube

def U_at_T_points(ucube,vcube):

    ulons = ucube.coord('longitude').points
    vlats = vcube.coord('latitude').points

    # derive a mask on the T points
    mask_U = ucube.data.mask.copy()
    mask_U_m1 = np.roll(mask_U,1,axis=-1)
    mask_V = vcube.data.mask.copy()
    mask_V_m1 = np.roll(mask_V,1,axis=-2)
    mask_T = ( mask_U & mask_U_m1 ) & ( mask_V & mask_V_m1 )

    # this unmasks all points...
    ucube.data[ucube.data.mask] = 0.0
    vcube.data[vcube.data.mask] = 0.0

    # Average U field (and longitudes) in the x-direction.
    u_m1 = np.roll(ucube.data, 1,axis=-1)
    ulons_m1 = np.roll(ulons, 1,axis=-1)
    ucubeT = ucube.copy()
    ucubeT.data = 0.5 * ( ucube.data + u_m1 )
    ucubeT.data.mask = mask_T

    # Average V field (and latitudes) in the y-direction.
    v_m1 = np.roll(vcube.data, 1,axis=-2)
    vlats_m1 = np.roll(vlats, 1,axis=-2)
    vcubeT = vcube.copy()
    vcubeT.data = 0.5 * ( vcube.data + v_m1 )
    vcubeT.data.mask = mask_T

    # Set approximate lat/lon coordinates.
    ucubeT.coord('longitude').points = 0.5 * ( ulons + ulons_m1 )
    ucubeT.coord('latitude').points = 0.5 * ( vlats + vlats_m1 )
    vcubeT.coord('longitude').points =  ucubeT.coord('longitude').points
    vcubeT.coord('latitude').points =  ucubeT.coord('latitude').points

    return ucubeT,vcubeT

def U_at_T_points_wrapper(filenamestem=None,ufield=None,vfield=None,outfile=None):

    if ufield is None:
        ufield = 'sea_water_x_velocity'
    if vfield is None:
        vfield = 'sea_water_y_velocity'

    ucube = iris.load_cube(filenamestem+'_grid-U.nc',ufield)
    vcube = iris.load_cube(filenamestem+'_grid-V.nc',vfield)   
    
    ucubeT,vcubeT = U_at_T_points(ucube,vcube)

    iris.FUTURE.netcdf_no_unlimited = True
    iris.save([ucubeT,vcubeT],outfile)

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("filenamestem", help="Input filename stem (for grid-U and grid-V")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",default=None,
                   help="output file")
    parser.add_argument("-u", "--ufield", action="store",dest="ufield",default=None,
                   help="standard name of u-component of field")
    parser.add_argument("-v", "--vfield", action="store",dest="vfield",default=None,
                   help="standard name of v-component of field")

    args = parser.parse_args()

    U_at_T_points_wrapper(args.filenamestem,ufield=args.ufield,vfield=args.vfield,outfile=args.outfile)
