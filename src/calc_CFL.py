#! /usr/bin/env python
'''
Routine to calculate CFL parameters from 3D velocity field components and a coordinates file. 

Nov 2021 : umask, vmask options. DS.

@author: Dave Storkey
'''

import netCDF4 as nc
import numpy as np

def calc_CFL(timestep=None,filename_coordinates=None,filename_velocities=None,filename_out=None):

    file_coordinates = nc.Dataset(filename_coordinates,'r')
    e1u = file_coordinates.variables['e1u']
    e2v = file_coordinates.variables['e2v']

    file_velocities = nc.Dataset(filename_velocities,'r')
    indim = file_velocities.dimensions
    uu=None
    vv=None
    for uvar in 'un','uo','vozocrtx':
        try:
            uu = file_velocities.variables[uvar]
        except KeyError:
            pass
        else:
            break 
    for vvar in 'vn','vo','vomecrty':
        try:
            vv = file_velocities.variables[vvar]
        except KeyError:
            pass
        else:
            break 

    if not uu and not vv:
        raise Exception("Can't find U or V velocity components in file "+filename_velocities)

    file_out = nc.Dataset(filename_out,'w')
    for dim in indim:
        file_out.createDimension(dim, len(indim[dim]))

    if uu:
        if "_FillValue" in uu.ncattrs():
            file_out.createVariable('cfl_u',datatype='f',dimensions=uu.dimensions,fill_value=uu._FillValue)
        else:
            file_out.createVariable('cfl_u',datatype='f',dimensions=uu.dimensions)
        file_out.variables['cfl_u'][:] = 0.0
        # factor of 2 to account for leapfrog timestep
        file_out.variables['cfl_u'][:,:,:] = abs(2*uu[:,:,:]*timestep/e1u[:,:])

    if vv:
        if "_FillValue" in vv.ncattrs():
            file_out.createVariable('cfl_v',datatype='f',dimensions=vv.dimensions,fill_value=vv._FillValue)
        else:
            file_out.createVariable('cfl_v',datatype='f',dimensions=vv.dimensions)
        file_out.variables['cfl_v'][:] = 0.0
        # factor of 2 to account for leapfrog timestep
        file_out.variables['cfl_v'][:,:,:] = abs(2*vv[:,:,:]*timestep/e2v[:,:])
        
    file_coordinates.close()
    file_velocities.close()
    file_out.close()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--timestep", action="store",dest="timestep",type=int,default=None,
                    help="model timestep in seconds")
    parser.add_argument("-c", "--coordinates", action="store",dest="filename_coordinates",default=None,
                    help="name of coordinates file")
    parser.add_argument("-v", "--velocities", action="store",dest="filename_velocities",default=None,
                    help="name of file with bottom fields")
    parser.add_argument("-o", "--out", action="store",dest="filename_out",default=None,
                    help="name of file to write CFL fields to")
 
    args = parser.parse_args()

    calc_CFL(timestep=args.timestep, filename_coordinates=args.filename_coordinates,
                 filename_velocities=args.filename_velocities,
                 filename_out=args.filename_out)

