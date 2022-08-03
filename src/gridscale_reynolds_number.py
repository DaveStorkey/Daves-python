#! /usr/bin/env python

'''
Routine to calculate gridscale Reynolds number given an input
horizontal velocity field, field of horizontal grid spacings and
a field of biharmonic viscosity. 

Ref: Griffies and Hallberg (2000), Mon. Weath. Review vol 128 pp 2935-2946

@author: Dave Storkey
@date: August 2018
'''

import numpy as np
import numpy.ma as ma
import netCDF4 as nc

def gridscale_reynolds_number(infiles=None, fieldnames=None, meshfile=None, grid=None, outfile=None, viscosity=None):

    # infiles: one or two input files.
    # infields: the horizontal velocity field and biharmonic viscosity field in that order.

    if type(infiles) is not list:
        infiles = [infiles]
    if len(infiles) == 1:
        velfile = infiles[0]
        viscfile = infiles[0]
    elif len(infiles) == 2:
        velfile = infiles[0]
        viscfile = infiles[1]
    else:
        raise Exception("Error: must specify one or two input files.")

    if type(fieldnames) is not list:
        fieldnames = [fieldnames]
    if len(fieldnames) == 1 and viscosity is not None:
        velname = fieldnames[0]
    elif len(fieldnames) == 2:
        velname = fieldnames[0]
        viscname = fieldnames[1]
    else:
        raise Exception("Error: must specify name of velocity field in fieldnames and name of viscosity field if -v option not set.")

    if grid is None:
        if velname in ["vozocrtx","uo","un","ub"]:
            grid = "U"
        elif velname in ["vomecrty","vo","vn","vb"] :
            grid = "V"
        else:
            raise Exception("Error. Grid not specified and could not guess from velocity field name.")
    else:
        grid = grid.upper()
 
    with nc.Dataset(velfile,'r') as data_in:
        model_dims = data_in.dimensions
        nx = len(model_dims['x'])
        ny = len(model_dims['y'])
        for depthname in ['z', 'nav_lev', 'depthu', 'depthv']:
            try:
                nz = len(model_dims[depthname])
            except KeyError:
                pass
            else:
                break
        else:
            raise Exception("Error: could not find depth coordinate (z, nav_lev, depthu or depthv) in velocity input file.")   
        vel =  data_in.variables[velname][:]

    with nc.Dataset(meshfile,'r') as data_in:
        model_dims = data_in.dimensions
        nx2 = len(model_dims['x'])
        ny2 = len(model_dims['y'])
        nz2 = len(model_dims['z'])
        nav_lon = data_in.variables['nav_lon'][:]
        nav_lat = data_in.variables['nav_lat'][:]
        depth = data_in.variables['gdept_0'][:]
        if grid == "U":
            dx =  data_in.variables['e1u'][:]
        elif grid == "V":
            dx =  data_in.variables['e2v'][:]
        else:
            raise Exception("Error. Unrecognised grid.")
    # sanity check
    if nx2 != nx or ny2 != ny or nz2 != nz:
        raise Exception("Error. Inconsistent dimension lengths between velocity and mesh mask input files.")    

    if viscosity is None:
        with nc.Dataset(viscfile,'r') as data_in:
            model_dims = data_in.dimensions
            nx1 = len(model_dims['x'])
            ny1 = len(model_dims['y'])
            nz1 = len(model_dims['z'])
            visc =  data_in.variables[viscname][:]
        # sanity check:
        if nx1 != nx or ny1 != ny or nz1 != nz:
            raise Exception("Error. Inconsistent dimension lengths between velocity and viscosity input files.")    
    else:
        visc = ma.zeros((nz,ny,nx))
        visc[:,:,:] = viscosity * np.cos(nav_lat[:,:]*np.pi/180.0)

    reynolds = np.abs(vel) * dx * dx * dx / np.abs(visc)

    with nc.Dataset(outfile,'w') as data_out:
        data_out.createDimension('x', nx)
        data_out.createDimension('y', ny)
        data_out.createDimension('z', nz)
        data_out.createVariable('nav_lon',datatype='f',dimensions=('y','x'))
        data_out.variables['nav_lon'][:] = nav_lon[:]
        data_out.variables['nav_lon'].standard_name = 'longitude'
        data_out.variables['nav_lon'].units = 'degrees_east'
        data_out.createVariable('nav_lat',datatype='f',dimensions=('y','x'))
        data_out.variables['nav_lat'][:] = nav_lat[:]
        data_out.variables['nav_lat'].standard_name = 'latitude'
        data_out.variables['nav_lat'].units = 'degrees_north'
        data_out.createVariable('depth',datatype='f',dimensions=('z'))
        data_out.variables['depth'][:] = depth[:]
        data_out.variables['depth'].standard_name = 'depth'
        data_out.variables['depth'].units = 'm'
        data_out.createVariable('gridRe',datatype='f',dimensions=('z','y','x'))
        data_out.variables['gridRe'][:] = reynolds[:]
        data_out.variables['gridRe'].long_name = 'gridscale Reynolds number'
        data_out.variables['gridRe'].coordinates = 'nav_lat nav_lon'


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", action="store",dest="infiles",nargs="+",
                         help="name(s) of one or two input files")
    parser.add_argument("-f", "--fieldnames", action="store",dest="fieldnames",nargs="+",
                         help="names of input fields: horizontal velocity and biharmonic viscosity")
    parser.add_argument("-m", "--mesh", action="store",dest="meshfile",
                         help="name of mesh mask file")
    parser.add_argument("-g", "--grid", action="store",dest="grid",
                         help="grid type of velocity field (U or V). Will try to guess if not supplied.")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",
                         help="name of output file")
    parser.add_argument("-v", "--viscosity", action="store",dest="viscosity",type=float,
                         help="scalar biharmonic viscosity (to be scaled with latitude)")

    args = parser.parse_args()

    gridscale_reynolds_number(infiles=args.infiles, fieldnames=args.fieldnames, meshfile=args.meshfile, 
                              grid=args.grid, outfile=args.outfile, viscosity=args.viscosity)
