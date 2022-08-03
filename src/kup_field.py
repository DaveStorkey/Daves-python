#! /usr/bin/env python
'''
Take a 3D NEMO field and turn it into a "k_up" field, ie. indexed in the vertical
by model level index above the topography, as in Le Sommer et al, Ocean Modelling (2009).

Use the bottom_level field from domain_cfg.nc or calculated from mask in input field. 
NB. If from domain_cfg.nc need to convert bottom_level indexing from fortran to C/python counting. 

To do : Include a "kup" coordinate variable??

Created Feb 2019

Oct 2020 : Python 3. DS.

@author: Dave Storkey
'''

import netCDF4 as nc
import numpy as np
import bottom_field as bf


def kup_field(infield=None,bottom_level=None):

    # NB. Assumes bottom_level counting is C/python, ie. top level = level 0, not 1.

    in_ones = np.ones(infield.shape)

    if len(infield.shape) == 4:
        tt = (np.array(range(infield.shape[0]) * in_ones.transpose())).transpose().astype(int)
        jj = (np.array(range(infield.shape[-2])) * in_ones.transpose(0,1,3,2)).transpose(0,1,3,2).astype(int)
        ii = (np.array(range(infield.shape[-1])) * in_ones).astype(int)
        bottom_level_broadcast = np.array([np.maximum(0,bottom_level-i) for i in range(infield.shape[-3])]).transpose(1,0,2,3).astype(int)
        print("tt.shape : ",tt.shape)
    elif len(infield.shape) == 3:
        jj = (np.array(range(infield.shape[-2])) * in_ones.transpose(0,2,1)).transpose(0,2,1).astype(int)
        ii = (np.array(range(infield.shape[-1])) * in_ones).astype(int)
        bottom_level_broadcast = np.array([np.maximum(0,bottom_level-i) for i in range(infield.shape[-3])]).astype(int)
    else:
        raise Exception("Error : can only handle 3 or 4 dimensions.")


    print("infield.shape : ",infield.shape)
    print("bottom_level_broadcast.shape : ",bottom_level_broadcast.shape)
    print("jj.shape : ",jj.shape)
    print("ii.shape : ",ii.shape)

    if len(infield.shape) == 4:
        outfield = infield[tt,bottom_level_broadcast,jj,ii]
    elif len(infield.shape) == 3:
        outfield = infield[bottom_level_broadcast,jj,ii]

    print("outfield.shape : ",outfield.shape)
    return outfield

def kup_field_wrapper(infile=None,fieldnames=None,outfile=None):

    data_in = nc.Dataset(infile,'r')
    if outfile is None:
        outfile = infile.replace(".nc","_kup.nc")
    data_out = nc.Dataset(outfile,'w')

    if fieldnames is None:
        fieldnames = data_in.variables.keys()

    for fieldname in fieldnames:
        print('Working on field : '+fieldname)
        try:
            field = data_in.variables[fieldname]
        except KeyError:
            raise Exception("Error: could not find "+fieldname+" in "+infile)

        found_depth = False
        outdims = []
        for dim in field.dimensions:
            if dim in ['z','depth','deptht','depthu','depthv','depthw']:
                # vertical dimension of output field is k_up
                found_depth = True
                if 'k_up' not in data_out.dimensions.keys():
                    data_out.createDimension('k_up',len(data_in.dimensions[dim]))
                outdims.append('k_up')
            else:
                if dim not in data_out.dimensions.keys():
                    data_out.createDimension(dim,len(data_in.dimensions[dim]))
                outdims.append(dim)
        if not found_depth or len(field.dimensions) < 3:
            print('Variable '+fieldname+' has no depth dimension. Skipping.')
            continue
        else:
            data_out.createVariable(fieldname,datatype='f',dimensions=outdims,fill_value=field._FillValue)

        # calculate bottom level index from mask in input field. 
        bottom_level, bottom_field = bf.bottom_field(infield=field[:])
 
        outfield = data_out.variables[fieldname]  
        print("outfield[:].shape : ",outfield[:].shape)
        outfield[:] = kup_field(infield=field[:],bottom_level=bottom_level)

        # Copy over the attributes adjusting the coordinates attribute as appropriate.
        attrdict={}
        for attr in field.ncattrs():
            attrdict[attr.encode('ascii')] = field.getncattr(attr)
        # Don't include the _FillValue attribute in the dictionary as this requires special treatment in createVariable.
        del attrdict['_FillValue']

        outfield.setncatts(attrdict)
        for depthvar in ['deptht','depthu','depthv','depthw','depth']:
            if depthvar in outfield.coordinates:
                outfield.coordinates = outfield.coordinates.replace(depthvar,'k_up')
                break
        else:
            raise Exception("Error: input field doesn't have a depth coordinate specified.")

    # Copy in all the coordinate and bounds variables and their attributes

    for var in 'nav_lat','nav_lon','area':
        invar=data_in.variables[var]
        data_out.createVariable(var,datatype='f',dimensions=invar.dimensions)
        data_out.variables[var][:] = invar[:]
        attrdict={}
        for attr in invar.ncattrs():
            attrdict[attr.encode('ascii')] = invar.getncattr(attr)
        data_out.variables[var].setncatts(attrdict)

    data_out.createVariable('k_up',datatype='f',dimensions=['k_up'])
    data_out.variables['k_up'].long_name="Model level index above bathymetry"    
    data_out.variables['k_up'].units="1"
    data_out.variables['k_up'].axis = "Z" 
    data_out.variables['k_up'].positive = "up" 
    data_out.variables['k_up'][:] = np.arange(len(data_out.variables['k_up']))+1

    data_in.close()
    data_out.close()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", action="store",dest="infile",
                    help="input file")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",
                    help="output file")
    parser.add_argument("-f", "--fieldnames", action="store",dest="fieldnames",nargs="+",
                    help="name of 3D field")

    args = parser.parse_args()

    kup_field_wrapper(infile=args.infile,fieldnames=args.fieldnames,outfile=args.outfile)
