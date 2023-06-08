#! /usr/bin/env python
'''
Extract bottom field and write to a separate file. 
NB. Requires input file to be masked. 

Created June 2016

Feb 2019 : Put the actual calculation in an inner routine
           which can be called from other python modules.
           Add option to output bottom level index as well
           as bottom field. DS. 

Oct 2020 : Bug fix typo which was giving the wrong bottom level
           in some cases. DS. 

June 2023 : Allow option to output bottom field projected onto masked
            3D grid. (To pass to cdftransport to calculate the bottom
            field component of a transport through a section). DS.

@author: Dave Storkey
'''

import netCDF4 as nc
import numpy as np
import numpy.ma as ma

def bottom_field(infield=None,wgrid=False):

    if len(infield.shape) == 4:
        bottom_level = ma.zeros([infield.shape[i] for i in (0,2,3)])
        bottom_field = ma.zeros([infield.shape[i] for i in (0,2,3)])
        if wgrid:
            # Special treatment required for W-grid fields, because field is not masked at the top of the bathymetry
            # and at least for W itself takes zero values here which aren't very interesting. We want to see the values
            # just above the bathymetry.
            for lev in [z+1 for z in range(len(infield[0,:,0,0])-2)]:
                bottom_level[:,:,:] = np.where( (~infield[:,lev,:,:].mask & infield[:,lev+1,:,:].mask), 
                                                  lev-1, bottom_level )
                bottom_field[:,:,:] = np.where( (~infield[:,lev,:,:].mask & infield[:,lev+1,:,:].mask), 
                                                  infield[:,lev-1,:,:], bottom_field )
        else:
            for lev in range(len(infield[0,:,0,0])-1):
                bottom_level[:,:,:] = np.where( (~infield[:,lev,:,:].mask & infield[:,lev+1,:,:].mask), 
                                                  lev, bottom_level )
                bottom_field[:,:,:] = np.where( (~infield[:,lev,:,:].mask & infield[:,lev+1,:,:].mask), 
                                                 infield[:,lev,:,:], bottom_field )
        
    elif len(infield.shape) == 3:
        bottom_level = ma.zeros([infield.shape[i] for i in (1,2)])
        bottom_field = ma.zeros([infield.shape[i] for i in (1,2)])
        if wgrid:
            # Special treatment required for W-grid fields, because field is not masked at the top of the bathymetry
            # and at least for W itself takes zero values here which aren't very interesting. We want to see the values
            # just above the bathymetry.
            for lev in [z+1 for z in range(len(infield[:,0,0])-2)]:
                bottom_level[:,:] = np.where( (~infield[lev,:,:].mask & infield[lev+1,:,:].mask), 
                                                  lev-1, bottom_level )
                bottom_field[:,:] = np.where( (~infield[lev,:,:].mask & infield[lev+1,:,:].mask), 
                                                 infield[lev-1,:,:], bottom_field )
        else:
            for lev in range(len(infield[:,0,0])-1):
                bottom_level[:,:] = np.where( (~infield[lev,:,:].mask & infield[lev+1,:,:].mask), 
                                                  lev, bottom_level )
                bottom_field[:,:] = np.where( (~infield[lev,:,:].mask & infield[lev+1,:,:].mask), 
                                                 infield[lev,:,:], bottom_field )
    else:
        raise Exception("ERROR : can only handle 4D or 3D fields!")

    return bottom_level,bottom_field


def bottom_field_wrapper(filename,fieldname,wgrid=False,level=False,out3D=False):

    infile = nc.Dataset(filename,'r')
    field = infile.variables[fieldname]
    if out3D:
        masked_ones = ma.ones(field.shape)
        masked_ones.mask = field[:].mask
    indim = field.dimensions
    indimlen = { dim : len(infile.dimensions[dim]) for dim in indim }
    if out3D:
        outfile = nc.Dataset(filename.replace(".nc","_bottom3D.nc"),'w')
    else:
        outfile = nc.Dataset(filename.replace(".nc","_bottom.nc"),'w')

    if fieldname in ['vovecrtz']:
        wgrid=True

    if len(field.dimensions) == 4:
        if out3D:
            for dim in indim:
                outfile.createDimension(dim,indimlen[dim])
            outfile.createVariable(fieldname,datatype='f',dimensions=list(outfile.dimensions.keys()),fill_value=field._FillValue)
            bot_field3D = outfile.variables[fieldname]
            bot_field2D = ma.zeros(tuple(indimlen[dim] for dim in [indim[0],indim[2],indim[3]]))
        else:
            for dim in indim[0],indim[2],indim[3]:
                outfile.createDimension(dim,len(infile.dimensions[dim]))
            outfile.createVariable(fieldname,datatype='f',dimensions=list(outfile.dimensions.keys()),fill_value=field._FillValue)
            bot_field2D = outfile.variables[fieldname]
    elif len(field.dimensions) == 3:
        if out3D:
            for dim in indim:
                outfile.createDimension(dim,indimlen[dim])
            outfile.createVariable(fieldname,datatype='f',dimensions=list(outfile.dimensions.keys()),fill_value=field._FillValue)
            bot_field3D = outfile.variables[fieldname]
            bot_field2D = ma.zeros(tuple(indimlen[dim] for dim in [indim[1],indim[2]]))
        else:
            for dim in indim[1],indim[2]:
                outfile.createDimension(dim,len(infile.dimensions[dim]))
            outfile.createVariable(fieldname,datatype='f',dimensions=list(outfile.dimensions.keys()),fill_value=field._FillValue)
            bot_field2D = outfile.variables[fieldname]
    else:
        raise Exception("ERROR : can only handle 4D or 3D fields!")

    if out3D:
        bot_field = bot_field3D
    else: 
        bot_field = bot_field2D

    if level:
        outfile.createVariable('bottom_level',datatype='i',dimensions=list(outfile.dimensions.keys())[-2:],fill_value=field._FillValue)
        bot_level = outfile.variables['bottom_level']
    else:
        bot_level = ma.zeros(bot_field2D[:].shape)

    bot_level[:], bot_field2D[:] = bottom_field(infield=infile.variables[fieldname][:],wgrid=wgrid)

    if out3D:
        if len(field.dimensions) == 4:
            bot_field3D[:] = (bot_field2D[:]*masked_ones[:].swapaxes(0,1)).swapaxes(0,1)
        else:
            bot_field3D[:] = bot_field2D[:]*masked_ones[:]

    # Copy over the attributes adjusting the coordinates attribute as appropriate.
    attrdict={}
    for attr in field.ncattrs():
        attrdict[attr.encode('ascii')] = field.getncattr(attr)
    # Don't include the _FillValue attribute in the dictionary as this requires special treatment in createVariable.
    print('field.ncattrs() : ', [attr.encode('ascii') for attr in field.ncattrs()])
    if b'_FillValue' in [attr.encode('ascii') for attr in field.ncattrs()]:
        print('hello!')
        del attrdict[b'_FillValue']

    bot_field.setncatts(attrdict)
    if not out3D:
        # remove the depth coordinate variable from the coordinate attribute:
        coords_list = bot_field.coordinates.split()
        if len(coords_list) == 3:
            bot_field.coordinates=" ".join(coords_list[1:])
        elif len(coords_list) == 4:
            bot_field.coordinates=" ".join([coords_list[0]]+coords_list[2:])

    if level:
        bot_level.coordinates = bot_field.coordinates

    # Copy in all the coordinate and bounds variables and their attributes

    if out3D:
        for depvar in 'depth','deptht','depthu','depthv':
            try:
                invar=infile.variables[depvar]
            except KeyError:
                pass
            else:
                break
        else:
            raise Exception('Could not find depth coordinate in input file.')
        varlist = ['nav_lat','nav_lon',depvar]
    else:
        varlist = ['nav_lat','nav_lon']

    for var in varlist:
        invar=infile.variables[var]
        outfile.createVariable(var,datatype='f',dimensions=invar.dimensions)
        outfile.variables[var][:] = invar[:]
        attrdict={}
        for attr in invar.ncattrs():
            attrdict[attr.encode('ascii')] = invar.getncattr(attr)
        outfile.variables[var].setncatts(attrdict)

         


    infile.close()
    outfile.close()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="name of input file")
    parser.add_argument("fieldname", help="name of field to extract")
    parser.add_argument("-W", "--wgrid", action="store_true",dest="wgrid",
                    help="W-grid field")
    parser.add_argument("-L", "--level", action="store_true",dest="level",
                    help="Output index of bottom level as well as bottom field value")
    parser.add_argument("--3D", action="store_true",dest="out3D",
                    help="Output bottom field projected onto masked 3D field.")

    args = parser.parse_args()

    bottom_field_wrapper(args.filename,args.fieldname,wgrid=args.wgrid,level=args.level,out3D=args.out3D)
