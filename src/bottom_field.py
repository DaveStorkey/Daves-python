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
import sys

def uvmask_from_tmask(tmask):
    # Function to calculate umask and vmask fields given a tmask field. 
    # These are the type of masks used in NEMO, so 1=sea and 0=land.
    # The umask and vmask fields output from this function are *not*
    # the same as the umask and vmask fields used internally in NEMO or
    # written to the mesh-mask file because here bathymetry edge points
    # are masked whereas for the standard NEMO fields edge points are
    # unmasked and the corresponding normal velocities set to zero.

    umask = np.minimum(tmask[:],np.roll(tmask[:],-1,-1))
    vmask = np.minimum(tmask[:],np.roll(tmask[:],-1,-2))

    return umask,vmask

def bottom_field(infield=None,tmask=None,grid=None):

    if grid is None:
        grid="T"

    if grid == "U" or grid == "V":
        if tmask is not None :
            umask,vmask = uvmask_from_tmask(tmask)
        else:
            raise Exception("bottom_field: Error - for U-grid and V-grid fields must supply a tmask field.")

    if grid == "U":
        fieldmask = umask
    elif grid == "V":
        fieldmask = vmask

    if len(infield.shape) == 4:
        bottom_level = ma.zeros([infield.shape[i] for i in (0,2,3)])
        bottom_field = ma.zeros([infield.shape[i] for i in (0,2,3)])
        if grid == "W":
            # Special treatment required for W-grid fields, because field is not masked at the top of the bathymetry
            # and at least for W itself takes zero values here which aren't very interesting. We want to see the values
            # just above the bathymetry.
            for lev in [z+1 for z in range(len(infield[0,:,0,0])-2)]:
                bottom_level[:,:,:] = np.where( ~infield[:,lev,:,:].mask & infield[:,lev+1,:,:].mask, 
                                                  lev-1, bottom_level )
                bottom_field[:,:,:] = np.where( ~infield[:,lev,:,:].mask & infield[:,lev+1,:,:].mask, 
                                                  infield[:,lev-1,:,:], bottom_field )
        elif grid == "U" or grid == "V":
            # For U-grid and V-grid fields need to use the non-standard derived umask/vmask fields.
            for lev in range(len(infield[0,:,0,0])-1):
                bottom_level[:,:,:] = np.where( fieldmask[0,lev  ,:,:] & ~fieldmask[0,lev+1,:,:], 
                                                  lev, bottom_level )
                bottom_field[:,:,:] = np.where( fieldmask[0,lev  ,:,:] & ~fieldmask[0,lev+1,:,:], 
                                                  infield[:,lev,:,:], bottom_field )
        else:
            for lev in range(len(infield[0,:,0,0])-1):
                bottom_level[:,:,:] = np.where( ~infield[:,lev,:,:].mask & infield[:,lev+1,:,:].mask, 
                                                  lev, bottom_level )
                bottom_field[:,:,:] = np.where( ~infield[:,lev,:,:].mask & infield[:,lev+1,:,:].mask, 
                                                 infield[:,lev,:,:], bottom_field )
        
    elif len(infield.shape) == 3:
        bottom_level = ma.zeros([infield.shape[i] for i in (1,2)])
        bottom_field = ma.zeros([infield.shape[i] for i in (1,2)])
        if grid == "W":
            # Special treatment required for W-grid fields, because field is not masked at the top of the bathymetry
            # and at least for W itself takes zero values here which aren't very interesting. We want to see the values
            # just above the bathymetry.
            for lev in [z+1 for z in range(len(infield[:,0,0])-2)]:
                bottom_level[:,:] = np.where( ~infield[lev,:,:].mask & infield[lev+1,:,:].mask,
                                                  lev-1, bottom_level )
                bottom_field[:,:] = np.where( ~infield[lev,:,:].mask & infield[lev+1,:,:].mask, 
                                                 infield[lev-1,:,:], bottom_field )
        elif grid == "U" or grid == "V":
            # For U-grid fields need to use the non-standard derived umask/vmask fields.
            for lev in range(len(infield[:,0,0])-1):
                bottom_level[:,:] = np.where( fieldmask[0,lev  ,:,:] & ~fieldmask[0,lev+1,:,:],
                                                  lev, bottom_level )
                bottom_field[:,:] = np.where( fieldmask[0,lev  ,:,:] & ~fieldmask[0,lev+1,:,:],
                                                 infield[lev,:,:], bottom_field )
        else: # T-grid default
            for lev in range(len(infield[:,0,0])-1):
                bottom_level[:,:] = np.where( ~infield[lev,:,:].mask & infield[lev+1,:,:].mask, 
                                                  lev, bottom_level )
                bottom_field[:,:] = np.where( ~infield[lev,:,:].mask & infield[lev+1,:,:].mask, 
                                                 infield[lev,:,:], bottom_field )
    else:
        raise Exception("ERROR : can only handle 4D or 3D fields!")

    return bottom_level,bottom_field


def bottom_field_wrapper(filename,fieldname,maskfile=None,grid=None,level=False,out3D=False):

    if grid is None:
        grid="T"

    if maskfile is not None:
        with nc.Dataset(maskfile,"r") as inmask:
            try:
                tmask = inmask.variables["tmask"][:]
            except KeyError:
                raise Exception("Error: Could not find tmask field in "+maskfile)
    elif grid == "U" or grid == "V":
        raise Exception("Error: mask file (containing tmask field) required for U- or V-grid fields.")
    else:
        tmask=None

    if out3D and (grid == "U" or grid == "V"):
        umask,vmask = uvmask_from_tmask(tmask)

    infile = nc.Dataset(filename,"r")
    field = infile.variables[fieldname]
    if out3D:
        masked_ones = ma.ones(field.shape)
        masked_ones.mask = field[:].mask
        if grid == "U":
            masked_ones.mask[:] = 1 - umask[:]
        elif grid == "V":
            masked_ones.mask[:] = 1 - vmask[:]
    indim = field.dimensions
    indimlen = { dim : len(infile.dimensions[dim]) for dim in indim }
    if out3D:
        outfile = nc.Dataset(filename.replace(".nc","_bottom3D.nc"),"w")
    else:
        outfile = nc.Dataset(filename.replace(".nc","_bottom.nc"),"w")

    if fieldname in ["vovecrtz"]:
        wgrid=True

    if len(field.dimensions) == 4:
        if out3D:
            for dim in indim:
                outfile.createDimension(dim,indimlen[dim])
            outfile.createVariable(fieldname,datatype="f",dimensions=list(outfile.dimensions.keys()),fill_value=field._FillValue)
            bot_field3D = outfile.variables[fieldname]
            bot_field2D = ma.zeros(tuple(indimlen[dim] for dim in [indim[0],indim[2],indim[3]]))
        else:
            for dim in indim[0],indim[2],indim[3]:
                outfile.createDimension(dim,len(infile.dimensions[dim]))
            outfile.createVariable(fieldname,datatype="f",dimensions=list(outfile.dimensions.keys()),fill_value=field._FillValue)
            bot_field2D = outfile.variables[fieldname]
    elif len(field.dimensions) == 3:
        if out3D:
            for dim in indim:
                outfile.createDimension(dim,indimlen[dim])
            outfile.createVariable(fieldname,datatype="f",dimensions=list(outfile.dimensions.keys()),fill_value=field._FillValue)
            bot_field3D = outfile.variables[fieldname]
            bot_field2D = ma.zeros(tuple(indimlen[dim] for dim in [indim[1],indim[2]]))
        else:
            for dim in indim[1],indim[2]:
                outfile.createDimension(dim,len(infile.dimensions[dim]))
            outfile.createVariable(fieldname,datatype="f",dimensions=list(outfile.dimensions.keys()),fill_value=field._FillValue)
            bot_field2D = outfile.variables[fieldname]
    else:
        raise Exception("ERROR : can only handle 4D or 3D fields!")

#    if out3D and (grid == "U" or grid == "V"):
#        outmask = nc.Dataset("uvmask.nc","w")
#        dims=["z","y","x"]
#        dimlens=umask.squeeze().shape
#        for dim,dimlen in zip(dims,dimlens):
#            outmask.createDimension(dim,dimlen)
#        outmask.createVariable("fieldmask",datatype="b",dimensions=list(outmask.dimensions.keys()))
#        outmask["fieldmask"][:] = ~field[:].mask.squeeze()
#        outmask.createVariable("onesmask",datatype="b",dimensions=list(outmask.dimensions.keys()))
#        outmask["onesmask"][:] = masked_ones[:].mask
#        outmask.createVariable("umask",datatype="b",dimensions=list(outmask.dimensions.keys()))
#        outmask["umask"][:] = umask[:].squeeze()
#        outmask.createVariable("vmask",datatype="b",dimensions=list(outmask.dimensions.keys()))
#        outmask["vmask"][:] = vmask[:].squeeze()

    if out3D:
        bot_field = bot_field3D
    else: 
        bot_field = bot_field2D

    if level:
        outfile.createVariable("bottom_level",datatype="i",dimensions=list(outfile.dimensions.keys())[-2:],fill_value=field._FillValue)
        bot_level = outfile.variables["bottom_level"]
    else:
        bot_level = ma.zeros(bot_field2D[:].shape)

    bot_level[:], bot_field2D[:] = bottom_field(infield=infile.variables[fieldname][:],grid=grid,tmask=tmask)

    if out3D:
        if len(field.dimensions) == 4:
            bot_field3D[:] = (bot_field2D[:]*masked_ones[:].swapaxes(0,1)).swapaxes(0,1)
        else:
            bot_field3D[:] = bot_field2D[:]*masked_ones[:]

    # Copy over the attributes adjusting the coordinates attribute as appropriate.
    attrdict={}
    for attr in field.ncattrs():
        attrdict[attr.encode("ascii")] = field.getncattr(attr)
    # Don"t include the _FillValue attribute in the dictionary as this requires special treatment in createVariable.
    print("field.ncattrs() : ", [attr.encode("ascii") for attr in field.ncattrs()])
    if b"_FillValue" in [attr.encode("ascii") for attr in field.ncattrs()]:
        print("hello!")
        del attrdict[b"_FillValue"]

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

    varlist=[]
    for var in "time_centered","time","time_counter","time_instant":
        try:
            invar=infile.variables[var]
        except KeyError:
            pass
        else:
            break
    else:
        raise Exception("Could not find time coordinate in input file.")
    varlist+=[var]
    for var in "nav_lon","nav_lon_grid_T","nav_lon_grid_U","nav_lon_grid_V":
        try:
            invar=infile.variables[var]
        except KeyError:
            pass
        else:
            break
    else:
        raise Exception("Could not find longitude coordinate in input file.")
    varlist+=[var]
    for var in "nav_lat","nav_lat_grid_T","nav_lat_grid_U","nav_lat_grid_V":
        try:
            invar=infile.variables[var]
        except KeyError:
            pass
        else:
            break
    else:
        raise Exception("Could not find latitude coordinate in input file.")
    varlist+=[var]
    if out3D:
        for var in "depth","deptht","depthu","depthv":
            try:
                invar=infile.variables[var]
            except KeyError:
                pass
            else:
                break
        else:
            raise Exception("Could not find depth coordinate in input file.")
        varlist+=[var]
        for var in "e3t","e3u","e3v","thkcello":
            try:
                invar=infile.variables[var]
            except KeyError:
                pass
            else:
                break
        else:
            raise Exception("Could not find thickness variable in input file.")
        varlist+=[var]

    for var in varlist:
        invar=infile.variables[var]
        attrdict={}
        for attr in invar.ncattrs():
            attrdict[attr.encode("ascii")] = invar.getncattr(attr)
        if b"_FillValue" in attrdict:
            outfile.createVariable(var,datatype="f",dimensions=invar.dimensions,fill_value=invar._FillValue)
            del attrdict[b"_FillValue"]
        else:
            outfile.createVariable(var,datatype="f",dimensions=invar.dimensions)
        outfile.variables[var][:] = invar[:]
        outfile.variables[var].setncatts(attrdict)

    infile.close()
    outfile.close()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="name of input file")
    parser.add_argument("fieldname", help="name of field to extract")
    parser.add_argument("-M", "--maskfile", action="store",dest="maskfile",
                    help="Name of file containing tmask field (usually mesh mask file).")
    parser.add_argument("-G", "--grid", action="store",dest="grid",
                    help="Which grid is field on? T (default), U, V, W. If U or V, meshmask file required.")
    parser.add_argument("-L", "--level", action="store_true",dest="level",
                    help="Output index of bottom level as well as bottom field value")
    parser.add_argument("--3D", action="store_true",dest="out3D",
                    help="Output bottom field projected onto masked 3D field.")

    args = parser.parse_args()

    bottom_field_wrapper(args.filename,args.fieldname,grid=args.grid,maskfile=args.maskfile,level=args.level,out3D=args.out3D)
