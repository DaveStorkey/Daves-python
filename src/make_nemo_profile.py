#! /usr/bin/env python

'''
Routine to generate horizontally averaged vertical profiles of NEMO fields.
Could in future do different kinds of aggregation.

Assumes that all the fields in the input file are on the same grid. 

Split out of the original version of plot_nemo_profile.py.

@author: Dave Storkey
@date: March 2019
'''
import netCDF4 as nc
import numpy as np
import numpy.ma as ma

def make_nemo_profile(infile=None, outfile=None, meshfile=None, fieldnames=None,  
                      maskfile=None, maskname=None, depthvar=None):

    if outfile is None:
        outfile=infile.replace('.nc','_profile.nc')

    infields=[]
    with nc.Dataset(infile,'r') as indata:
        for count, fieldname in enumerate(fieldnames):
            infields.append(indata.variables[fieldname][0])
            print 'min/max infields[-1][:] : ', ma.min(infields[-1][:]), ma.max(infields[-1][:])
            print 'ma.count_masked(infields[-1][:]) : ',ma.count_masked(infields[-1][:])
            if ma.count_masked(infields[-1][:]) == 0:
                raise Exception("Input field from "+datafile+"is not masked (has no _FillValues). All input fields must be masked.")

            fieldattr={}
            for attr in indata.variables[fieldname].ncattrs():
                fieldattr[attr.encode('ascii')] = indata.variables[fieldname].getncattr(attr)
            if '_FillValue' in fieldattr.keys():
                del fieldattr['_FillValue']

            if count == 0:
                # get depth coordinate from first file
                if depthvar is not None:
                    try:
                        depth = indata.variables[depthvar][:]
                    except KeyError:
                        raise Exception('Error: could not find coordinate variable '+depthvar+' in file.')
                else:
                    for depthvar in ['depth','deptht','depthu','depthv','depthw']:
                        try:
                            depth=indata.variables[depthvar][:]
                        except KeyError:
                            pass
                        else:
                            break
                    else:
                        raise Exception('Error: Could not find depth variable in file.')
                depthattr={}
                for attr in indata.variables[depthvar].ncattrs():
                    depthattr[attr.encode('ascii')] = indata.variables[depthvar].getncattr(attr)

                # get cell dimensions from first file if mesh_mask not specified 
                if meshfile is None:
                    try:
                        cell_area = indata.variables['area'][:]
                    except KeyError:
                        raise Exception('Error: Could not find cell area in '+datafile)
                    for thickname in ['thkcello', 'e3t', 'e3u', 'e3v', 'e3f', 'e3w']:
                        try:
                            cell_thick_3d = indata.variables[thickname][:]
                        except KeyError:
                            pass
                        else:
                            break
                    else:
                        raise Exception('Error: Could not find cell thickness in '+datafile)

    if maskfile is not None:
        if maskname is None:
            maskname='mask'
        with nc.Dataset(maskfile,'r') as indata:
            mask = indata.variables[maskname][:]
    else:
        # very important to do a copy here
        mask = infields[0][:].copy()
        mask[:] = 1.0

    if meshfile is not None:
        with nc.Dataset(meshfile,'r') as incoords:
            cell_area = incoords.variables['e1t'][:] * incoords.variables['e2t'][:]
            cell_thick = incoords.variables['e3t_0'][:]
        ones3D = ma.ones(infields[0].shape)
        # broadcast 1D vertical scale factor field to 3D
        cell_thick_3d = ma.transpose(cell_thick[:] * ma.transpose(ones3D[:]))

    # apply 3D mask 
    cell_volume_3d = cell_area[:,:] * cell_thick_3d[:,:,:]
    cell_volume_3d.mask = infields[0].mask
    # xy_averaged_volume is different on different depth levels because 
    # of the different mask at different depths. 
    xy_averaged_volume = ma.average(ma.average(cell_volume_3d[:],axis=-1),axis=-1)
    # remove any redundant time dimension picked up from mesh mask file.
    if len(xy_averaged_volume.shape) == 2:
        xy_averaged_volume = xy_averaged_volume[0]

    with nc.Dataset(outfile,'w') as outdata:
        outdata.createDimension(depthvar,len(depth))
        outdata.createVariable(depthvar,datatype='f',dimensions=depthvar,fill_value=-1.0e+20)
        outdata.variables[depthvar][:] = depth[:]
        outdata.variables[depthvar].setncatts(depthattr)
        for fieldname,field in zip(fieldnames,infields):
            print 'field.shape : ',field.shape
            x_averaged_field = ma.average( mask[:] * field[:,:,:] * cell_volume_3d[:,:,:],axis=-1 )
            print 'x_averaged_field.shape : ',x_averaged_field.shape
            xy_averaged_field = ma.average(x_averaged_field[:,:],axis=-1)    
            print 'xy_averaged_field.shape : ',xy_averaged_field.shape
            profile_out = xy_averaged_field[:] / (xy_averaged_volume[:])
            if len(profile_out.shape) == 2:
                # remove degenerate time dimension.
                profile_out = profile_out[0]
            outdata.createVariable(fieldname,datatype='f',dimensions=depthvar,fill_value=-1.0e+20)
            outdata.variables[fieldname][:] = profile_out[:]            
            outdata.variables[fieldname].setncatts(fieldattr)


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", action="store",dest="infile",default=None,
                    help="name of input data file with 3D fields")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",nargs="+",default=None,
                    help="name of output data file(s) to write profiles to")
    parser.add_argument("-f", "--fieldnames", action="store",dest="fieldnames",nargs="+",default=None,
                    help="name of field(s) to average")
    parser.add_argument("-z", "--depthvar", action="store",dest="depthvar",default=None,
                    help="name of (non-standard) depth coordinate in file")
    parser.add_argument("-m", "--meshfile", action="store",dest="meshfile",default=None,
                    help="name of coordinates file (not required if cell area and thicknesses in data file).")
    parser.add_argument("-k", "--maskfile", action="store",dest="maskfile",default=None,
                    help="name of optional mask file")
    parser.add_argument("-K", "--maskname", action="store",dest="maskname",default=None,
                    help="name of optional mask field")

    args = parser.parse_args()

    make_nemo_profile( infile=args.infile, outfile=args.outfile, meshfile=args.meshfile, 
                       fieldnames=args.fieldnames, maskfile=args.maskfile, maskname=args.maskname,
                       depthvar=args.depthvar)



