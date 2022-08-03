#! /usr/bin/env python

'''
Routine to plot horizontally averaged vertical profiles of NEMO fields.

Will either read in 3D fields, do the spatial averaging and plot the profiles 
(and optionally write the profiles to netcdf files) or read in and plot previously
calculated profiles.

Feb 2019 : Generalise it to read all the common depth coordinate names and also 
           have an option to input a non-standard depth coordinate name. Also 
           option to plot the depth coordinate increasing upwards. 

To do: Try rewriting using Iris and see if cube.collapsed gives the right answers; in particular
       does it cope with the masking correctly?

@author: Dave Storkey
@date: March 2018
'''
import matplotlib
# Hack for running in batch mode on SPICE. 
# Note this disables plt.show()
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import textwrap

def plot_nemo_profile(datafiles=None, profilefiles=None, meshfile=None, fieldname=None, outfile=None, 
                      maskfile=None, maskname=None, depthvar=None, xmin=None, xmax=None, depthmax=None, title=None,
                      labels=None, fieldlabel=None, legendloc=None, colours=None, outprofiles=None):

    if datafiles is not None and profilefiles is not None:
        raise Exception("Can only specify input 3D fields or profiles, not both.")

    fig, ax = plt.subplots()
    if colours is not None:
        ax.set_color_cycle(colours)

    print 'fieldname is ',fieldname

    if datafiles is not None:

        nfiles = len(datafiles)
        infields=[]
        for count, datafile in enumerate(datafiles):
            print 'datafile is ',datafile
            with nc.Dataset(datafile,'r') as indata:
                infields.append(indata.variables[fieldname][0])
                print 'min/max infields[-1][:] : ', ma.min(infields[-1][:]), ma.max(infields[-1][:])
                print 'ma.count_masked(infields[-1][:]) : ',ma.count_masked(infields[-1][:])
                if ma.count_masked(infields[-1][:]) == 0:
                    raise Exception("Input field from "+datafile+"is not masked (has no _FillValues). All input fields must be masked.")
                print 'ma.min(infields[-1][:]), ma.max(infields[-1][:]) : ',ma.min(infields[-1][:]),ma.max(infields[-1][:])

                if count == 0:
                    fieldattr={}
                    for attr in indata.variables[fieldname].ncattrs():
                        fieldattr[attr.encode('ascii')] = indata.variables[fieldname].getncattr(attr)
                    if '_FillValue' in fieldattr.keys():
                        del fieldattr['_FillValue']
                    # use long name and units attributes to label x-axis if not specified
                    if fieldlabel is None:
                        try:
                            fieldlabel = fieldattr['long_name']+' ('+fieldattr['units']+')'
                        except KeyError:
                            pass
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

    elif profilefiles is not None:

        nfiles = len(profilefiles)
        infields=[]
        for count, profilefile in enumerate(profilefiles):
            with nc.Dataset(profilefile,'r') as indata:
                infields.append(indata.variables[fieldname][:])
                if count == 0:
                    fieldattr={}
                    for attr in indata.variables[fieldname].ncattrs():
                        fieldattr[attr.encode('ascii')] = indata.variables[fieldname].getncattr(attr)
                    # use long name and units attributes to label x-axis if not specified
                    if fieldlabel is None:
                        try:
                            fieldlabel = fieldattr['long_name']+' ('+fieldattr['units']+')'
                        except KeyError:
                            pass
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

    else:
        raise Exception("Must specify one (and only one) of datafiles or profilefiles")

    field_min = None
    field_max = None

    for count, field in enumerate(infields):

        if datafiles is not None:
            print 'field.shape : ',field.shape
            x_averaged_field = ma.average( mask[:] * field[:,:,:] * cell_volume_3d[:,:,:],axis=-1 )
            print 'x_averaged_field.shape : ',x_averaged_field.shape
            xy_averaged_field = ma.average(x_averaged_field[:,:],axis=-1)    
            print 'xy_averaged_field.shape : ',xy_averaged_field.shape
            profile_to_plot = xy_averaged_field[:] / (xy_averaged_volume[:])
            if len(profile_to_plot.shape) == 2:
                # remove degenerate time dimension.
                profile_to_plot = profile_to_plot[0]
            if outprofiles is not None: 
                with nc.Dataset(outprofiles[count],'w') as outdata:
                    outdata.createDimension(depthvar,len(depth))
                    outdata.createVariable(depthvar,datatype='f',dimensions=depthvar,fill_value=-1.0e+20)
                    outdata.variables[depthvar][:] = depth[:]
                    outdata.createVariable(fieldname,datatype='f',dimensions=depthvar,fill_value=-1.0e+20)
                    outdata.variables[depthvar].setncatts(depthattr)
                    outdata.variables[fieldname][:] = profile_to_plot[:]            
                    outdata.variables[fieldname].setncatts(fieldattr)
        else:
            profile_to_plot = field

        if field_min is None:
            field_min = profile_to_plot.min()
        else:
            field_min=min(profile_to_plot.min(),field_min)

        if field_max is None:
            field_max = profile_to_plot.max()
        else:
            field_max=max(profile_to_plot.max(),field_max)

        print 'depth.shape : ',depth.shape
        print 'profile_to_plot.shape : ',profile_to_plot.shape
        plt.plot(profile_to_plot[:],depth[:])
        
    if xmin is None:
        xmin = field_min
    if xmax is None:
        xmax = field_max
    plt.gca().set_xlim([xmin,xmax])

    if 'positive' in depthattr.keys() and depthattr['positive'] == "up":
        # depth coord increasing upwards
        if depthmax is not None:
            plt.gca().set_ylim([0.0,depthmax])
    else:
        # default to depth coord increasing downwards
        if depthmax is not None:
            plt.gca().set_ylim([depthmax,0.0])
        else:
            plt.gca().invert_yaxis()

    if fieldlabel is not None:
        plt.gca().set_xlabel(fieldlabel)
    if "depth" in depthvar:
        plt.gca().set_ylabel('depth (m)')
    else:
        if 'units' in depthattr.keys() and depthattr['units'] != '1':
            plt.gca().set_ylabel(depthvar+' ('+depthunits+')')
        else:
            plt.gca().set_ylabel(depthvar)
    if labels is not None:
        if legendloc is None:
            legendloc=4
        plt.legend(labels[0:nfiles],loc=legendloc,fontsize='medium')

    if title is not None:
        plt.gca().set_title(textwrap.fill(title,60))

    if outfile is not None:
        matplotlib.rcParams['font.size'] =8
        plt.savefig(outfile,dpi=200)
    else:
        plt.show()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--datafiles", action="store",dest="datafiles",nargs="+",default=None,
                    help="name of input data file(s) with 3D fields")
    parser.add_argument("-p", "--profilefiles", action="store",dest="profilefiles",nargs="+",default=None,
                    help="name of input data file(s) with profiles")
    parser.add_argument("-q", "--outprofiles", action="store",dest="outprofiles",nargs="+",default=None,
                    help="name of output data file(s) to write profiles to")
    parser.add_argument("-m", "--meshfile", action="store",dest="meshfile",default=None,
                    help="name of coordinates file (not required if cell area and thicknesses in data file).")
    parser.add_argument("-k", "--maskfile", action="store",dest="maskfile",default=None,
                    help="name of optional mask file")
    parser.add_argument("-K", "--maskname", action="store",dest="maskname",default=None,
                    help="name of optional mask field")
    parser.add_argument("-f", "--fieldname", action="store",dest="fieldname",default=None,
                    help="name of field to plot")
    parser.add_argument("-F", "--fieldlabel", action="store",dest="fieldlabel",default=None,
                    help="label for field")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",default=None,
                    help="plot to output file - filetype by extension. If unset plot to GUI.")
    parser.add_argument("-z", "--depthvar", action="store",dest="depthvar",default=None,
                    help="name of (non-standard) depth coordinate in file")
    parser.add_argument("-x", "--xmin", action="store",dest="xmin",type=float,default=None,
                    help="start of x-axis (min field value to plot)")
    parser.add_argument("-X", "--xmax", action="store",dest="xmax",type=float,default=None,
                    help="end of x-axis (max field value to plot)")
    parser.add_argument("-D", "--depthmax", action="store",dest="depthmax",type=float,default=None,
                    help="maximum depth for plot")
    parser.add_argument("-l", "--labels", action="store",dest="labels",nargs="+",default=None,
                    help="labels for legend")
    parser.add_argument("-L", "--legendloc", action="store",dest="legendloc",default=None,type=int,
                    help="location of legend - see matplotlib legend command")
    parser.add_argument("-c", "--colours", action="store",dest="colours",nargs="+",default=None,
                    help="colours for lines")
    parser.add_argument("-t", "--title", action="store",dest="title",
                    help="title for plot")

    args = parser.parse_args()

    plot_nemo_profile( datafiles=args.datafiles, profilefiles=args.profilefiles, meshfile=args.meshfile, 
                       fieldname=args.fieldname, outfile=args.outfile, outprofiles=args.outprofiles,
                       xmin=args.xmin, xmax=args.xmax, depthmax=args.depthmax, fieldlabel=args.fieldlabel, 
                       labels=args.labels, legendloc=args.legendloc, maskfile=args.maskfile, maskname=args.maskname,
                       colours=args.colours, title=args.title, depthvar=args.depthvar)



