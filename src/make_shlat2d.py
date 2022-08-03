#! /usr/bin/env python

'''
Routine to create an shlat2d field based on the old
hardwiring of fmask points for low res global models
in NEMO 3.6 and earlier. 

@author: Dave Storkey
@date: Feb 2019

Apr 2020 : Added visfield and highlatps options. DS.
'''

import numpy as np
import numpy.ma as ma
import netCDF4 as nc

def make_shlat2d(infile=None, outfile=None, visfield=None,
                 eorca1=None, highlatps=None, southps=None,
                 antarcpenps=None):

    with nc.Dataset(infile,'r') as data_in:
        model_dims = data_in.dimensions
        nx = len(model_dims['x'])
        ny = len(model_dims['y'])
        nav_lon = data_in.variables['nav_lon'][:]
        nav_lat = data_in.variables['nav_lat'][:]
        if visfield:
            try:
                surfmask = data_in.variables['tmask'][0][0][:]
            except KeyError:
                raise Exception('Error: for visfield option need tmask in input file.')

    shlat2d = ma.zeros((ny,nx))

    if eorca1:
        # Hardwired indices for eORCA1 copied from dommsk.F90 @ NEMO 3.6
        # taking into account that fortran counts from 1 and python counts from 0
        # and also that fortran array(a:b) is equivalent to python array [a-1:b]

        ii0 = 282 - 1    ;   ii1 = 283         # Gibraltar Strait 
        ij0 = 241 - 1    ;   ij1 = 241     ;   shlat2d[ij0:ij1,ii0:ii1] = 2.  

        ii0 = 314 - 1    ;   ii1 = 315         # Bhosporus Strait 
        ij0 = 248 - 1    ;   ij1 = 248     ;   shlat2d[ij0:ij1,ii0:ii1] = 2.  

        ii0 =  48 - 1    ;   ii1 =  48         # Makassar Strait (Top) 
        ij0 = 189 - 1    ;   ij1 = 190     ;   shlat2d[ij0:ij1,ii0:ii1] = 3.  

        ii0 =  44 - 1    ;   ii1 =  44         # Lombok Strait 
        ij0 = 164 - 1    ;   ij1 = 165     ;   shlat2d[ij0:ij1,ii0:ii1] = 2.  

        ii0 =  53 - 1    ;   ii1 =  53         # Ombai Strait 
        ij0 = 164 - 1    ;   ij1 = 165     ;   shlat2d[ij0:ij1,ii0:ii1] = 2.  

        ii0 =  56 - 1    ;   ii1 =  56         # Timor Passage 
        ij0 = 164 - 1    ;   ij1 = 165     ;   shlat2d[ij0:ij1,ii0:ii1] = 2.  

        ii0 =  58 - 1    ;   ii1 =  58         # West Halmahera Strait 
        ij0 = 181 - 1    ;   ij1 = 182     ;   shlat2d[ij0:ij1,ii0:ii1] = 3.  

        ii0 =  55 - 1    ;   ii1 =  55         # East Halmahera Strait 
        ij0 = 181 - 1    ;   ij1 = 182     ;   shlat2d[ij0:ij1,ii0:ii1] = 3.  

        print("Number of nonzero values in shlat2d : ",np.count_nonzero(shlat2d))

    elif highlatps is not None:
        # Free slip at low latitudes and partial slip (shlat=1) at high latitudes. 
        # Cut-over latitude given by value of highlatps. 
        shlat2d[:] = ma.where( ma.absolute(nav_lat) > highlatps, 1.0, 0.0)
        
    elif southps is not None:
        # Free slip at low latitudes and partial slip (shlat=1) at high latitudes only 
        # in the southern hemisphere. 
        # Cut-over latitude given by value of southps. 
        shlat2d[:] = ma.where( -nav_lat > southps, 1.0, 0.0)
        
    elif antarcpenps is not None:
        # Partial slip (shlat=1) just around tip of Antarctic peninsula
        # Cut-over latitude given by value of highlatps. 
        shlat2d[:] = ma.where( ( ( nav_lat > -66.0 ) & ( nav_lat < -59.5 ) ) & 
                               ( ( nav_lon - nav_lat - 66.0 > -70.0 ) & ( nav_lon - nav_lat - 66.0 < -43.0 ) ), 1.0, 0.0)
        
    if visfield:
        shlat2d_vis = shlat2d[:] + surfmask[:]

    with nc.Dataset(outfile,'w') as data_out:
        data_out.createDimension('x', nx)
        data_out.createDimension('y', ny)
        data_out.createVariable('nav_lon',datatype='f',dimensions=('y','x'))
        data_out.variables['nav_lon'][:] = nav_lon[:]
        data_out.variables['nav_lon'].standard_name = 'longitude'
        data_out.variables['nav_lon'].units = 'degrees_east'
        data_out.createVariable('nav_lat',datatype='f',dimensions=('y','x'))
        data_out.variables['nav_lat'][:] = nav_lat[:]
        data_out.variables['nav_lat'].standard_name = 'latitude'
        data_out.variables['nav_lat'].units = 'degrees_north'
        data_out.createVariable('shlat2d',datatype='f',dimensions=('y','x'))
        data_out.variables['shlat2d'][:] = shlat2d[:]
        data_out.variables['shlat2d'].long_name = 'shlat topographic boundary condition'
        data_out.variables['shlat2d'].coordinates = 'nav_lat nav_lon'
        if visfield:
            data_out.createVariable('shlat2d_vis',datatype='f',dimensions=('y','x'))
            data_out.variables['shlat2d_vis'][:] = shlat2d_vis[:]
            data_out.variables['shlat2d_vis'].long_name = 'shlat topographic boundary condition visualisation'
            data_out.variables['shlat2d_vis'].coordinates = 'nav_lat nav_lon'

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", action="store",dest="infile",default=None,
                         help="name of input file (for dimensions and lat/lon)")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",default=None,
                         help="name of shlat2d output file")
    parser.add_argument("-v", "--visfield", action="store_true",dest="visfield",default=None,
                         help="write visualisation field to output (requires tmask or Bathymetry in input file)")
    parser.add_argument("--eorca1", action="store_true",dest="eorca1",default=None,
                         help="create shlat2d field for eORCA1 to match old hardwiring")
    parser.add_argument("--antarcpenps", action="store_true",dest="antarcpenps",default=None,
                         help="partial slip just around tip of Antarctic peninsula")
    parser.add_argument("--highlatps", action="store",dest="highlatps",type=float,default=None,
                         help="free slip at |lat| < highlatps and partial slip (shlat=1) at |lat| > highlatps")
    parser.add_argument("--southps", action="store",dest="southps",type=float,default=None,
                         help="free slip at -lat < highlatps and partial slip (shlat=1) at -lat > southps")

    args = parser.parse_args()

    make_shlat2d(infile=args.infile, outfile=args.outfile, visfield=args.visfield,
                 eorca1=args.eorca1, highlatps=args.highlatps, southps=args.southps, antarcpenps=args.antarcpenps)
