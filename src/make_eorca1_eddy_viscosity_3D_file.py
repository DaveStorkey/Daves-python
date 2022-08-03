#! /usr/bin/env python

'''
Routine to create "eddy_viscosity_3D.nc" file required at NEMO 4.0
for eORCA1. 

@author: Dave Storkey
@date: Nov 2017
'''
import netCDF4 as nc
import numpy as np

def make_eddy_viscosity_3D(ahm0=None, aht0=None, ahmcoef_file=None, domcfg_file=None):

    with nc.Dataset(domcfg_file,'r') as domcfg_in:
        lon = domcfg_in.variables['nav_lon'][:]
        lat = domcfg_in.variables['nav_lat'][:]
        gdept_1d = domcfg_in.variables['gdept_1d'][0][:]

    with nc.Dataset(ahmcoef_file,'r') as ahmcoef_in:
        ahmcoef = ahmcoef_in.variables['icof'][:]

    zcoft = 0.01 * ahmcoef[:]
    print 'gdept_1d : ',gdept_1d[:]

    nx = ahmcoef.shape[1]
    ny = ahmcoef.shape[0]
    nz = gdept_1d.shape[0]

    outfile=nc.Dataset('eddy_viscosity_3D.nc','w')
    outfile.createDimension('x',size=nx)
    outfile.createDimension('y',size=ny)
    outfile.createDimension('z',size=nz)
    outfile.createVariable('nav_lon',datatype='f',dimensions=('y','x'))
    outfile.createVariable('nav_lat',datatype='f',dimensions=('y','x'))
    outfile.variables['nav_lon'][:]=lon[:]
    outfile.variables['nav_lat'][:]=lat[:]
    outfile.createVariable('ahmt_3d',datatype='f',dimensions=('z','y','x'))
    outfile.createVariable('ahmf_3d',datatype='f',dimensions=('z','y','x'))
    outfile.variables['ahmt_3d'][:]=ahm0
    outfile.variables['ahmf_3d'][:]=ahm0

    ahmt_3d = outfile.variables['ahmt_3d']
    lat3d = lat[:] * np.ones(ahmt_3d.shape)

    # horizontal variation
    ahmt_3d[:] = np.where(np.absolute(lat3d[:]) <= 2.5, aht0, ahm0)
    ahmt_3d[:] = np.where( (np.absolute(lat3d[:]) > 2.5) & (np.absolute(lat3d[:]) <= 20.0), 
            aht0 + 0.5*(ahm0-aht0)*(1-np.cos(np.radians(np.absolute(lat[:])-2.5)*180.0/17.5)), ahmt_3d[:] )
    ahmt_3d[:] = zcoft[:] * ahm0 + (1.-zcoft[:]) * ahmt_3d[:]

    # vertical variation
    zcoef = 1.0 + np.rint(9.0*(gdept_1d[:]-800.0)/(3000.0-800.0))
    zcoef[:] = np.clip(zcoef[:],1.0,10.0)
    print 'zcoef : ',zcoef[:]
    ahmtop_3d = ahmt_3d[0,:,:] * np.ones(ahmt_3d.shape)
    ahmt_3d[:,:,:] = np.transpose( np.minimum( ahm0, zcoef[:] * np.transpose(ahmtop_3d[:,:,:]) ) )

    outfile.close()


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--ahm0", action="store",dest="ahm0",type=float,default=None,
                    help="ahm0")
    parser.add_argument("-t", "--aht0", action="store",dest="aht0",type=float,default=None,
                    help="aht0")
    parser.add_argument("-d", "--domcfg", action="store",dest="domcfg_file",default=None,
                    help="domcfg file")
    parser.add_argument("-c", "--ahmcoef", action="store",dest="ahmcoef_file",default=None,
                    help="ahmcoef file")

    args = parser.parse_args()

    make_eddy_viscosity_3D(ahm0=args.ahm0, aht0=args.aht0, ahmcoef_file=args.ahmcoef_file, 
                           domcfg_file=args.domcfg_file)
