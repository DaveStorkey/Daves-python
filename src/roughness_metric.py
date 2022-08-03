#! /usr/bin/env python
'''
Play with roughness metrics based on FFTs

Dave Storkey
Feb 2021
'''

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

def roughness(array_in, ncutoff=None, view=None, vmax=None):

    fftsizeX=np.int(array_in.shape[1]/2)
    fftsizeY=np.int(array_in.shape[0]/2)

    fft2d = np.abs(np.fft.rfft2(array_in))
    # keep only positive wavenumbers
    xx = np.arange(fftsizeX)
    yy = np.arange(fftsizeY)
    fft2d = fft2d[0:fftsizeY,0:fftsizeX]

    if view:
        plt.pcolormesh(xx,yy,fft2d,vmin=0.0,vmax=vmax)
        plt.colorbar()
        plt.show()

    xx2d = np.ones((fftsizeY,fftsizeX))*xx
    yy2d = (np.ones((fftsizeY,fftsizeX)).transpose()*yy).transpose()

    # factors of 2.pi here??
    xxi = fftsizeX*2/xx2d
    yyi = fftsizeY*2/yy2d
    ncutoff = 5.0

    fft2d_smallscale = np.where(np.sqrt(xxi*xxi+yyi*yyi) < ncutoff,fft2d,0.0)

    if view:
        plt.pcolormesh(xx,yy,fft2d_smallscale,vmin=0.0)
        plt.colorbar()
        plt.show()

    metric = np.average(fft2d_smallscale) 

    return(metric)


def roughness_metric(filename,outfile=None,i0=None,j0=None,squaresize=None,vmax=None,map=None):

    with nc.Dataset(filename,'r') as data_in:
        bathy = data_in.variables['Bathymetry']
        try:
            fillvalue = bathy._FillValue
        except AttributeError:
            fillvalue = None
        bathy_dims = bathy.dimensions
        bathy_dimlens = [len(data_in.dimensions[bdim]) for bdim in bathy_dims]
        bathy = bathy[:]

    if len(bathy.shape) == 3:
        # get rid of degenerate time dimension if necessary
        bathy=bathy[0]
        bathy_dims=bathy_dims[1:]
        bathy_dimlens=bathy_dimlens[1:]

    if map:
        # Produce 2D map of roughness metric
        if squaresize is None:
            squaresize=11
        halfsquare = np.int(squaresize/2)
        (ny,nx) = bathy.shape
      
        metric = np.zeros((ny,nx))
        for iy in halfsquare + np.arange(ny-squaresize-1):
            print('iy : ',iy)
            for ix in halfsquare + np.arange(nx-squaresize-1):
                metric[iy,ix] = roughness(bathy[iy-halfsquare:iy+halfsquare,ix-halfsquare:ix+halfsquare],ncutoff=3.0)

        if outfile:
            with nc.Dataset(outfile,'w') as dataout:
                for dim,dimlen in zip(bathy_dims,bathy_dimlens):
                    dataout.createDimension(dim,dimlen)
                if fillvalue is not None:
                    dataout.createVariable('roughness',datatype='f',dimensions=list(dataout.dimensions.keys()),
                                            fill_value=fillvalue)
                else:
                    dataout.createVariable('roughness',datatype='f',dimensions=list(dataout.dimensions.keys()))
                dataout.variables['roughness'][:] = metric[:]            

        plt.pcolormesh(metric)
        plt.colorbar()
        plt.show()

    else:
        # Look at FFTs and roughness metric for specific point
        if i0 is None:
            i0 = 210
        if j0 is None:
            j0 = 150
        if squaresize is None:
            squaresize=100

        metric = roughness(bathy[j0:j0+squaresize,i0:i0+squaresize],ncutoff=5.0,
                           view=True,vmax=vmax)    
        print('metric value : ',metric)


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="name of input file")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",
                    help="output filename for mapped roughness")
    parser.add_argument("--i0", action="store",dest="i0",type=int,
                    help="starting i-index")
    parser.add_argument("--j0", action="store",dest="j0",type=int,
                    help="starting j-index")
    parser.add_argument("-M", "--map", action="store_true",dest="map",
                    help="true => produce 2D map of roughness metric")
    parser.add_argument("-Q", "--squaresize", action="store",dest="squaresize",type=int,
                    help="size of sample square in gridpoints")
    parser.add_argument("-V", "--vmax", action="store",dest="vmax",type=float,
                    help="maximum value of FFT to plot")

    args = parser.parse_args()

    roughness_metric(args.filename,outfile=args.outfile,i0=args.i0,j0=args.j0,squaresize=args.squaresize,
                     vmax=args.vmax, map=args.map )
