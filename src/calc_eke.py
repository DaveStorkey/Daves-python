#! /usr/bin/env python

'''
Created in Oct 2013

@author: Dave Storkey
'''

import sys
import iris
import numpy as np

def calc_eke(infileU,infileV,outfile):

    u=iris.load_cube(infile,'')

    rho = tem.copy()
    rho.data[:,:,:,:] = 0.0
    rho.var_name='sigma_theta'
    rho.standard_name='sea_water_sigma_theta'
    rho.long_name='sea_water_sigma_theta'
    rho.units='kg/m3'

    print 'tem.data.shape : ',tem.data.shape
    print 'tem.data.shape[-1] : ',tem.data.shape[-1]

    partsize = 1 + tem.data.shape[-1]/int(nsplit)
    print 'partsize is ',partsize

    istart = 0
    iend = partsize

    for split in range(int(nsplit)):
        print 'istart is : ',istart
        print 'iend is : ',iend
        rho_split = calc_sigma0(tem.data[0,:,:,istart:iend],sal.data[0,:,:,istart:iend])
        rho.data[0,:,:,istart:iend] = rho_split
        istart = iend
        iend = iend + partsize    
        if iend > tem.data.shape[-1]: iend = tem.data.shape[-1]

    iris.save(rho,outfile)

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("infileU", help="name of U-grid input file")
    parser.add_argument("infileV", help="name of V-grid input file")
    parser.add_argument("outfile", help="name of output file")

    args = parser.parse_args()

    calc_eke(args.infileU, args.infileV, args.outfile)
