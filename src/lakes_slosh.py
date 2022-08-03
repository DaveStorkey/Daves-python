#! /usr/bin/env python
###!/opt/python/gnu/2.7.9/bin/python ## for supercomputer
'''
Move water between a group of lakes or inland seas in order to equalise the SSH values. 
Done as a correction for an SSH drift which has resulted from defining a single closea
rectangle for a group of lakes/inland seas in NEMO.

Created March 2016

@author: Dave Storkey
'''

import netCDF4 as nc
import numpy as np
import numpy.ma as ma

def lakes_slosh(filename,coordinates=None,varlist=None,north=None,south=None,east=None,west=None):

    ncfile=nc.Dataset(filename, mode='r+')

    for var in varlist:

        print ''
        ssh_in = ncfile.variables[var]

        # No masking in restarts or coordinates files. Assume points where SSH=0.0 are land points
        # and mask accordingly:
        ssh_masked = ma.MaskedArray(ssh_in[0])
        ssh_masked.mask=False
        ssh_masked.mask[np.where(ssh_masked==0.0)]=True

        coordsfile=nc.Dataset(coordinates, mode='r')
        e1t_masked = ma.MaskedArray(coordsfile.variables['e1t'])
        e1t_masked.mask=False
        e1t_masked.mask[np.where(ssh_masked==0.0)]=True
        e2t_masked = ma.MaskedArray(coordsfile.variables['e2t'])
        e2t_masked.mask=False
        e2t_masked.mask[np.where(ssh_masked==0.0)]=True
  
        ssh_cutouts=[]
        e1t_cutouts=[]
        e2t_cutouts=[]
        for j2,j1,i2,i1 in zip(north,south,east,west):
            # Convert from fortran to python counting. Note that ranges in python array slicing
            # *don't* include the second value specified. Hence we only "convert" the west and
            # south values.
            i1=i1-1; i2=i2; j1=j1-1; j2=j2; 
            ssh_cutouts.append(ssh_masked[j1:j2,i1:i2])
            e1t_cutouts.append(e1t_masked[j1:j2,i1:i2])
            e2t_cutouts.append(e2t_masked[j1:j2,i1:i2])
    
        area_list=[]
        sshmean_list=[]
        for ssh,e1t,e2t in zip(ssh_cutouts,e1t_cutouts,e2t_cutouts):
            # Calculate area for each lake
            area_list.append( ma.sum(e1t*e2t) )                 
            # Calculate area-weighted mean SSH for each lake    
            sshmean_list.append( ma.sum(ssh*e1t*e2t)/ma.sum(e1t*e2t) )

        # Calculate new mean SSH as area-weighted average of individual mean SSH values:
        sshnew = sum([sshmean*area for (sshmean, area) in zip(sshmean_list, area_list)])/sum(area_list)

        # Add/subtract relevant quantity to SSH field in each lake
        for j2,j1,i2,i1,ssh,sshmean in zip(north,south,east,west,ssh_cutouts,sshmean_list):
            print 'sshmean, sshnew : ',sshmean,sshnew
            # convert from fortran to python counting (see note above for similar procedure):
            i1=i1-1; i2=i2; j1=j1-1; j2=j2; 
            ssh = ssh + sshnew - sshmean
            ssh.mask=False
            ssh_in[0,j1:j2,i1:i2] = ssh

    ncfile.close()


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="name of input file")
    parser.add_argument("-c", "--coordinates", action="store",dest="coords",
                    help="coordinates file (or file containing e1t oe e2t")
    parser.add_argument("-v", "--var", action="store",dest="varlist",nargs='+',
                    help="name of SSH variables in file")
    parser.add_argument("-N", "--north", action="store",dest="north",type=int,default=None,nargs='+',
                    help="northmost j-value of rectangles")
    parser.add_argument("-S", "--south", action="store",dest="south",type=int,default=None,nargs='+',
                    help="southmost j-value of rectangles")
    parser.add_argument("-E", "--east", action="store",dest="east",type=int,default=None,nargs='+',
                    help="eastmost i-value of rectangles")
    parser.add_argument("-W", "--west", action="store",dest="west",type=int,default=None,nargs='+',
                    help="westmost i-value of rectangles")

    args = parser.parse_args()

    lakes_slosh(args.filename,coordinates=args.coords,varlist=args.varlist,
                north=args.north,south=args.south,east=args.east,west=args.west)        
    
