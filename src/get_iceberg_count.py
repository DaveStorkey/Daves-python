#!/usr/bin/env python

# Extract the number of icebergs in each grid cell in an iceberg restart file
# and write out as a field. 
#
# Adapted from a routine by Dan Copsey. 
# DS. Nov 2021
#

import netCDF4 as nc
import numpy as np
import pdb

def get_iceberg_count(infile=None,outfile=None):

    rootgrp = nc.Dataset(infile, "r")

    try:
        nimpp, njmpp = rootgrp.DOMAIN_position_first
    except AttributeError:
        nimpp, njmpp = 0, 0

    stored_heat = rootgrp.variables['stored_heat'][:]
    i_array = rootgrp.variables['xi'][:]
    j_array = rootgrp.variables['yj'][:]
    lon_array = rootgrp.variables['lon'][:]
    lat_array = rootgrp.variables['lat'][:]

    print('Stored heat shape = ',stored_heat.shape)

    ny = stored_heat.shape[0]
    nx = stored_heat.shape[1]
    iceberg_count = np.zeros((ny,nx))
    nn = i_array.shape[0]

    for n in range(nn):
        gi = int(np.rint(i_array[n]))
        gj = int(np.rint(j_array[n]))

        i = gi - nimpp - 2
        j = gj - njmpp - 2

        if i == 8 and j == 1 and n < 100:
            print('Lon = {0}, Lat = {1}'.format(lon_array[n], lat_array[n]))

        try:
            iceberg_count[j,i] = iceberg_count[j,i]+1
        except:
            pdb.set_trace()

    with nc.Dataset(outfile,"w") as dataout:
        for dim,dimlen in zip(("y","x"),(ny,nx)):
            dataout.createDimension(dim,dimlen)
        dataout.createVariable("n_icebergs",datatype='i',dimensions=list(dataout.dimensions.keys()),
                                fill_value=-9999)
        dataout.variables["n_icebergs"][:] = iceberg_count[:]
                    
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", action="store",dest="infile",
                    help="name of input file")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",
                    help="name of output file")

    args = parser.parse_args()

    get_iceberg_count(infile=args.infile,outfile=args.outfile)
