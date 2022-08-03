#! /usr/bin/env python

'''
Routine to calculate average trend field given the start and end
state from a model run. 

Apr 2021 : Add U and V trends and update to python 3.

@author: Dave Storkey
@date: Oct 2016
'''

import netCDF4 as nc
import numpy as np
import numpy.ma as ma

def model_trend(initfile=None,initrestart=None,endrestart=None,years=None,days=None,
                thickwgt=None,extensive=None):

    if extensive:
        tem_factor = 1026.0 * 3991.86795711963
        sal_factor = 1026.0 * 0.001
    else:
        tem_factor = 1.0
        sal_factor = 1.0

    if initrestart is not None:
        startfile = nc.Dataset(initrestart,'r')
        temname_start = 'tn'
        salname_start = 'sn'
        uvelname_start = 'un'
        vvelname_start = 'vn'
        thicknameT_start = 'fse3t_n'
        thicknameU_start = 'fse3u_n'
        thicknameV_start = 'fse3v_n'
    elif initfile is not None:
        startfile = nc.Dataset(initfile,'r')
        temname_start = 'votemper'
        salname_start = 'vosaline'
        uvelname_start = 'vozocrtx'
        vvelname_start = 'vomecrty'
        thicknameT_start = 'vovvle3t'
        thicknameU_start = 'vovvle3u'
        thicknameV_start = 'vovvle3v'
    else:
        raise Exception('initial file not specified')

    if endrestart is not None:
        endfile = nc.Dataset(endrestart,'r')
        temname_end = 'tn'
        salname_end = 'sn'
        uvelname_end = 'un'
        vvelname_end = 'vn'
        thicknameT_end = 'fse3t_n'
        thicknameU_end = 'fse3u_n'
        thicknameV_end = 'fse3v_n'
    else:
        raise Exception('end file not specified')
 
    # if keywords days or years set, set scale factor to 
    # convert to change per second
    scale_factor = 1.0
    if days is not None:
        scale_factor = 1.0/(86400*days)
    elif years is not None:
        scale_factor = 1.0/(31536000*years)

    tem_start = startfile.variables[temname_start]
    sal_start = startfile.variables[salname_start]
    uvel_start = startfile.variables[uvelname_start]
    vvel_start = startfile.variables[vvelname_start]
    if thickwgt or extensive:
        thickT_start = startfile.variables[thicknameT_start]
        thickU_start = startfile.variables[thicknameU_start]
        thickV_start = startfile.variables[thicknameV_start]

    tem_end = endfile.variables[temname_end]
    sal_end = endfile.variables[salname_end]
    uvel_end = endfile.variables[uvelname_end]
    vvel_end = endfile.variables[vvelname_end]
    if thickwgt or extensive:
        thickT_end = endfile.variables[thicknameT_end]
        thickU_end = endfile.variables[thicknameU_end]
        thickV_end = endfile.variables[thicknameV_end]

    tem_trend = ma.zeros(tem_start[:].shape)
    if thickwgt or extensive:
        print('scale_factor : ',scale_factor)
        print('tem_factor : ',tem_factor)
        tem_trend[:] = ( tem_end[:]*thickT_end[:] - tem_start[:]*thickT_start[:] ) * scale_factor * tem_factor
    else:
        tem_trend[:] = ( tem_end[:] - tem_start[:] ) * scale_factor

    sal_trend = ma.zeros(sal_start[:].shape)
    if thickwgt or extensive:
        sal_trend[:] = ( sal_end[:]*thickT_end[:] - sal_start[:]*thickT_start[:] ) * scale_factor * sal_factor
    else:
        sal_trend[:] = ( sal_end[:] - sal_start[:] ) * scale_factor

    uvel_trend = ma.zeros(uvel_start[:].shape)
    if thickwgt or extensive:
        uvel_trend[:] = ( uvel_end[:]*thickU_end[:] - uvel_start[:]*thickU_start[:] ) * scale_factor 
    else:
        uvel_trend[:] = ( uvel_end[:] - uvel_start[:] ) * scale_factor

    vvel_trend = ma.zeros(vvel_start[:].shape)
    if thickwgt or extensive:
        vvel_trend[:] = ( vvel_end[:]*thickV_end[:] - vvel_start[:]*thickV_start[:] ) * scale_factor 
    else:
        vvel_trend[:] = ( vvel_end[:] - vvel_start[:] ) * scale_factor

    dims = startfile.dimensions
    if extensive:
        trend_file = nc.Dataset("model_trend_heatsalt.nc","w")
        temname_out = "opottemptend"
        salname_out = "osalttend"
        uvelname_out = uvelname_start+"_thick"
        vvelname_out = vvelname_start+"_thick"
    elif thickwgt:
        trend_file = nc.Dataset("model_trend_thickwgt.nc","w")
        temname_out = temname_start+"_thick"
        salname_out = salname_start+"_thick"
        uvelname_out = uvelname_start+"_thick"
        vvelname_out = vvelname_start+"_thick"
    else:
        trend_file = nc.Dataset("model_trend.nc","w")
        temname_out = temname_start
        salname_out = salname_start
        uvelname_out = uvelname_start
        vvelname_out = vvelname_start
    for dim in dims.keys():
        trend_file.createDimension(dim, len(dims[dim])) 

    
    trend_file.createVariable(temname_out,datatype='f',dimensions=tem_start.dimensions)
    trend_file.variables[temname_out][:] = tem_trend[:]

    trend_file.createVariable(salname_out,datatype='f',dimensions=sal_start.dimensions)
    trend_file.variables[salname_out][:] = sal_trend[:]

    trend_file.createVariable(uvelname_out,datatype='f',dimensions=uvel_start.dimensions)
    trend_file.variables[uvelname_out][:] = uvel_trend[:]

    trend_file.createVariable(vvelname_out,datatype='f',dimensions=vvel_start.dimensions)
    trend_file.variables[vvelname_out][:] = vvel_trend[:]

    startfile.close()
    endfile.close()
    trend_file.close()
    

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--initfile", action="store",dest="initfile",default=None,
                    help="name of init file from NEMO (nn_istate=1)")
    parser.add_argument("-I", "--initrestart", action="store",dest="initrestart",default=None,
                    help="name of initial restart file from NEMO")
    parser.add_argument("-E", "--endrestart", action="store",dest="endrestart",default=None,
                    help="name of end restart file from NEMO")
    parser.add_argument("-d", "--days", action="store",type=float,dest="days",default=None,
                    help="number of days in time interval")
    parser.add_argument("-y", "--years", action="store",type=float,dest="years",default=None,
                    help="number of years in time interval")
    parser.add_argument("-t", "--thickwgt", action="store_true",dest="thickwgt",default=None,
                    help="calculated thickness-weighted difference fields")
    parser.add_argument("-x", "--extensive", action="store_true",dest="extensive",default=None,
                    help="calculated difference fields of extensive quantities: heat and salt per unit area")

    args = parser.parse_args()

    model_trend( initfile=args.initfile, initrestart=args.initrestart, endrestart=args.endrestart, 
                 days=args.days, years=args.years, thickwgt=args.thickwgt, extensive=args.extensive ) 

