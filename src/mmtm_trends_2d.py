#! /usr/bin/env python
'''
Routine to do sanity checks of momentum trend diagnostics
from NEMO. This version for the 2D trends. 

Hardwired to expect output from a run with vector invariant
momentum advection for now.

Created August 2021

Apr 2021 : Add total trend diagnostic minus model trend residual as option. 

Aug 2021 : Generalise to do V component as well.

@author: Dave Storkey
'''

import netCDF4 as nc
import numpy as np
import numpy.ma as ma

def mmtm_trends(filename,components=None,xcpt=None,weighted=None,metrics_only=None,trendfile=None):

    if metrics_only:
        trends_file = nc.Dataset(filename,'r')
    else:
        trends_file = nc.Dataset(filename,'r+')

    if weighted:
        for grd in ['u','v']:
            try:
                total_trend = trends_file.variables[grd+'trd_tot2d_h'+grd]
            except KeyError:
                pass
            else:
                break
        else:
            raise Exception('Error : could not find total trend diagnostic utrd_tot2d_hu or vtrd_tot2d_hv in file '+filename)
    else:
        for grd in ['u','v']:
            try:
                total_trend = trends_file.variables[grd+'trd_tot2d']
            except KeyError:
                pass
            else:
                break
        else:
            raise Exception('Error : could not find total trend diagnostic utrd_tot2d or vtrd_tot2d in file '+filename)

    if weighted:
        wgt="_h"+grd
    else:
        wgt=""

    if metrics_only:

        utrd={}
        utrd['tot'] = trends_file.variables[grd+'trd_tot2d'+wgt][:]
        utrd['res'] = trends_file.variables[grd+'trd_res2d'+wgt][:]

    else:

        # 'atf' trend not included in standard list because 'tot' trend in NEMO doesn't include this.
        cpt_list = ['spg2d', 'pvo2d', 'frc2d', 'tau2d', 'tfr2d', 'bfr2d']

        if components is None:
            components = cpt_list
            if xcpt is not None:
                for cpt in xcpt:
                    try:
                        components.remove(cpt)
                    except ValueError:
                        raise Exception('Error: unrecognised xcpt : '+cpt)
        else:
            for cpt in components:
                if cpt not in cpt_list:
                    raise Exception('Error: unrecognised component : '+cpt)

        print('File : '+filename)
        print('Summing the following components : '+' '.join(components))

        utrd={}
        utrd['tot'] = trends_file.variables[grd+'trd_tot2d'+wgt]

        utrd['sum'] = ma.array(utrd['tot'][:]) 
        utrd['sum'][:] = 0.0
        utrd['res'] = ma.array(utrd['tot'][:]) 
        utrd['res'][:] = 0.0
        for cpt in components:
            utrd[cpt] = trends_file.variables[grd+'trd_'+cpt+wgt]
            utrd['sum'][:] = utrd['sum'][:] + utrd[cpt][:] 

        utrd['res'][:] = utrd['tot'][:] - utrd['sum'][:]
    
        if trendfile is not None:
            trendfiledata = nc.Dataset(trendfile,'r')
            utrd['model'] = ma.array(utrd['tot'][:])
            utrd['res2']  = ma.array(utrd['tot'][:])
            utrd['model'] = trendfiledata.variables['vozocrtx']
            utrd['res2'][:] = utrd['tot'][:] - utrd['model'][:]

    utrd_res_metric = ma.average(ma.absolute(utrd['res'][:,:-1,:-1]))
    utrd_tot_metric = ma.average(ma.absolute(utrd['tot'][:,:-1,:-1]))
    utrd_res_metric_fracerr = utrd_res_metric/utrd_tot_metric

    print('Results excluding EW wrap and northfold points: ')
    print('utrd_res_metric : ',utrd_res_metric)
    print('utrd_tot_metric : ',utrd_tot_metric)
    print('utrd_res_metric_fracerr : ',utrd_res_metric_fracerr)

    if not metrics_only:

        # Read in attributes from total trend field
        attrdict={}
        for attr in utrd['tot'].ncattrs():
            attrdict[attr] = utrd['tot'].getncattr(attr)
        # Don't include the _FillValue attribute in the dictionary as this requires special treatment in createVariable.
        if '_FillValue' in attrdict.keys():
            print('Found _FillValue in attrdict.keys()')
            del attrdict['_FillValue']
            fillvalue=True
        else:
            print("Didn't find _FillValue in attrdict.keys()")
            fillvalue=False

        output_list = ['sum','res']
        if trendfile is not None:
            output_list = output_list+['model','res2']

        for cpt in output_list:
            if not grd+'trd_'+cpt+wgt in trends_file.variables.keys():
                if fillvalue:
                    trends_file.createVariable(grd+'trd_'+cpt+wgt,datatype='f',dimensions=utrd['tot'].dimensions,fill_value=utrd['tot']._FillValue)
                else:
                    trends_file.createVariable(grd+'trd_'+cpt+wgt,datatype='f',dimensions=utrd['tot'].dimensions)
                utrd_out = trends_file.variables[grd+'trd_'+cpt+wgt]
                utrd_out.setncatts(attrdict)
                if cpt == 'sum':
                    utrd_out.long_name='Sum of U trends: '+' '.join(components)
                elif cpt == 'model':
                    utrd_out.long_name='i-trend: model end minus model start'
                elif cpt == 'res2':
                    utrd_out.long_name='Total U trend minus Model U trend'
                else:
                    utrd_out.long_name='Total U trend minus sum of U trends'
            else:
                utrd_out = trends_file.variables[grd+'trd_'+cpt+wgt]

            utrd_out[:] = utrd[cpt][:]

    trends_file.close()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="name of input file")
    parser.add_argument("-C","--components",action="store",dest="components",nargs="+", 
                    help="required components chosen from : hpg, spg, keg, rvo, pvo, zad, ldf, zdf, tau, atf, bta") 
    parser.add_argument("-X","--xcpt",action="store",dest="xcpt",nargs="+", 
                    help="components to omit from sum, chosen from : hpg, spg, keg, rvo, pvo, zad, ldf, zdf, tau, atf, bta") 
    parser.add_argument("-m","--metrics_only",action="store_true",dest="metrics_only")
    parser.add_argument("-W","--weighted",action="store_true",dest="weighted",
                    help="read thickness-weighted input trends *_hu, *_hv")
    parser.add_argument("-T","--trendfile",action="store",dest="trendfile",
                    help="name of file with model trends")

    args = parser.parse_args()

    mmtm_trends(args.filename, components=args.components, xcpt=args.xcpt, metrics_only=args.metrics_only, 
                weighted=args.weighted, trendfile=args.trendfile)
