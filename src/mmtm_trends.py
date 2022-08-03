#! /usr/bin/env python
'''
Routine to do sanity checks of momentum trend diagnostics
from NEMO.

Hardwired to expect output from a run with vector invariant
momentum advection for now.

Created August 2019

Apr 2021 : Add total trend diagnostic minus model trend residual as option. 

Aug 2021 : Generalise to do V component as well.

@author: Dave Storkey
'''

import netCDF4 as nc
import numpy as np
import numpy.ma as ma

def mmtm_trends(filename,components=None,xcpt=None,metrics_only=None,advsum=None,html=None,trendfile=None):

    if metrics_only:
        trends_file = nc.Dataset(filename,'r')
    else:
        trends_file = nc.Dataset(filename,'r+')

    for trdname in ['utrd','vtrd']:
        try:
            total_trend = trends_file.variables[trdname+'_tot']
        except KeyError:
            pass
        else:
            break
    else:
        raise Exception('Error : could not find total trend diagnostic utrd_tot or vtrd_tot in file '+filename)

    if metrics_only:

        utrd={}
        for thickname in ['e3u','e3v','thkcello']:
            try:
                e3u = trends_file.variables[thickname][:]
            except KeyError:
                pass
            else:
                break
        else:
            raise Exception("Error: couldn't find cell thickness variable.")
        utrd['tot'] = trends_file.variables[trdname+'_tot'][:]
        utrd['res'] = trends_file.variables[trdname+'_res'][:]

    else:

        # 'atf' trend not included in standard list because 'tot' trend in NEMO doesn't include this.
        cpt_list = ['hpg', 'spg', 'spg2d', 'keg', 'rvo', 'pvo', 'pvo2d', 'zad', 'ldf', 'zdf', 'bta', 'bfr', 'bfr2d']

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
        for thickname in ['e3u','e3v','thkcello']:
            try:
                e3u = trends_file.variables[thickname][:]
            except KeyError:
                pass
            else:
                break
        else:
            raise Exception("Error: couldn't find cell thickness variable.")
        utrd['tot'] = trends_file.variables[trdname+'_tot']

        utrd['sum'] = ma.array(utrd['tot'][:]) 
        utrd['sum'][:] = 0.0
        utrd['res'] = ma.array(utrd['tot'][:]) 
        utrd['res'][:] = 0.0
        umask = ma.ones(utrd['tot'][:].shape[-3:])
        print(trdname+'[tot][:].shape[-3:] : ',utrd['tot'][:].shape[-3:])
        print('umask.shape : ',umask.shape)
        umask[:] = 1.0
        umask[-1,:,:] = 0.0
        for cpt in components:
            utrd[cpt] = trends_file.variables[trdname+'_'+cpt]
            # note with numpy broadcasting rules this should add a 
            # 2D cpt array to each level of the 3D sum array. 
            utrd['sum'][:] = utrd['sum'][:] + utrd[cpt][:] 

        # mask out bottom level
        utrd['sum'][:] = utrd['sum'][:] * umask[:]

        if advsum:
            utrd['adv'] = ma.array(utrd['tot'][:]) 
            utrd['adv'][:] = utrd['rvo'][:] + utrd['keg'][:] + utrd['zad'][:]

        utrd['res'][:] = utrd['tot'][:] - utrd['sum'][:]
    
        if trendfile is not None:
            trendfiledata = nc.Dataset(trendfile,'r')
            utrd['model'] = ma.array(utrd['tot'][:])
            utrd['res2']  = ma.array(utrd['tot'][:])
            utrd['model'] = trendfiledata.variables['vozocrtx']
            utrd['res2'][:] = utrd['tot'][:] - utrd['model'][:]

    if len(utrd['res'].shape) == 3:
        axes = (1,2)
        zlength = utrd['res'].shape[0]
    elif len(utrd['res'].shape) == 4:
        axes = (0,2,3)
        zlength = utrd['res'].shape[1]

    utrd_res_metric = ma.average(ma.absolute(utrd['res'][:]),axis=axes,weights=e3u[:])
    utrd_tot_metric = ma.average(ma.absolute(utrd['tot'][:]),axis=axes,weights=e3u[:])
    utrd_res_metric_fracerr = utrd_res_metric/utrd_tot_metric
    per_level_weights = ma.average(e3u[:],axis=axes)
    if html:
        print('<center>')
        print('<p><b>volume average absolute residual : ', 
               ma.average(utrd_res_metric[:],weights=per_level_weights),' m2/s2<br>')
        print('<p><b>volume average absolute total trend : ', 
               ma.average(utrd_tot_metric[:],weights=per_level_weights),' m2/s2<br>')
        print('as fractional error compared to volume average absolute total trend : ',
               ma.average(utrd_res_metric_fracerr[:],weights=per_level_weights),'<br><br>')
        print('Following table gives the same quantities by model level:')
        print('</b></p> ')
        print('<table BORDER=1 COLS=3 WIDTH="30%" NOSAVE >')
        print('<tr>')
        print('<th><b>level</b></th><th><b>residual (m2/s2)</b></th><th><b>fractional err</b></th>')
        for z in range(zlength):
            print('</tr><tr>')
            print('<td>',z,'</td><td>',utrd_res_metric[z],'</td><td>',utrd_tot_metric[z],'</td><td>',utrd_res_metric_fracerr[z],'</td>')
        print('</tr>')
        print('</table></center>')

    else:
        print('sum(abs(residual thickness-weighted average of components minus total trend diagnostic)) : ',
               ma.average(utrd_res_metric[:],weights=per_level_weights) )
        print('sum(abs(total trend diagnostic)) : ',
               ma.average(utrd_tot_metric[:],weights=per_level_weights) )
        print('as fractional error compared to weighted average(abs(utrd_tot)) : ',
               ma.average(utrd_res_metric_fracerr[:],weights=per_level_weights) )
        print(' ')
        print('level, residual, fractional err : ')
        for z in range(zlength):
            print(z,':',utrd_res_metric[z],utrd_tot_metric[z],utrd_res_metric_fracerr[z])

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
        if advsum:
            output_list = output_list+['adv']
        if trendfile is not None:
            output_list = output_list+['model','res2']

        for cpt in output_list:
            if not trdname+'_'+cpt in trends_file.variables.keys():
                if fillvalue:
                    trends_file.createVariable(trdname+'_'+cpt,datatype='f',dimensions=utrd['tot'].dimensions,fill_value=utrd['tot']._FillValue)
                else:
                    trends_file.createVariable(trdname+'_'+cpt,datatype='f',dimensions=utrd['tot'].dimensions)
                utrd_out = trends_file.variables[trdname+'_'+cpt]
                utrd_out.setncatts(attrdict)
                if cpt == 'sum':
                    utrd_out.long_name='Sum of U trends: '+' '.join(components)
                elif cpt == 'adv':
                    utrd_out.long_name='i-trend: total advection = rvo + keg + zad'
                elif cpt == 'model':
                    utrd_out.long_name='i-trend: model end minus model start'
                elif cpt == 'res2':
                    utrd_out.long_name='Total U trend minus Model U trend'
                else:
                    utrd_out.long_name='Total U trend minus sum of U trends'
            else:
                utrd_out = trends_file.variables[trdname+'_'+cpt]

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
    parser.add_argument("-a","--advsum",action="store_true",dest="advsum",
                    help="sum advective trends to give total advective trend")
    parser.add_argument("-H","--HTML",action="store_true",dest="html",
                    help="HTML formatted output")
    parser.add_argument("-T","--trendfile",action="store",dest="trendfile",
                    help="name of file with model trends")

    args = parser.parse_args()

    mmtm_trends(args.filename, components=args.components, xcpt=args.xcpt, metrics_only=args.metrics_only, 
                advsum=args.advsum, html=args.html, trendfile=args.trendfile)
