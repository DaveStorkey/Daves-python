#! /usr/bin/env python

'''
Routine to calculate topostrophy from model velocity fields. 

@author: Dave Storkey
@date: August 2018
'''

import iris
import numpy as np
import numpy.ma as ma
import netCDF4 as nc
import U_at_T_points as UaT

def bin_and_weight_field(field_in, weights=None, nbins=55, total_depths=None, cell_depths=None):

    bin_width = 5500.0/float(nbins)
    print 'nbins, bin_width : ',nbins, bin_width
    bins_low_bound = np.arange(nbins)*bin_width
    bins_high_bound = bins_low_bound+bin_width
    bins_centre = 0.5 * ( bins_low_bound + bins_high_bound )


    numerator = np.zeros((nbins,nbins))
    denominator = np.zeros((nbins,nbins))

    # Loop over bins for total depth of the water column:
    for count_total, (bin_low_total,bin_high_total) in enumerate(zip(bins_low_bound,bins_high_bound)):
        print 'count_total : ',count_total
        # Loop over bins for individual cell depth:
        for count_cell, (bin_low_cell,bin_high_cell) in enumerate(zip(bins_low_bound,bins_high_bound)):
            print 'count_cell : ',count_cell
            bin_indices = ma.where( ( (bin_low_total < total_depths) & (bin_high_total > total_depths) ) &
                                    ( (bin_low_cell  < cell_depths)  & (bin_high_cell  > cell_depths) ) )
#            if count_total == 14 and count_cell == 11:
#                field_out = ma.array(field_in,mask=False,copy=True)
#                field_out = ma.where( ( ( (bin_low_total < total_depths) & (bin_high_total > total_depths) ) &
#                                      ( (bin_low_cell  < cell_depths)  & (bin_high_cell  > cell_depths) ) ) &
#                                      ~field_in.mask, field_in, 0.0)
#                print 'min/max field_out : ',np.min(field_out),np.max(field_out)
#                field_out[np.where(np.isnan(field_out))] = 100.0
#                weights_out = ma.array(weights,mask=False,copy=True)
#                weights_out = ma.where( ( ( (bin_low_total < total_depths) & (bin_high_total > total_depths) ) &
#                                        ( (bin_low_cell  < cell_depths)  & (bin_high_cell  > cell_depths) ) ) &
#                                        ~weights.mask, weights, 0.0)
#            print 'min/max field_in[bin_indices] : ',ma.min(field_in[bin_indices]),ma.max(field_in[bin_indices])
#            print 'len(ma.where(np.isnan(field_in[bin_indices]))[0]) :',len(ma.where(np.isnan(field_in[bin_indices]))[0])
            numerator[count_cell,count_total] = ma.sum(field_in[bin_indices] * weights[bin_indices])
            denominator[count_cell,count_total] = ma.sum(weights[bin_indices])
            print 'numerator = ',numerator[count_cell,count_total]
            print 'denominator = ',denominator[count_cell,count_total]

    return bins_centre, numerator, denominator
   
def topostrophy(filenamestem=None, bathyfile=None, domcfg=None, utopofile=None, topofile=None, nbins=None):

    # 1. Read in bathymetry and horizontal grid cell dimensions
    #    and calculate field of (normalised) u_topo.

    with nc.Dataset(domcfg,'r') as domcfgfile:
        # get rid of size 1 time dimension as we read fields in.
        e1u = domcfgfile.variables['e1u'][0][:]
        e2v = domcfgfile.variables['e2v'][0][:]
        e1t = domcfgfile.variables['e1t'][0][:]
        e2t = domcfgfile.variables['e2t'][0][:]
        e3t_0 = domcfgfile.variables['e3t_0'][0][:]

    bathy = iris.load_cube(bathyfile,'sea_floor_depth')
    lats = bathy.coord('latitude').points
    ff = np.sin( lats[:] * np.pi / 180.0 )

    # dH/dx
    bathy_xm1 = np.roll(bathy.data[:,:], 1,axis=1)
    bathy_xp1 = np.roll(bathy.data[:,:],-1,axis=1)
    e1u_xm1 = np.roll(e1u, 1,axis=1)
    dHdx = (bathy_xp1 - bathy_xm1)/(e1u_xm1 + e1u)

    # dH/dy
    bathy_ym1 = np.roll(bathy.data[:,:], 1,axis=0)
    bathy_yp1 = np.roll(bathy.data[:,:],-1,axis=0)
    e2v_ym1 = np.roll(e2v, 1,axis=0)
    dHdy = (bathy_yp1 - bathy_ym1)/(e2v_ym1 + e2v)

    abs_grad_H = np.sqrt(dHdx*dHdx + dHdy*dHdy)

    ffout = bathy.copy()
    ffout.var_name = 'ff'
    ffout.standard_name = 'coriolis_parameter'
    ffout.units = 's-1'
    ffout.data[:] = ff[:]

    absgradH = bathy.copy()
    absgradH.var_name = 'AbsGradH'
    absgradH.standard_name = 'sea_floor_depth'
    absgradH.units = 'm'
    absgradH.data[:] = abs_grad_H

    dHdxout = bathy.copy()
    dHdxout.var_name = 'dHdx'
    dHdxout.standard_name = 'sea_floor_depth'
    dHdxout.units = '1'
    dHdxout.data[:] = dHdx[:]

    dHdyout = bathy.copy()
    dHdyout.var_name = 'dHdy'
    dHdyout.standard_name = 'sea_floor_depth'
    dHdyout.units = '1'
    dHdyout.data[:] = dHdy[:]

    utopo = bathy.copy()
    utopo.var_name = 'uo'
    utopo.standard_name = 'sea_water_x_velocity'
    utopo.units = 'm/s'
    utopo.data[:] = 0.0
    vtopo = bathy.copy()
    vtopo.var_name = 'vo'
    vtopo.standard_name = 'sea_water_y_velocity'
    vtopo.units = 'm/s'
    vtopo.data[:] = 0.0

    utopo.data[:] = ma.where( np.abs(ff) * abs_grad_H != 0.0, (ff * dHdy)/(np.abs(ff) * abs_grad_H),0.0)
    vtopo.data[:] = ma.where( np.abs(ff) * abs_grad_H != 0.0, (-ff * dHdx)/(np.abs(ff) * abs_grad_H), 0.0)

    # 2. Read in model velocity field and calculate dot product
    #    and cross product of normalised model velocity field with
    #    u_topo.

    u_model = iris.load_cube(filenamestem+'_grid-U.nc','sea_water_x_velocity')[0]
    v_model = iris.load_cube(filenamestem+'_grid-V.nc','sea_water_y_velocity')[0]

    uT,vT = UaT.U_at_T_points(u_model,v_model)

    # normalise:
    uv_mag = np.sqrt(uT.data*uT.data + vT.data*vT.data)
    uT.data = uT.data/uv_mag
    vT.data = vT.data/uv_mag

    # form the dot product and cross product with u_topo:

    u_dot = uT.copy()
    u_dot.data[:] = 0.0    

    u_dot.data[:] = uT.data[:]*utopo.data[:] + vT.data[:]*vtopo.data[:]

    u_cross = uT.copy()
    u_cross.data[:] = 0.0    

    u_cross.data[:] = uT.data[:]*vtopo.data[:] - vT.data[:]*utopo.data[:]

    iris.FUTURE.netcdf_no_unlimited = True
    iris.save([ffout,dHdxout,dHdyout,absgradH,utopo,vtopo,u_dot,u_cross],utopofile)

    # 3. Topostrophy calculation:

    cell_depths_T = np.zeros(e3t_0.shape)
    cell_depths_T[0] = e3t_0[0] 
    for level in np.arange(cell_depths_T.shape[0]-1) + 1:
        cell_depths_T[level] = cell_depths_T[level-1] + e3t_0[level]
    bathy3D = np.ones(e3t_0.shape)
    bathy3D = bathy.data * bathy3D
    cell_volume_T = ma.array(e1t[:]*e2t[:]*e3t_0[:],mask=u_dot.data.mask,copy=True)
    weights = abs_grad_H * np.abs(ff)

    (nz,ny,nx) = u_dot.shape
    depths,numerator,denominator = bin_and_weight_field(u_dot.data[:nz,1:ny,1:nx],
                                                 weights=weights[1:ny,1:nx]*cell_volume_T[:nz,1:ny,1:nx],
                                                 total_depths=bathy3D[:nz,1:ny,1:nx],
                                                 cell_depths=cell_depths_T[:nz,1:ny,1:nx], nbins=nbins)

#    with nc.Dataset('debug_file.nc','w') as debug_out:
#        debug_out.createDimension('x',bathy3D.shape[-1])
#        debug_out.createDimension('y',bathy3D.shape[-2])
#        debug_out.createDimension('z',bathy3D.shape[-3])
#        debug_out.createVariable('bathy3D',datatype='f',dimensions=('z','y','x'))
#        debug_out.variables['bathy3D'][:] = bathy3D[:]
#        debug_out.createVariable('cell_depths_T',datatype='f',dimensions=('z','y','x'))
#        debug_out.variables['cell_depths_T'][:] = cell_depths_T[:]
#        debug_out.createVariable('weights',datatype='f',dimensions=('y','x'))
#        debug_out.variables['weights'][:] = weights[:]
#        debug_out.createVariable('field_out',datatype='f',dimensions=('z','y','x'))
#        debug_out.variables['field_out'][:nz,1:ny,1:nx] = field_out[:]
#        debug_out.createVariable('weights_out',datatype='f',dimensions=('z','y','x'))
#        debug_out.variables['weights_out'][:nz,1:ny,1:nx] = weights_out[:]

    topostrophy = ma.where(denominator != 0.0, numerator/denominator, 0.0)

    with nc.Dataset(topofile,'w') as topo_out:
        topo_out.createDimension('x',topostrophy.shape[-1])
        topo_out.createDimension('y',topostrophy.shape[-2])
        topo_out.createVariable('depthx',datatype='f',dimensions=('x'))
        topo_out.variables['depthx'][:] = depths[:]
        topo_out.variables['depthx'].units = 'm'
        topo_out.variables['depthx'].standard_name = 'depth_below_geoid'
        topo_out.createVariable('depthy',datatype='f',dimensions=('y'))
        topo_out.variables['depthy'][:] = depths[:]
        topo_out.variables['depthy'].units = 'm'
        topo_out.variables['depthy'].standard_name = 'depth_below_geoid'
        topo_out.createVariable('numerator',datatype='f',dimensions=('y','x'))
        topo_out.variables['numerator'][:] = numerator[:]
        topo_out.variables['numerator'].coordinates = 'depthy depthx'
        topo_out.createVariable('denominator',datatype='f',dimensions=('y','x'))
        topo_out.variables['denominator'][:] = denominator[:]
        topo_out.variables['denominator'].coordinates = 'depthy depthx'
        topo_out.createVariable('topostrophy',datatype='f',dimensions=('y','x'))
        topo_out.variables['topostrophy'][:] = topostrophy[:]
        topo_out.variables['topostrophy'].coordinates = 'depthy depthx'

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("filenamestem", help="filename stem of input files (grid T, U, V)")
    parser.add_argument("-b", "--bathyfile", action="store",dest="bathyfile",default=None,
                   help="file containing model bathymetry in metres")
    parser.add_argument("-d", "--domcfg", action="store",dest="domcfg",default=None,
                   help="domain_cfg file")
    parser.add_argument("-u", "--utopofile", action="store",dest="utopofile",default=None,
                   help="output file for u_topo")
    parser.add_argument("-t", "--topofile", action="store",dest="topofile",default=None,
                   help="output file for topostrophy")
    parser.add_argument("-n", "--nbins", action="store",dest="nbins",type=int,default=55,
                   help="number of depth bins for topostrophy calculation")

    args = parser.parse_args()

    topostrophy(args.filenamestem, bathyfile=args.bathyfile, domcfg=args.domcfg, utopofile=args.utopofile,
                topofile=args.topofile, nbins=args.nbins)
