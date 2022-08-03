#! /usr/bin/env python
"""
    Script to generate an extended ORCA grid to include coverage of the
    Ross and Ronne-Filchner ice shelves. 
"""

__author__ = 'Dave Storkey (dave.storkey@metoffice.gov.uk)'
__date__ = '$Date: December 2013 $'

import sys
import math
import numpy as np
import netCDF4

####################################################################################

def write_coords(fileout,datatype,glamt=None,glamu=None,glamv=None,glamf=None,
                                  gphit=None,gphiu=None,gphiv=None,gphif=None,
                                  e1t =None,e1u =None,e1v =None,e1f =None,
                                  e2t =None,e2u =None,e2v =None,e2f =None,
                                  aniso =False, ortho=None ):

    dataset_out = netCDF4.Dataset(fileout, mode='w')

    for var in [glamt, glamu, glamv, glamf,
                gphit, gphiu, gphiv, gphif, 
                e1t , e1u , e1v , e1f, 
                e2t , e2u , e2v , e2f, ortho ]: 
        try:
            dimy = var.shape[0]
        except AttributeError:
            pass
        else:
            dimx = var.shape[1]
            break
    else:
        raise AttributeError

    dataset_out.createDimension('y',dimy)     
    dataset_out.createDimension('x',dimx)     

    for [var,nvar] in zip( [ glamt , glamu , glamv , glamf , 
                             gphit , gphiu , gphiv , gphif , 
                             e1t  , e1u  , e1v  , e1f  ,  
                             e2t  , e2u  , e2v  , e2f, ortho ], 
                           ['glamt','glamu','glamv','glamf',
                            'gphit','gphiu','gphiv','gphif', 
                            'e1t','e1u','e1v','e1f', 
                            'e2t','e2u','e2v','e2f', 'ortho' ] ):
        if var is not None:
            dataset_out.createVariable(nvar,datatype,dimensions=('y','x'))
            var_out = dataset_out.variables[nvar]
            var_out[:,:] = var[:,:]

    if aniso and e1t is not None and e2t is not None:
        dataset_out.createVariable('anisotropy',datatype,dimensions=('y','x'))
        var_out = dataset_out.variables['anisotropy']
        var_out[:,:] = e2t[:,:]/(np.maximum(e1t[:,:],0.0001))

    dataset_out.close()

####################################################################################

def fix_range(lon_in, lon_range):

    lon_out = np.where( lon_in <= lon_range[0], lon_in + 360.0, lon_in )    
    lon_out = np.where( lon_out >= lon_range[1], lon_out - 360.0, lon_out )    

    return lon_out

####################################################################################

def e1_calc(lam, phi, jpfac):

# radius of the earth:
    ra = 6371229.0
    deg2rad = math.pi/180.0

    nj = lam.shape[0]
    ni = lam.shape[2]
    e1 = np.zeros([nj,ni],dtype=np.float64)

    cosp = np.cos( phi[:,0,:,0] * deg2rad )
    dldi =     lam[:,0,:,-2] - 8.0*lam[:,0,:,-1] \
         + 8.0*lam[:,0,:,+1] -     lam[:,0,:,+2] 
# grid spacing in i-space is 1/jpfac
    dldi = deg2rad * dldi * jpfac/12.0
    dpdi =     phi[:,0,:,-2] - 8.0*phi[:,0,:,-1] \
         + 8.0*phi[:,0,:,+1] -     phi[:,0,:,+2] 
# grid spacing in i-space is 1/jpfac
    dpdi = deg2rad * dpdi * jpfac/12.0

    e1 = ra * np.sqrt( dldi*dldi*cosp*cosp + dpdi*dpdi )        

    return e1

####################################################################################

def e2_calc(lam, phi, jpfac):

# radius of the earth:
    ra = 6371229.0
    deg2rad = math.pi/180.0

    nj = lam.shape[0]
    ni = lam.shape[2]
    e2 = np.zeros([nj,ni])

    cosp = np.cos( phi[:,0,:,0] * deg2rad )

    dldj =     lam[:,-2,:,0] - 8.0*lam[:,-1,:,0] \
         + 8.0*lam[:,+1,:,0] -     lam[:,+2,:,0] 
# grid spacing in j-space is 1/jpfac
    dldj = deg2rad * dldj * jpfac/12.0

    dpdj =     phi[:,-2,:,0] - 8.0*phi[:,-1,:,0] \
         + 8.0*phi[:,+1,:,0] -     phi[:,+2,:,0] 
# grid spacing in j-space is 1/jpfac
    dpdj = deg2rad * dpdj * jpfac/12.0

    e2 = ra * np.sqrt( dldj*dldj*cosp*cosp + dpdj*dpdj )        

#    print 'phi[0,0,0,0], phi[1,0,0,0], phi[2,0,0,0] :',phi[0,0,0,0], phi[1,0,0,0], phi[2,0,0,0]
#    print 'cosp[0,0] :',cosp[0,0]
#    print 'dldj[0,0]/deg2rad :',dldj[0,0]/deg2rad
#    print 'dpdj[0,0]/deg2rad :',dpdj[0,0]/deg2rad
#    print 'e2[0,0] :',e2[0,0]

    return e2

####################################################################################

def orthogonality_check(lam, phi, e1, e2, jpfac):

# radius of the earth:
    ra = 6371229.0
    deg2rad = math.pi/180.0
    rad2deg = 180.0/math.pi

    nj = lam.shape[0]
    ni = lam.shape[2]
    ortho = np.zeros([nj,ni])

    sinp = np.sin( phi[:,0,:,0] * deg2rad )
    cosp = np.cos( phi[:,0,:,0] * deg2rad )

    dldi =     lam[:,0,:,-2] - 8.0*lam[:,0,:,-1] \
         + 8.0*lam[:,0,:,+1] -     lam[:,0,:,+2] 
# grid spacing in i-space is 1/jpfac
    dldi = deg2rad * dldi * jpfac/12.0

    dpdi =     phi[:,0,:,-2] - 8.0*phi[:,0,:,-1] \
         + 8.0*phi[:,0,:,+1] -     phi[:,0,:,+2] 
# grid spacing in i-space is 1/jpfac
    dpdi = deg2rad * dpdi * jpfac/12.0

    dldj =     lam[:,-2,:,0] - 8.0*lam[:,-1,:,0] \
         + 8.0*lam[:,+1,:,0] -     lam[:,+2,:,0] 
# grid spacing in j-space is 1/jpfac
    dldj = deg2rad * dldj * jpfac/12.0

    dpdj =     phi[:,-2,:,0] - 8.0*phi[:,-1,:,0] \
         + 8.0*phi[:,+1,:,0] -     phi[:,+2,:,0] 
# grid spacing in j-space is 1/jpfac
    dpdj = deg2rad * dpdj * jpfac/12.0

    ortho = rad2deg * np.arccos( (ra*ra) * (cosp*cosp*dldi*dldj + sinp*sinp*dpdi*dpdj) / (e1[:,:]*e2[:,:]) )

    return ortho

####################################################################################
####################################################################################

def extend_orca(base_coords,north_coords,lat_cut=None,jpfac=10.0):

############################################################
# 1. Specify bounds of Ross and Ronne-Filchner Ice Shelves
#    allowing for some land points at the edges. 
############################################################

# Specify longitudes in the range [0,360] (because the zero meridion
# doesn't bisect either of our ice shelves but the 180 degree
# meridion does bisect the Ross Ice Shelf). 

    lam_west   = { 'ross':  150.0  , 'ronne-filchner': 270.0}
    lam_centre = { 'ross':  185.0  , 'ronne-filchner': 302.5}
    lam_east   = { 'ross':  220.0  , 'ronne-filchner': 335.0}
    phi_south  = { 'ross':  -86.0  , 'ronne-filchner': -84.0}

###################################################################
# 2. Read in base coordinates and coordinates for northern segment
###################################################################

# Base Grid (to be extended):

    base_dataset = netCDF4.Dataset(base_coords, mode='r')

# This is a netcdf variable object (useful for getting things like dimensions):
    glamt_base = base_dataset.variables['glamt']

# These are simple numpy arrays (useful for numpy processing):
    lamt_base = fix_range( np.squeeze(np.copy(base_dataset.variables['glamt'])), [0.0,360.0] )
    lamu_base = fix_range( np.squeeze(np.copy(base_dataset.variables['glamu'])), [0.0,360.0] )
    lamv_base = fix_range( np.squeeze(np.copy(base_dataset.variables['glamv'])), [0.0,360.0] )
    lamf_base = fix_range( np.squeeze(np.copy(base_dataset.variables['glamf'])), [0.0,360.0] )

    phit_base = np.squeeze(np.copy(base_dataset.variables['gphit']))
    phiu_base = np.squeeze(np.copy(base_dataset.variables['gphiu']))
    phiv_base = np.squeeze(np.copy(base_dataset.variables['gphiv']))
    phif_base = np.squeeze(np.copy(base_dataset.variables['gphif']))

    e1t_base = np.squeeze(np.copy(base_dataset.variables['e1t']))
    e1u_base = np.squeeze(np.copy(base_dataset.variables['e1u']))
    e1v_base = np.squeeze(np.copy(base_dataset.variables['e1v']))
    e1f_base = np.squeeze(np.copy(base_dataset.variables['e1f']))

    e2t_base = np.squeeze(np.copy(base_dataset.variables['e2t']))
    e2u_base = np.squeeze(np.copy(base_dataset.variables['e2u']))
    e2v_base = np.squeeze(np.copy(base_dataset.variables['e2v']))
    e2f_base = np.squeeze(np.copy(base_dataset.variables['e2f']))

    lat_south = phit_base[0,0]

# North Grid (the grid from which we extract the northern segment):

    north_dataset = netCDF4.Dataset(north_coords, mode='r')

# These are netcdf variable objects (useful for getting things like dimensions):
    glamt_north_grid = north_dataset.variables['glamt']
    glamu_north_grid = north_dataset.variables['glamu']
    glamv_north_grid = north_dataset.variables['glamv']
    glamf_north_grid = north_dataset.variables['glamf']

# These are simple numpy arrays with dimensions of size 1 removed 
# (useful for numpy processing):
    lamt_north_grid = fix_range( np.squeeze(np.copy(north_dataset.variables['glamt'])), [0.0,360.0] )
    lamu_north_grid = fix_range( np.squeeze(np.copy(north_dataset.variables['glamu'])), [0.0,360.0] )
    lamv_north_grid = fix_range( np.squeeze(np.copy(north_dataset.variables['glamv'])), [0.0,360.0] )
    lamf_north_grid = fix_range( np.squeeze(np.copy(north_dataset.variables['glamf'])), [0.0,360.0] )

    phit_north_grid = np.squeeze(np.copy(north_dataset.variables['gphit']))
    phiu_north_grid = np.squeeze(np.copy(north_dataset.variables['gphiu']))
    phiv_north_grid = np.squeeze(np.copy(north_dataset.variables['gphiv']))
    phif_north_grid = np.squeeze(np.copy(north_dataset.variables['gphif']))

###################################################################
# 3. Extract and scale the northern segment so that it fits into 
#    the southern hole in the base grid. 
###################################################################

# Number of extra "halo" points to allow at the north and south edges
# of the extracted grid so we can do interpolation later:

    n_halo = 3

# Extract the distorted part of the ORCA grid (north of lat_cut):

    if lat_cut is not None:
        ind_cut = np.where(phit_north_grid >= lat_cut)
        j_cut = ind_cut[0].min()
        print 'j_cut : ',j_cut

# Include n_halo points south of lat_cut to facilitate 4th-order 
# interpolation in the j-direction later. So the line of points
# which will be matched to the base grid southern line is
# the (n_halo_1)th line of the T-U arrays. 

    lamt_north_circle = lamt_north_grid[j_cut-n_halo-1:-1,:]
    lamu_north_circle = lamu_north_grid[j_cut-n_halo-1:-1,:]
    lamv_north_circle = lamv_north_grid[j_cut-n_halo-1:-1,:]
    lamf_north_circle = lamf_north_grid[j_cut-n_halo-1:-1,:]

    phit_north_circle = phit_north_grid[j_cut-n_halo-1:-1,:]
    phiu_north_circle = phiu_north_grid[j_cut-n_halo-1:-1,:]
    phiv_north_circle = phiv_north_grid[j_cut-n_halo-1:-1,:]
    phif_north_circle = phif_north_grid[j_cut-n_halo-1:-1,:]

    lat_cut_exact = phit_north_circle[n_halo,0]

# Scale the latitudes and reverse the sign so that the extracted part of the grid 
# sits within the southern hole in the base grid. Also need to scale the scale 
# factors by the same amount. 

    scale_factor = (90.0+lat_south)/(90.0-lat_cut_exact)
    print 'scale_factor : ',scale_factor

    lamt_south_circle = lamt_north_circle
    lamu_south_circle = lamu_north_circle
    lamv_south_circle = lamv_north_circle
    lamf_south_circle = lamf_north_circle

    phit_south_circle = ( (90.0-phit_north_circle)*scale_factor ) - 90.0
    phiu_south_circle = ( (90.0-phiu_north_circle)*scale_factor ) - 90.0
    phiv_south_circle = ( (90.0-phiv_north_circle)*scale_factor ) - 90.0
    phif_south_circle = ( (90.0-phif_north_circle)*scale_factor ) - 90.0

# Find the i-value of the north pole. (Actually there are 2 i-values of course - use
# the first one). Find the latitude that this corresponds to for the northern edge
# of this piece of grid. 

    ind_south = np.unravel_index(phit_south_circle.argmin(),phit_south_circle.shape)
    lam_centre_south_circle = lamt_south_circle[n_halo,ind_south[1]]

    print 'ind_south : ',ind_south
    print 'phit_south_circle[ind_south] : ',phit_south_circle[ind_south]
    print 'lam_centre_south_circle : ',lam_centre_south_circle

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# LOOP ON ICE SHELVES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Initialise dictionaries of grid segments.

    ii_seg = {}
    ni_seg = {}
    nj_seg = {}

    lamt_seg = {}
    lamu_seg = {}
    lamv_seg = {}
    lamf_seg = {}

    phit_seg = {}
    phiu_seg = {}
    phiv_seg = {}
    phif_seg = {}

    e1t_seg = {}
    e1u_seg = {}
    e1v_seg = {}
    e1f_seg = {}

    e2t_seg = {}
    e2u_seg = {}
    e2v_seg = {}
    e2f_seg = {}

# If you just loop on lam_west.keys() here and you only have one
# dictionary entry it will loop on the individual characters of 
# of the key string...

    for ii_shelf in range(len(lam_west)):

        iceshelf = lam_west.keys()[ii_shelf]

        print '==================================================='
        print '>>> Processing grid for',iceshelf,'ice shelf'
        print '==================================================='

###################################################################
# 4. Rotate the southern segment so that we have the isotropic
#    part over our ice shelf. 
###################################################################

        rotation_angle = lam_centre[iceshelf] - lam_centre_south_circle
        print 'rotation_angle : ',rotation_angle

        phit_rotated = phit_south_circle
        phiu_rotated = phiu_south_circle
        phiv_rotated = phiv_south_circle
        phif_rotated = phif_south_circle

        lamt_rotated = fix_range( lamt_south_circle + rotation_angle, [0.0,360.0] )
        lamu_rotated = fix_range( lamu_south_circle + rotation_angle, [0.0,360.0] )
        lamv_rotated = fix_range( lamv_south_circle + rotation_angle, [0.0,360.0] )
        lamf_rotated = fix_range( lamf_south_circle + rotation_angle, [0.0,360.0] )

#        write_coords('coords_'+iceshelf+'_rotated.nc',north_dataset,glamt_north_grid,lamt_rotated,phit_rotated)

#############################################################################
# 5. Create an intermediate grid segment ("seg1") for this ice shelf to hold 
#    a grid interpolated in the i-direction but not the j-direction.
#############################################################################

        lamt_source = np.squeeze( lamt_rotated[n_halo] )
        lamu_source = np.squeeze( lamu_rotated[n_halo] )

# For each point on the target grid we require extra points on a finer grid
# to calculate the e1 and e2 scale factors via finite differencing. Use
# 4th-order differencing so require 4 extra points or 5 points in total:

        npoints = 5

# Determine dimensions of intermediate grid segment

# "&" is a bitwise "and" operator in the following expression:
        ii_seg[iceshelf] = np.where( ( lamt_base[0] >= lam_west[iceshelf] ) & 
                                     ( lamt_base[0] <= lam_east[iceshelf] ) )

        lamt_target = np.squeeze( lamt_base[0,ii_seg[iceshelf]] )
        lamu_target = np.squeeze( lamu_base[0,ii_seg[iceshelf]] )

        dlam_target_fine = (lamt_target[1]-lamt_target[0])/jpfac
        print 'dlam_target_fine :',dlam_target_fine

        ni_seg1 = len(lamt_target)       

# "&" is a bitwise "and" operator in the following expression:
        nj_seg1 = np.where( ( phit_rotated >= phi_south[iceshelf] ) &  
                            ( lamt_rotated >= lam_west[iceshelf] )  &
                            ( lamt_rotated <= lam_east[iceshelf] ) )[0].max() + 1 
# Allow n_halo extra grid points at the southern edge so we can do 4th-order
# interpolation in the j-direction later. Already have n_halo extra points
# at the northern edge.
        nj_seg1 = nj_seg1 + n_halo

        print 'ni_seg1 : ',ni_seg1
        print 'nj_seg1 : ',nj_seg1

        lamt_seg1 = np.zeros([nj_seg1,ni_seg1,npoints])
        lamu_seg1 = np.zeros([nj_seg1,ni_seg1,npoints])
        lamv_seg1 = np.zeros([nj_seg1,ni_seg1,npoints])
        lamf_seg1 = np.zeros([nj_seg1,ni_seg1,npoints])

        phit_seg1 = np.zeros([nj_seg1,ni_seg1,npoints])
        phiu_seg1 = np.zeros([nj_seg1,ni_seg1,npoints])
        phiv_seg1 = np.zeros([nj_seg1,ni_seg1,npoints])
        phif_seg1 = np.zeros([nj_seg1,ni_seg1,npoints])

        print 'lamt_rotated.shape : ',lamt_rotated.shape
        print 'lamt_source.shape : ',lamt_source.shape
        print 'lamu_rotated.shape : ',lamu_rotated.shape
        print 'lamu_source.shape : ',lamu_source.shape

###################################################################
# 6. Work out the source grid points and interpolation coefficients 
#    for the interpolation in the i-direction along the first j-line
#    (at the join of the two grids). 
###################################################################

# For each target grid point (T/V or U/F), 4th order interpolation 
# in the i-direction will always use 2 source T/V points and 
# 2 source U/F points. 

#----------------------------------------------
# indices of source points (counting from zero)
#----------------------------------------------
        itv1 = np.zeros([ni_seg1,npoints,2],dtype=np.uint16)
        itv2 = np.zeros([ni_seg1,npoints,2],dtype=np.uint16)
        iuf1 = np.zeros([ni_seg1,npoints,2],dtype=np.uint16)
        iuf2 = np.zeros([ni_seg1,npoints,2],dtype=np.uint16)

        for ii in range(ni_seg1):
            for ip in np.arange(npoints)-2:
                lamt_tg = lamt_target[ii] + ip * dlam_target_fine
                lamu_tg = lamu_target[ii] + ip * dlam_target_fine
                
# T/V target points:
#------------------
# nearest T/V source point:
                itv1[ii,ip,0] = np.argmin(abs(lamt_tg - lamt_source))    
# next-nearest T/V source point:
                itv2[ii,ip,0] = itv1[ii,ip,0] + math.copysign(1, lamt_tg-lamt_source[itv1[ii,ip,0]])             
# nearest U/F source point:
                iuf1[ii,ip,0] = np.argmin(abs(lamt_tg - lamu_source))    
# next-nearest U/F source point:
                iuf2[ii,ip,0] = iuf1[ii,ip,0] + math.copysign(1, lamt_tg-lamu_source[iuf1[ii,ip,0]])             

# U/F target points:
#------------------
# nearest T/V source point:
                itv1[ii,ip,1] = np.argmin(abs(lamu_tg - lamt_source))    
# next-nearest T/V source point:
                itv2[ii,ip,1] = itv1[ii,ip,1] + math.copysign(1, lamu_tg-lamt_source[itv1[ii,ip,1]])             
# nearest U/F source point:
                iuf1[ii,ip,1] = np.argmin(abs(lamu_tg - lamu_source))    
# next-nearest U/F source point:
                iuf2[ii,ip,1] = iuf1[ii,ip,1] + math.copysign(1, lamu_tg-lamu_source[iuf1[ii,ip,1]])             

        print 'min/max itv1 : ',itv1.min(),itv1.max()
        print 'min/max itv2 : ',itv2.min(),itv2.max()
        print 'min/max iuf1 : ',iuf1.min(),iuf1.max()
        print 'min/max iuf2 : ',iuf2.min(),iuf2.max()

#-----------------------------
# interpolation coefficients
#-----------------------------
        ctv1 = np.zeros([ni_seg1,npoints,2])
        ctv2 = np.zeros([ni_seg1,npoints,2])
        cuf1 = np.zeros([ni_seg1,npoints,2])
        cuf2 = np.zeros([ni_seg1,npoints,2])

        for ii in range(ni_seg1):
            for ip in np.arange(npoints)-2:

# T/V target points:
#------------------
                lt_tgt = lamt_target[ii] + ip * dlam_target_fine
                lt_src1 = lamt_source[itv1[ii,ip,0]]
                lt_src2 = lamt_source[itv2[ii,ip,0]]
                lu_src1 = lamu_source[iuf1[ii,ip,0]]
                lu_src2 = lamu_source[iuf2[ii,ip,0]]

                ctv1[ii,ip,0] = ( (lt_tgt -lt_src2)*(lt_tgt -lu_src1)*(lt_tgt -lu_src2) ) \
                              / ( (lt_src1-lt_src2)*(lt_src1-lu_src1)*(lt_src1-lu_src2) )
                ctv2[ii,ip,0] = ( (lt_tgt -lt_src1)*(lt_tgt -lu_src1)*(lt_tgt -lu_src2) ) \
                              / ( (lt_src2-lt_src1)*(lt_src2-lu_src1)*(lt_src2-lu_src2) )
                cuf1[ii,ip,0] = ( (lt_tgt -lt_src1)*(lt_tgt -lt_src2)*(lt_tgt -lu_src2) ) \
                              / ( (lu_src1-lt_src1)*(lu_src1-lt_src2)*(lu_src1-lu_src2) )
                cuf2[ii,ip,0] = ( (lt_tgt -lt_src1)*(lt_tgt -lt_src2)*(lt_tgt -lu_src1) ) \
                              / ( (lu_src2-lt_src1)*(lu_src2-lt_src2)*(lu_src2-lu_src1) )

# U/F target points:
#------------------
                lu_tgt = lamu_target[ii] + ip * dlam_target_fine
                lt_src1 = lamt_source[itv1[ii,ip,1]]
                lt_src2 = lamt_source[itv2[ii,ip,1]]
                lu_src1 = lamu_source[iuf1[ii,ip,1]]
                lu_src2 = lamu_source[iuf2[ii,ip,1]]

                ctv1[ii,ip,1] = ( (lu_tgt -lt_src2)*(lu_tgt -lu_src1)*(lu_tgt -lu_src2) ) \
                              / ( (lt_src1-lt_src2)*(lt_src1-lu_src1)*(lt_src1-lu_src2) )
                ctv2[ii,ip,1] = ( (lu_tgt -lt_src1)*(lu_tgt -lu_src1)*(lu_tgt -lu_src2) ) \
                              / ( (lt_src2-lt_src1)*(lt_src2-lu_src1)*(lt_src2-lu_src2) )
                cuf1[ii,ip,1] = ( (lu_tgt -lt_src1)*(lu_tgt -lt_src2)*(lu_tgt -lu_src2) ) \
                              / ( (lu_src1-lt_src1)*(lu_src1-lt_src2)*(lu_src1-lu_src2) )
                cuf2[ii,ip,1] = ( (lu_tgt -lt_src1)*(lu_tgt -lt_src2)*(lu_tgt -lu_src1) ) \
                              / ( (lu_src2-lt_src1)*(lu_src2-lt_src2)*(lu_src2-lu_src1) )

###################################################################
# 7. Use the indices and coefficients calculated in section 6 to 
#    perform interpolation in the i-direction for each T-U line 
#    and each V-F line on the "seg1" grid. 
###################################################################

        for ij in range(nj_seg1):
 
            lamt_seg1[ij,:,:] = ctv1[:,:,0] * lamt_rotated[ij,itv1[:,:,0]] \
                              + ctv2[:,:,0] * lamt_rotated[ij,itv2[:,:,0]] \
                              + cuf1[:,:,0] * lamu_rotated[ij,iuf1[:,:,0]] \
                              + cuf2[:,:,0] * lamu_rotated[ij,iuf2[:,:,0]] 

            lamu_seg1[ij,:,:] = ctv1[:,:,1] * lamt_rotated[ij,itv1[:,:,1]] \
                              + ctv2[:,:,1] * lamt_rotated[ij,itv2[:,:,1]] \
                              + cuf1[:,:,1] * lamu_rotated[ij,iuf1[:,:,1]] \
                              + cuf2[:,:,1] * lamu_rotated[ij,iuf2[:,:,1]] 

            lamv_seg1[ij,:,:] = ctv1[:,:,0] * lamv_rotated[ij,itv1[:,:,0]] \
                              + ctv2[:,:,0] * lamv_rotated[ij,itv2[:,:,0]] \
                              + cuf1[:,:,0] * lamf_rotated[ij,iuf1[:,:,0]] \
                              + cuf2[:,:,0] * lamf_rotated[ij,iuf2[:,:,0]] 

            lamf_seg1[ij,:,:] = ctv1[:,:,1] * lamv_rotated[ij,itv1[:,:,1]] \
                              + ctv2[:,:,1] * lamv_rotated[ij,itv2[:,:,1]] \
                              + cuf1[:,:,1] * lamf_rotated[ij,iuf1[:,:,1]] \
                              + cuf2[:,:,1] * lamf_rotated[ij,iuf2[:,:,1]] 

            phit_seg1[ij,:,:] = ctv1[:,:,0] * phit_rotated[ij,itv1[:,:,0]] \
                              + ctv2[:,:,0] * phit_rotated[ij,itv2[:,:,0]] \
                              + cuf1[:,:,0] * phiu_rotated[ij,iuf1[:,:,0]] \
                              + cuf2[:,:,0] * phiu_rotated[ij,iuf2[:,:,0]] 

            phiu_seg1[ij,:,:] = ctv1[:,:,1] * phit_rotated[ij,itv1[:,:,1]] \
                              + ctv2[:,:,1] * phit_rotated[ij,itv2[:,:,1]] \
                              + cuf1[:,:,1] * phiu_rotated[ij,iuf1[:,:,1]] \
                              + cuf2[:,:,1] * phiu_rotated[ij,iuf2[:,:,1]] 

            phiv_seg1[ij,:,:] = ctv1[:,:,0] * phiv_rotated[ij,itv1[:,:,0]] \
                              + ctv2[:,:,0] * phiv_rotated[ij,itv2[:,:,0]] \
                              + cuf1[:,:,0] * phif_rotated[ij,iuf1[:,:,0]] \
                              + cuf2[:,:,0] * phif_rotated[ij,iuf2[:,:,0]] 

            phif_seg1[ij,:,:] = ctv1[:,:,1] * phiv_rotated[ij,itv1[:,:,1]] \
                              + ctv2[:,:,1] * phiv_rotated[ij,itv2[:,:,1]] \
                              + cuf1[:,:,1] * phif_rotated[ij,iuf1[:,:,1]] \
                              + cuf2[:,:,1] * phif_rotated[ij,iuf2[:,:,1]] 

###################################################################
# 8. Create the arrays to hold the final grid segments.
#    Create the new grid in the j-direction based on the ratio
#    of the i-direction resolutions between the two grids along 
#    the join. 
###################################################################

        ni_seg[iceshelf] = ni_seg1
        resolution_ratio = (lamt_source[1]-lamt_source[0])/(lamt_target[1]-lamt_target[0])
# take off the n_halo extra grid points at each end of the source grid before calculated nj_seg...
        nj_seg[iceshelf] = int((nj_seg1-2*n_halo) * resolution_ratio) + 1

        print 'resolution_ratio : ',resolution_ratio
        print 'ni_seg[iceshelf] : ',ni_seg[iceshelf]
        print 'nj_seg[iceshelf] : ',nj_seg[iceshelf]

        lamt_seg[iceshelf] = np.zeros([nj_seg[iceshelf],npoints,ni_seg[iceshelf],npoints])
        lamu_seg[iceshelf] = np.zeros([nj_seg[iceshelf],npoints,ni_seg[iceshelf],npoints])
        lamv_seg[iceshelf] = np.zeros([nj_seg[iceshelf],npoints,ni_seg[iceshelf],npoints])
        lamf_seg[iceshelf] = np.zeros([nj_seg[iceshelf],npoints,ni_seg[iceshelf],npoints])

        phit_seg[iceshelf] = np.zeros([nj_seg[iceshelf],npoints,ni_seg[iceshelf],npoints])
        phiu_seg[iceshelf] = np.zeros([nj_seg[iceshelf],npoints,ni_seg[iceshelf],npoints])
        phiv_seg[iceshelf] = np.zeros([nj_seg[iceshelf],npoints,ni_seg[iceshelf],npoints])
        phif_seg[iceshelf] = np.zeros([nj_seg[iceshelf],npoints,ni_seg[iceshelf],npoints])

###################################################################
# 9. Work out the source grid points and interpolation coefficients 
#    for the interpolation in the j-direction.
###################################################################

# Source points in j-space:

        jt_source = np.arange(nj_seg1,dtype=np.float64)
        jv_source = jt_source + 0.5

# Generate a set of target points in *source* j-space based on the resolution ratio 
# calculated above, then work out interpolation coefficients for the grid 
# parameters. Note that zeroth T-U line on the target grid segment is the join-line
# between the two grids. 

        jt_target = n_halo + np.arange(nj_seg[iceshelf])/resolution_ratio
        jv_target = n_halo + ( np.arange(nj_seg[iceshelf]) + 0.5 )/resolution_ratio

        dj_target_fine = 1.0/(resolution_ratio*jpfac)
        print 'dj_target_fine :',dj_target_fine
        
# For each target grid point (T/U or V/F), 4th order interpolation 
# in the j-direction will always use 2 source T/U points and 
# 2 source V/F points. 

#----------------------------------------------
# indices of source points (counting from zero)
#----------------------------------------------
        jtu1 = np.zeros([nj_seg[iceshelf],npoints,2],dtype=np.uint16)
        jtu2 = np.zeros([nj_seg[iceshelf],npoints,2],dtype=np.uint16)
        jvf1 = np.zeros([nj_seg[iceshelf],npoints,2],dtype=np.uint16)
        jvf2 = np.zeros([nj_seg[iceshelf],npoints,2],dtype=np.uint16)

        for ij in range(nj_seg[iceshelf]):
            for ip in np.arange(npoints)-2:
                jt_tg = jt_target[ij] + ip * dj_target_fine
                jv_tg = jv_target[ij] + ip * dj_target_fine

# T/U target points:
#------------------
# nearest T/U source point:
                jtu1[ij,ip,0] = np.argmin(abs(jt_tg - jt_source))    
# next-nearest T/U source point:
                jtu2[ij,ip,0] = jtu1[ij,ip,0] + math.copysign(1, jt_tg-jt_source[jtu1[ij,ip,0]])             
# nearest V/F source point:
                jvf1[ij,ip,0] = np.argmin(abs(jt_tg - jv_source))    
# next-nearest V/F source point:
                jvf2[ij,ip,0] = jvf1[ij,ip,0] + math.copysign(1, jt_tg-jv_source[jvf1[ij,ip,0]])             

# V/F target points:
#------------------
# nearest T/U source point:
                jtu1[ij,ip,1] = np.argmin(abs(jv_tg - jt_source))    
# next-nearest T/U source point:
                jtu2[ij,ip,1] = jtu1[ij,ip,1] + math.copysign(1, jv_tg-jt_source[jtu1[ij,ip,1]])             
# nearest V/F source point:
                jvf1[ij,ip,1] = np.argmin(abs(jv_tg - jv_source))    
# next-nearest V/F source point:
                jvf2[ij,ip,1] = jvf1[ij,ip,1] + math.copysign(1, jv_tg-jv_source[jvf1[ij,ip,1]])             

#            print 'ij, jvf1[ij,1],  jvf2[ij,1] :',ij, jvf1[ij,1],  jvf2[ij,1]

        print 'min/max jtu1 : ',jtu1.min(),jtu1.max()
        print 'min/max jtu2 : ',jtu2.min(),jtu2.max()
        print 'min/max jvf1 : ',jvf1.min(),jvf1.max()
        print 'min/max jvf2 : ',jvf2.min(),jvf2.max()

#-----------------------------
# interpolation coefficients
#-----------------------------
        ctu1 = np.zeros([nj_seg[iceshelf],npoints,2])
        ctu2 = np.zeros([nj_seg[iceshelf],npoints,2])
        cvf1 = np.zeros([nj_seg[iceshelf],npoints,2])
        cvf2 = np.zeros([nj_seg[iceshelf],npoints,2])

        for ij in range(nj_seg[iceshelf]):
            for ip in np.arange(npoints)-2:

# T/U target points:
#------------------
                jt_tgt = jt_target[ij] + ip * dj_target_fine
                jt_src1 = jt_source[jtu1[ij,ip,0]]
                jt_src2 = jt_source[jtu2[ij,ip,0]]
                jv_src1 = jv_source[jvf1[ij,ip,0]]
                jv_src2 = jv_source[jvf2[ij,ip,0]]

                ctu1[ij,ip,0] = ( (jt_tgt -jt_src2)*(jt_tgt -jv_src1)*(jt_tgt -jv_src2) ) \
                              / ( (jt_src1-jt_src2)*(jt_src1-jv_src1)*(jt_src1-jv_src2) )
                ctu2[ij,ip,0] = ( (jt_tgt -jt_src1)*(jt_tgt -jv_src1)*(jt_tgt -jv_src2) ) \
                              / ( (jt_src2-jt_src1)*(jt_src2-jv_src1)*(jt_src2-jv_src2) )
                cvf1[ij,ip,0] = ( (jt_tgt -jt_src1)*(jt_tgt -jt_src2)*(jt_tgt -jv_src2) ) \
                              / ( (jv_src1-jt_src1)*(jv_src1-jt_src2)*(jv_src1-jv_src2) )
                cvf2[ij,ip,0] = ( (jt_tgt -jt_src1)*(jt_tgt -jt_src2)*(jt_tgt -jv_src1) ) \
                              / ( (jv_src2-jt_src1)*(jv_src2-jt_src2)*(jv_src2-jv_src1) )

# V/F target points:
#------------------
                jv_tgt = jv_target[ij] + ip * dj_target_fine
                jt_src1 = jt_source[jtu1[ij,ip,1]]
                jt_src2 = jt_source[jtu2[ij,ip,1]]
                jv_src1 = jv_source[jvf1[ij,ip,1]]
                jv_src2 = jv_source[jvf2[ij,ip,1]]

                ctu1[ij,ip,1] = ( (jv_tgt -jt_src2)*(jv_tgt -jv_src1)*(jv_tgt -jv_src2) ) \
                              / ( (jt_src1-jt_src2)*(jt_src1-jv_src1)*(jt_src1-jv_src2) )
                ctu2[ij,ip,1] = ( (jv_tgt -jt_src1)*(jv_tgt -jv_src1)*(jv_tgt -jv_src2) ) \
                              / ( (jt_src2-jt_src1)*(jt_src2-jv_src1)*(jt_src2-jv_src2) )
                cvf1[ij,ip,1] = ( (jv_tgt -jt_src1)*(jv_tgt -jt_src2)*(jv_tgt -jv_src2) ) \
                              / ( (jv_src1-jt_src1)*(jv_src1-jt_src2)*(jv_src1-jv_src2) )
                cvf2[ij,ip,1] = ( (jv_tgt -jt_src1)*(jv_tgt -jt_src2)*(jv_tgt -jv_src1) ) \
                              / ( (jv_src2-jt_src1)*(jv_src2-jt_src2)*(jv_src2-jv_src1) )

###################################################################
# 10. Use the coefficients calculated in section 9 to interpolate
#     in the j-direction and create the final grid segment for
#     this ice shelf.  
###################################################################

        for ii in range(ni_seg[iceshelf]):
            for ip in np.arange(npoints)-2:
 
                lamt_seg[iceshelf][:,:,ii,ip] = ctu1[:,:,0] * lamt_seg1[jtu1[:,:,0],ii,ip] \
                                              + ctu2[:,:,0] * lamt_seg1[jtu2[:,:,0],ii,ip] \
                                              + cvf1[:,:,0] * lamv_seg1[jvf1[:,:,0],ii,ip] \
                                              + cvf2[:,:,0] * lamv_seg1[jvf2[:,:,0],ii,ip] 

                lamu_seg[iceshelf][:,:,ii,ip] = ctu1[:,:,1] * lamt_seg1[jtu1[:,:,1],ii,ip] \
                                              + ctu2[:,:,1] * lamt_seg1[jtu2[:,:,1],ii,ip] \
                                              + cvf1[:,:,1] * lamv_seg1[jvf1[:,:,1],ii,ip] \
                                              + cvf2[:,:,1] * lamv_seg1[jvf2[:,:,1],ii,ip] 

                lamv_seg[iceshelf][:,:,ii,ip] = ctu1[:,:,0] * lamu_seg1[jtu1[:,:,0],ii,ip] \
                                              + ctu2[:,:,0] * lamu_seg1[jtu2[:,:,0],ii,ip] \
                                              + cvf1[:,:,0] * lamf_seg1[jvf1[:,:,0],ii,ip] \
                                              + cvf2[:,:,0] * lamf_seg1[jvf2[:,:,0],ii,ip] 

                lamf_seg[iceshelf][:,:,ii,ip] = ctu1[:,:,1] * lamu_seg1[jtu1[:,:,1],ii,ip] \
                                              + ctu2[:,:,1] * lamu_seg1[jtu2[:,:,1],ii,ip] \
                                              + cvf1[:,:,1] * lamf_seg1[jvf1[:,:,1],ii,ip] \
                                              + cvf2[:,:,1] * lamf_seg1[jvf2[:,:,1],ii,ip] 

                phit_seg[iceshelf][:,:,ii,ip] = ctu1[:,:,0] * phit_seg1[jtu1[:,:,0],ii,ip] \
                                              + ctu2[:,:,0] * phit_seg1[jtu2[:,:,0],ii,ip] \
                                              + cvf1[:,:,0] * phiv_seg1[jvf1[:,:,0],ii,ip] \
                                              + cvf2[:,:,0] * phiv_seg1[jvf2[:,:,0],ii,ip] 

                phiu_seg[iceshelf][:,:,ii,ip] = ctu1[:,:,1] * phit_seg1[jtu1[:,:,1],ii,ip] \
                                              + ctu2[:,:,1] * phit_seg1[jtu2[:,:,1],ii,ip] \
                                              + cvf1[:,:,1] * phiv_seg1[jvf1[:,:,1],ii,ip] \
                                              + cvf2[:,:,1] * phiv_seg1[jvf2[:,:,1],ii,ip] 

                phiv_seg[iceshelf][:,:,ii,ip] = ctu1[:,:,0] * phiu_seg1[jtu1[:,:,0],ii,ip] \
                                              + ctu2[:,:,0] * phiu_seg1[jtu2[:,:,0],ii,ip] \
                                              + cvf1[:,:,0] * phif_seg1[jvf1[:,:,0],ii,ip] \
                                              + cvf2[:,:,0] * phif_seg1[jvf2[:,:,0],ii,ip] 

                phif_seg[iceshelf][:,:,ii,ip] = ctu1[:,:,1] * phiu_seg1[jtu1[:,:,1],ii,ip] \
                                              + ctu2[:,:,1] * phiu_seg1[jtu2[:,:,1],ii,ip] \
                                              + cvf1[:,:,1] * phif_seg1[jvf1[:,:,1],ii,ip] \
                                              + cvf2[:,:,1] * phif_seg1[jvf2[:,:,1],ii,ip] 

######################################################################
# 11. Calculate the e1/e2 scale factors using 4th-order differencing
#     on the fine grid. 
######################################################################

        print 'phit_seg[iceshelf][0,0,0,0],phit_seg[iceshelf][1,0,0,0],phit_seg[iceshelf][2,0,0,0], :',\
               phit_seg[iceshelf][0,0,0,0],phit_seg[iceshelf][1,0,0,0],phit_seg[iceshelf][2,0,0,0]

        e1t_seg[iceshelf] = e1_calc(lamt_seg[iceshelf],phit_seg[iceshelf],jpfac)
        e1u_seg[iceshelf] = e1_calc(lamu_seg[iceshelf],phiu_seg[iceshelf],jpfac)
        e1v_seg[iceshelf] = e1_calc(lamv_seg[iceshelf],phiv_seg[iceshelf],jpfac)
        e1f_seg[iceshelf] = e1_calc(lamf_seg[iceshelf],phif_seg[iceshelf],jpfac)

        e2t_seg[iceshelf] = e2_calc(lamt_seg[iceshelf],phit_seg[iceshelf],jpfac)
        e2u_seg[iceshelf] = e2_calc(lamu_seg[iceshelf],phiu_seg[iceshelf],jpfac)
        e2v_seg[iceshelf] = e2_calc(lamv_seg[iceshelf],phiv_seg[iceshelf],jpfac)
        e2f_seg[iceshelf] = e2_calc(lamf_seg[iceshelf],phif_seg[iceshelf],jpfac)

        ortho = orthogonality_check(lamt_seg[iceshelf],phit_seg[iceshelf],e1t_seg[iceshelf],e2t_seg[iceshelf],jpfac)

        print 'ortho.shape : ',ortho.shape

        write_coords('coordsT_'+iceshelf+'.nc',glamt_base.dtype,glamt=lamt_seg[iceshelf][:,0,:,0],gphit=phit_seg[iceshelf][:,0,:,0],
                                                                e1t=e1t_seg[iceshelf], e2t=e2t_seg[iceshelf], ortho=ortho)

###################################################################
# 12. Create the extended grid and slot in the ice-shelf segments. 
###################################################################

# Dimensions of extended grid:
    nj_south = max(nj_seg.values())
    ni_ext = lamt_base.shape[1]
    nj_ext = lamt_base.shape[0] + nj_south

    print '=================='
    print '>>> Final grid <<<'
    print '=================='
    print 'ni_ext :',ni_ext
    print 'nj_ext :',nj_ext

    lamt_ext = np.zeros([nj_ext,ni_ext])
    lamu_ext = np.zeros([nj_ext,ni_ext])
    lamv_ext = np.zeros([nj_ext,ni_ext])
    lamf_ext = np.zeros([nj_ext,ni_ext])

    phit_ext = np.zeros([nj_ext,ni_ext])
    phiu_ext = np.zeros([nj_ext,ni_ext])
    phiv_ext = np.zeros([nj_ext,ni_ext])
    phif_ext = np.zeros([nj_ext,ni_ext])

    e1t_ext = np.zeros([nj_ext,ni_ext])
    e1u_ext = np.zeros([nj_ext,ni_ext])
    e1v_ext = np.zeros([nj_ext,ni_ext])
    e1f_ext = np.zeros([nj_ext,ni_ext])

    e2t_ext = np.zeros([nj_ext,ni_ext])
    e2u_ext = np.zeros([nj_ext,ni_ext])
    e2v_ext = np.zeros([nj_ext,ni_ext])
    e2f_ext = np.zeros([nj_ext,ni_ext])

# Set equal to the base coordinates over most of the domain:
# Switch longitude range back to [-180,+180]

    lamt_ext[nj_south:,:] = fix_range(lamt_base,[-180.0,+180.0])
    lamu_ext[nj_south:,:] = fix_range(lamu_base,[-180.0,+180.0])
    lamv_ext[nj_south:,:] = fix_range(lamv_base,[-180.0,+180.0])
    lamf_ext[nj_south:,:] = fix_range(lamf_base,[-180.0,+180.0])

    phit_ext[nj_south:,:] = phit_base
    phiu_ext[nj_south:,:] = phiu_base
    phiv_ext[nj_south:,:] = phiv_base
    phif_ext[nj_south:,:] = phif_base

    e1t_ext[nj_south:,:] = e1t_base
    e1u_ext[nj_south:,:] = e1u_base
    e1v_ext[nj_south:,:] = e1v_base
    e1f_ext[nj_south:,:] = e1f_base

    e2t_ext[nj_south:,:] = e2t_base
    e2u_ext[nj_south:,:] = e2u_base
    e2v_ext[nj_south:,:] = e2v_base
    e2f_ext[nj_south:,:] = e2f_base

# Slot in the ice shelf segments:

    for ii_shelf in range(len(lam_west)):

        iceshelf = lam_west.keys()[ii_shelf]

        for ij in range(nj_seg[iceshelf]):
            
            lamt_ext[nj_south-ij,ii_seg[iceshelf]] = fix_range( lamt_seg[iceshelf][ij,0,:,0],[-180.0,+180.0] )
            lamu_ext[nj_south-ij,ii_seg[iceshelf]] = fix_range( lamu_seg[iceshelf][ij,0,:,0],[-180.0,+180.0] )
            lamv_ext[nj_south-ij,ii_seg[iceshelf]] = fix_range( lamv_seg[iceshelf][ij,0,:,0],[-180.0,+180.0] )
            lamf_ext[nj_south-ij,ii_seg[iceshelf]] = fix_range( lamf_seg[iceshelf][ij,0,:,0],[-180.0,+180.0] )

            phit_ext[nj_south-ij,ii_seg[iceshelf]] = phit_seg[iceshelf][ij,0,:,0]
            phiu_ext[nj_south-ij,ii_seg[iceshelf]] = phiu_seg[iceshelf][ij,0,:,0]
            phiv_ext[nj_south-ij,ii_seg[iceshelf]] = phiv_seg[iceshelf][ij,0,:,0]
            phif_ext[nj_south-ij,ii_seg[iceshelf]] = phif_seg[iceshelf][ij,0,:,0]

            e1t_ext[nj_south-ij,ii_seg[iceshelf]] = e1t_seg[iceshelf][ij,:]
            e1u_ext[nj_south-ij,ii_seg[iceshelf]] = e1u_seg[iceshelf][ij,:]
            e1v_ext[nj_south-ij,ii_seg[iceshelf]] = e1v_seg[iceshelf][ij,:]
            e1f_ext[nj_south-ij,ii_seg[iceshelf]] = e1f_seg[iceshelf][ij,:]

            e2t_ext[nj_south-ij,ii_seg[iceshelf]] = e2t_seg[iceshelf][ij,:]
            e2u_ext[nj_south-ij,ii_seg[iceshelf]] = e2u_seg[iceshelf][ij,:]
            e2v_ext[nj_south-ij,ii_seg[iceshelf]] = e2v_seg[iceshelf][ij,:]
            e2f_ext[nj_south-ij,ii_seg[iceshelf]] = e2f_seg[iceshelf][ij,:]


    write_coords('coords_extended.nc',glamt_base.dtype,glamt=lamt_ext,glamu=lamu_ext,glamv=lamv_ext,glamf=lamf_ext,
                                                       gphit=phit_ext,gphiu=phiu_ext,gphiv=phiv_ext,gphif=phif_ext,
                                                       e1t=e1t_ext,e1u=e1u_ext,e1v=e1v_ext,e1f=e1f_ext,
                                                       e2t=e2t_ext,e2u=e2u_ext,e2v=e2v_ext,e2f=e2f_ext,
                                                       aniso=True)
    


    base_dataset.close()
    north_dataset.close()


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("base_coords", help="name of base coordinates file")
    parser.add_argument("north_coords", help="name of coordinates file from which to extract north section")
    parser.add_argument("-c", "--lat_cut", action="store",type=float,dest="lat_cut",
                    default=None,help="set cut off latitude")
    parser.add_argument("-f", "--jpfac", action="store",type=float,dest="jpfac",
                    default=10.0,help="resolution factor between target grid and grid for calculation of scale factors")

    args = parser.parse_args()

    extend_orca(args.base_coords, args.north_coords, lat_cut=args.lat_cut, jpfac=args.jpfac)
                        

