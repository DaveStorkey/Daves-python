'''
Created on Feb 12, 2013

This version August 2017

@author: Tim Graham and Dave Storkey
'''

import numpy as np
import netCDF4
import iris

def read_cube(filename,fieldname):
    '''
    Read a variable from a netcdf file as an Iris cube.
    Try matching varname to the standard_name first and if that fails
    try matching to the variable name in the file. 
    '''

    try:
        cube = iris.load_cube(filename, fieldname)
    except iris.exceptions.ConstraintMismatchError:
        try: 
            cube = iris.load_cube(filename, iris.Constraint(cube_func=lambda nemoname: nemoname.var_name == fieldname))
        except iris.exceptions.ConstraintMismatchError:
            raise Exception('Could not find '+fieldname+' in file '+filename)

    return cube

def read_nc_var(fname,var):
    '''
    Read a variable from a netcdf file.
    Inputs are filename and variable name
    '''
    fileid=netCDF4.Dataset(fname)
    data=fileid.variables[var][:].copy()
    fileid.close()
    return data

def copy_attributes(var_in,var_out):
    '''
    Copy attributes from one netCDF4 variable class instance
    to another. Exclude _FillValue attribute as this is created
    in the createVariable statement
    '''
    attrdict={}
    for attr in var_in.ncattrs():
        attrdict[attr.encode('ascii')] = var_in.getncattr(attr)
    del attrdict['_FillValue']
    var_out.setncatts(attrdict)
    

def wrap_longitude(longitudes,base=0.0):
    '''
    Convert points in longitude to range from base (default zero) to base+360 
    ''' 
    if isinstance(longitudes,list): 
        longitudes=np.array(longitudes)
        intype='list'
    elif isinstance(longitudes,tuple): 
        longitudes=np.array(longitudes)
        intype='tuple'
    else: intype=None
    wrapped_longs=((longitudes-base+720.) % 360.)+base
    
    if intype=='list': wrapped_longs=wrapped_longs.tolist()
    if intype=='tuple': wrapped_longs=tuple(wrapped_longs.tolist())
    
    return wrapped_longs

def mask_gen(lon,lat,region):
    '''
    Generate a numpy boolean array consisting of False within region and True outside of the region
    Inputs arguments:
    lon = 1D or 2D numpy array of longitude
    lat = 1D or 2D numpy array of latitude
    region = tuple or list of format [North, East, South, West] specifying the limits of the region
    Aug 2017 : allow None to be specified for some region elements
    '''
    
    #Check input types
    if not isinstance(lat,np.ndarray) or (lat.ndim not in [1,2]):
        raise RuntimeError("lat must be a 1 dimensional or 2dimensional numpy array")
    if not isinstance(lon,np.ndarray) or (lon.ndim not in [1,2]):
        raise RuntimeError("lon must be a 1 dimensional or 2 dimensional numpy array")
    if not isinstance(region,(tuple,list)) or (len(region) != 4):
        raise RuntimeError("Region must be 4 element tuple or list")
    if lon.shape != lat.shape:
        raise RuntimeError("lat and lon must have same shape")
    
    #Wrap longitude to be in range region[3]:region[3]+360
    lon=wrap_longitude(lon,base=region[3])
    
    # Allow None to be specified for regional elements
    if region[0] is None:
        region[0] = np.amax(lat)
    if region[1] is None:
        region[1] = np.amax(lon)
    if region[2] is None:
        region[2] = np.amin(lat)
    if region[3] is None:
        region[3] = np.amin(lon)

    #Search for lat and lon points within range and set to 1
    latmask=np.where(lat < region[0],1,0)*np.where(lat > region[2],1,0)
    lonmask=np.where(lon < region[1],1,0)*np.where(lon > region[3],1,0)
    
    #Combine latmask and lonmask
    if lon.ndim == 1:
        pass
    else:
        mask=np.where(latmask*lonmask == 0,True,False)

    return mask

def combine_masks(mask1,mask2):
    '''
    Combine two masks together.
    The output will be true whereever one of the masks is true
    Masks must be commutative but don't have to have same number of dimensions
    e.g. a 2D basin mask can be combined with a 3D land-sea mask
    '''
    #mask1=np.where(mask1,0,1)
    #mask2=np.where(mask2,0,1)
    #mask1=mask1*-1+1
    #mask2=mask2*-1+1
    #newmask=np.where(mask1*mask2 == 0,True,False)
    return np.array((mask1*-1+1)*(mask2*-1+1)*-1+1,dtype=bool)

def trim_cube(cube,dims=None):
    '''
    Trim an iris cube along all dimensions to the minimum and maximum indexes where valid data occurs
    This is useful for grids such as NEMO where we might want to mask based on lat/lon in northern hemisphere
    
    Returns a tuple of the modified cube and the slice object that can be used to slice further identical cubes
    '''
    cube_index = [slice(None)]*cube.ndim
    if not np.ma.isMA(cube.data):
        return cube,cube_index
    
    mask_inv=np.where(cube.data.mask,0,1) #It's easier to work with inverse of mask
    #Find the shape and number of dimensions
    shape=cube.shape
    ndims=cube.ndim
    if not dims: dims=range(ndims)
    
    for dim in dims:
        if shape[dim] == 1:
            #Nothing to do for a unitary dimension
            pass
        else:
            #Sum over remaining dimensions
            axes=range(ndims)#.remove(dim)
            if dim < 0:
                dim=ndims+dim
            axes.remove(dim)
            newsum=np.apply_over_axes(np.sum, mask_inv, axes).squeeze()
            nodata=np.where(newsum == 0,0,1)
            min_point=nodata.argmax() #The first point where there is data
            max_point=shape[dim]-nodata[::-1].argmax()
            if min_point != 0 or max_point != shape[dim]: #See if we need to extract data for this dimension
                cube_index[dim] = np.arange(min_point,max_point)
                
    return (cube[tuple(cube_index)], cube_index)

def rms_diff(cube1,cube2,coords,weights=None,mask=None):
    '''
    Purpose: Calculate the spatial RMS difference between two cubes
    Inputs: cube1 - IRIS cube 
            cube2 - IRIS cube
            coords - a list of coordinate names over which to calculate RMS to pass to iris.cube.collapsed
    Optional inputs: weights - numpy array of weights (e.g. for area weighted RMS)
                     mask - numpy array of ones (valid points) and zeros (masked points)
    Output: (Weighted) Cube of RMS difference between the cubes 
    '''
    assert cube1.shape == cube2.shape, "Both cubes must have same shape"
    
    if weights is not None:
        assert  cube1.shape == weights.shape, 'Weights must have same shape as cubes'
        
        #Mask the weights using mask from cube1
        miss=np.where(cube1.data.mask,0,1)
        weights*=miss
        
        if mask is not None:
            weights*=mask
        
        #Find out the dimensions of coords so we can collapse weights in same way
        rms=(((cube1-cube2)**2)*weights).collapsed(coords,iris.analysis.SUM)
        
        dims=[]
        for coordname in coords:
            coord=cube1.coord(coordname)
            d=cube1.coord_dims(coord)
            for dim in d: dims.append(dim)
        dims=list(set(dims)) #Make into a list of unique dims
        dims.sort()
        dims.reverse()
        ndims=cube1.ndim
        for dim in dims:
            #Work from highest numbered dimension first and back through array using negative dimension indices
            weights=weights.sum(axis=ndims-dim-2)
            ndims=ndims-1
        
        rms=(rms/weights)**0.5
    else:
        if mask is not None:
            cube1*=mask
            cube2*=mask
        rms=(cube1-cube2).collapsed(coords,iris.analysis.RMS)
        
    return rms

def grid_zonm(cube, avg_lat=True, mask=None, dim=-1, weights=None):
    '''
    Purpose: Calculate the mean along grid lines of a cube on an irregular grid
              rather than a true zonal mean
    Inputs: cube - an iris cube
    Optional Inputs: avg_lat - Make the X coordinate of the cube the average latitude of the full field
                     mask - A mask of zeros (invalid points) and ones (valid points)
                            Must be same shape as cube
                     dim - Dimension over which to average (default is last)
                     weights - An array of weights to be passed in to the averaging (e.g the width of the grid cells)
                               Passed directly to iris.cube.collapsed( )
    Output: A new IRIS cube
    
    Method: Add a new dim_coord to the input cube and collapse cube along this dimension
    '''
    
    #Make dim positive
    if dim < 0: dim = cube.ndim + dim
    dim_len=cube.shape[dim]
    
    xcoord=iris.iris.coords.DimCoord(np.arange(dim_len), long_name='X')
    cube=cube.copy()
    #Avoid modifying the cube in place
    cube.add_dim_coord(xcoord,dim)
    
    if avg_lat:
        lat=cube.coord('latitude').copy()
        lat_dims=cube.coord_dims(lat)
        lat_dims=list(lat_dims)
        lat_dims.remove(dim) #lat_dims should now be the dimension to add the 1D lat coord to (unless dim is < lat_dims)
        lat_dims=lat_dims[0]
        if lat_dims > dim: lat_dims -= 1
        latshape=list(lat.points.shape)
        # Could add bounds here?
        lat=lat.points.mean(axis=latshape.index(cube.shape[dim]))
        lat_coord=iris.coords.AuxCoord(lat,long_name='Latitude',standard_name='latitude',units='degrees')
        
    #Remove latitude and longitude coordinates
    cube.remove_coord(cube.coord('latitude'))
    cube.remove_coord(cube.coord('longitude'))
    
    if mask is not None:
        assert mask.shape==cube.shape, "Mask must be same shape as cube. Hint: Use iris.util.broadcast_weights"
        ind=np.where(mask==0)
        cube.data.mask[ind]=True
        
    cube=cube.collapsed('X',iris.analysis.MEAN,weights=weights)
    
    if avg_lat:
        cube.add_aux_coord(lat_coord,data_dims=lat_dims)
        
    return cube

#------------------------------------------------------------------------------------------------
def set_ORCA_halos(field, tuvf, isign, pole_type ): 

# @author : Mike Bell
#
# inputs
# ------
# field                 2D numpy array of floats on an ORCA grid including halos (halo values might be incorrect on input)  
# tuvf                  string: 't', 'u', 'v' or 'f' depending on the grid on which field is valid 
# isign                 +1. for true scalars (e.g. tracers, divergence); -1 for true vectors (velocities) or pseudo scalars (vorticity) 
#                       to set halos for a tracer field that has been averaged onto the U grid, use tuvf = 'u' and isign = +1.                     
# pole_type             'f' or 't' depending on the model configuration: ORCA2, ORCA025 and ORCA12 use a 't' point pivot; ORCA1 an 'f' type

# returns 
# -------
# field                 2D numpy array of floats on an ORCA grid with the halos set correctly   

# The derivation of this code from the Fortran is described in ../@Notes/lbc_nfd_generic_convert_to_python.txt 


   import sys
   import numpy as np

# check that the inputs for pole_type and tuvf are recognised
   if pole_type != 't' and pole_type != 'f' : 
      print ( ' set_ORCA_halos: pole_type is set incorrectly: it should be \'t\' or \'f\'; it is ', pole_type )   
      sys.exit()

   if tuvf != 't' and tuvf != 'u' and tuvf != 'v' and tuvf != 'f' and tuvf != 'w' : 
      print ( ' set_ORCA_halos: tuvf is set incorrectly: it should be \'t\', \'u\', \'v\' or \'f\'; it is ', tuvf )   
      sys.exit()

# set the cyclic boundary conditions 
   field[:, 0] = field[:,-2]
   field[:,-1] = field[:, 0]

   jpj, jpi = np.shape(field)     # jpi is required below 
   
   if pole_type == 't' :       

      if tuvf == 't' or tuvf == 'w' : 
         field[-1, 1:] = isign * field [ -3, -1:0:-1 ]   
         field[-1, 0]  = isign * field [ -3, 2 ]           
         field[-2, int(jpi/2):] = field[-2, int((jpi+1)/2): 0: -1]   
	 
      elif tuvf == 'u' :      # isign = +1. for vectors like velocity
         field[-1, :-1] = isign * field[-3, :0:-1]  	 
         field[-1, 0 ]  = isign * field[-3, 1]
         field[-1, -1]  = isign * field[-3, -2]
         field[-2, int(jpi/2)-1:-1] = - field[-2, int((jpi+1)/2):0:-1]   
	 
      elif tuvf == 'v' : 
         field[-2, 1:] = isign * field [ -3, -1:0:-1 ]     # rhs corrected from , -1:1:-1] to -1:0:-1] on 31/03/21 
         field[-1, 1:] = isign * field [ -4, -1:0:-1 ]     # same correction as above
         field[-1, 0 ] = isign * field [ -4, 2 ] 

      elif tuvf == 'f' : 
         field[-2, :-1] = isign * field[-3, -1:0:-1]
         field[-1, :-1] = isign * field[-4, -1:0:-1]
         field[-1,0]  = isign * field[-4,1]
         field[-1,-1] = isign * field[-4,-2]

   if pole_type == 'f' :  

      if tuvf == 't' or tuvf == 'w' : 
         field[-1, : ] = isign * field[ -2, ::-1 ]

      elif tuvf == 'u' : 
         field[-1, :-1] = isign * field[-2, -2::-1]  	 
         field[-1,-1]   = isign *  field[-2,-3]	 

      elif tuvf == 'v' : 
         field[-1, :] = isign * field[-3, ::-1] 
         field[-2, int(jpi/2):] = isign * field[-2, int((jpi-1)/2)::-1] 

      elif tuvf == 'f' : 
         field[-1, :-1] = isign * field[-3, -2::-1]
         field[-1,-1] = isign * field[-3, -3]
         field[-2, int(jpi/2):-1] = isign * field[-2, int((jpi-1)/2)-1::-1] 

   return field 

    
