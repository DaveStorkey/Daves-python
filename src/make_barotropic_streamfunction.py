#! /usr/bin/env python

'''
Routine to calculate the quasi barotropic streamfunction from the vertically
integrated zonal mass transport output from model.

Or alternatively from the velocity, thickness and e2u fields.  

@author: Dave Storkey and Tim Graham.
'''

import iris
import iris.analysis

def make_barotropic_streamfunction(infilename,coordsfilename=None,outfilename=None,use_vel=False):

    if "grid_U" in infilename or "grid-U" in infilename:
        print "Using grid U variables"
        dirn = "x"
        vely  = "u"
        nemo_vel = "vozocrtx"
    elif "grid_V" in infilename or "grid-V" in infilename:
        print "Using grid V variables"
        dirn = "y"
        vely  = "v"
        nemo_vel = "vomecrty"
    else:
        raise Exception("Can only handle files called *grid_U.nc, *grid-U.nc, *grid_V.nc or *grid-V.nc or I get confused...")

    if outfilename is None:
        outfilename = infilename.replace(".nc","_barostrmfn.nc")

    if not use_vel:
        try:
            masstr_vint = iris.load_cube(infilename,"vertical_integral_of_ocean_mass_"+dirn+"_transport")    
        except iris.exceptions.ConstraintMismatchError:
            # if we can't match the standard name try matching the internal variable name     
           try: 
                masstr_vint = iris.load_cube(infilename,iris.Constraint(cube_func=lambda nemoname: nemoname.var_name == vely+"_masstr_vint") )
           except iris.exceptions.ConstraintMismatchError:
                print "Can't find vertically integrated mass transport in file."
                print "Attempting to calculate from velocity and thickness fields."
                use_vel = True        

    if use_vel:
        # Try to read in e2u from the coordinates file.
        if coordsfilename is None:
            raise Exception("Coordinates file required to calculate streamfunction from velocity field.")

        e2 = iris.load_cube(coordsfilename,iris.Constraint(cube_func=lambda nemoname: nemoname.var_name == "e2"+vely) )

        # Read in velocities
        try:
            vel = iris.load_cube(infilename,"sea_water_"+dirn+"_velocity")    
        except iris.exceptions.ConstraintMismatchError:
            # if we can't match the standard name try matching the internal variable name     
            for varname in [nemo_vel,vely+"o"]:
                try:
                    vel = iris.load_cube(infilename,iris.Constraint(cube_func=lambda nemoname: nemoname.var_name == varname) )
                except iris.exceptions.ConstraintMismatchError:
                    pass
                else:
                    break
            else:
                raise Exception("Can't find velocity field in file.")

        try:
            thick = iris.load_cube(infilename,"cell_thickness")    
        except iris.exceptions.ConstraintMismatchError:
            # if we can't match the standard name try matching the internal variable name     
            for varname in ["e3"+vely,"thkcello"]:
                try:
                    thick = iris.load_cube(infilename,iris.Constraint(cube_func=lambda nemoname: nemoname.var_name == varname) )
                except iris.exceptions.ConstraintMismatchError:
                    pass
                else:
                    break
            else:
                raise Exception("Can't find cell thickness field in file.")

        # 1026. kg/m3 is the reference density of seawater as used in NEMO. 
        masstr_vint = ( vel * thick * e2 * 1026. ).collapsed('depth',iris.analysis.SUM)           
        masstr_vint.var_name = vely+"masstr_vint"
        #masstr_vint.standard_name = "vertical_integral_of_ocean_mass_"+dirn+"_transport"

    barostrmfn = ocean_quasi_barotropic_streamfunc(masstr_vint,velocity=vely)

    iris.save([masstr_vint,barostrmfn], outfilename)

def ocean_quasi_barotropic_streamfunc(incube,velocity=None):
    """
    Computes the quasi-barotropic streamfunction from vertical integral
    of zonal mass transport by doing a cumulative sum over the y axis and
    multiplying by -1.

    Parameters
    ----------
    cube: :class:`iris.cube.Cube`
        A cube containing vertically integrated zonal mass transport

    Returns
    -------
    : :class:`iris.cube.Cube`
       A cube of quasibarotropic streamfunction
    """

    if velocity == "u":
        dim = incube.coord_dims(incube.coord(axis='Y'))[0]
    elif velocity == "v":
        dim = incube.coord_dims(incube.coord(axis='X'))[0]
    data_array = incube.data.cumsum(axis=dim)
    # 1026 is the density - converting from mass transport to volume transport.
    # convert to Sverdrups while we're at it
    outcube = incube.copy(data=data_array * -1.e-06/1026.)
    outcube.var_name="barostrmfn"
    outcube.standard_name="ocean_barotropic_streamfunction"
    outcube.units="Sv"

    return outcube

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("infilename", help="name of input file")
    parser.add_argument("-c", "--coordsfile", action="store",dest="coordsfilename",default=None,
                    help="coordinates file. Required if using velocity field to calculate streamfunction.")
    parser.add_argument("-o", "--outfile", action="store",dest="outfilename",default=None,
                    help="output file. If unspecified write to original file.")
    parser.add_argument("-v", "--use_vel", action="store_true",dest="use_vel",default=None,
                    help="use velocities and thicknesses even if u/v_masstr_vint is in file.")

    args = parser.parse_args()

    make_barotropic_streamfunction(args.infilename,coordsfilename=args.coordsfilename,outfilename=args.outfilename,use_vel=args.use_vel)
