#! /usr/bin/env python
'''
Wrapper routine for plot_nemo_section.py to plot
trends in the zonal velocity (accelerations) along
sections in the equatorial Pacific (and Atlantic) 
overlaid on the zonal velocity. If no trend specified
just plot zonal velocity as colours and lines.

Dave Storkey
July 2019
'''

from plot_nemo_section import plot_nemo_section

def plot_utrd_section(filestem=None,filenames=None,utrd=None,section=None,cmap=None,mnfld=None,mxfld=None,
                      factor=None,toldeg=None):

    if filestem is not None:
        if utrd is None:
            filenames=filestem+"_grid-U.nc"
        elif 'utrd_adv' in utrd:
            filenames=[filestem+"_utrd_adv.nc",filestem+"_grid-U.nc"]
        elif 'utrd' in utrd:
            filenames=filestem+"_grid-U.nc"
        else:
            filenames=[filestem+"_mmtm_adv_accel.nc",filestem+"_grid-U.nc"]
    elif filenames is None:
        raise Exception("Error: must specify filestem or filenames.")
    else:
         filestem = filenames[-1].replace("_grid-U.nc","")

    if utrd is None:
        var_names=["uo","uo"]
        scientific=False
        if mnfld is None:
            mnfld=[-1.5,-1.5]
        if mxfld is None:
            mxfld=[1.5,1.5]
    elif 'over' in utrd:
        var_names=[utrd,"uo"]
        scientific=True
        if mnfld is None:
            mnfld=[-1.0e-02,-1.5]
        if mxfld is None:
            mxfld=[1.0e-02,1.5]
    else:
        var_names=[utrd,"uo"]
        scientific=True
        if mnfld is None:
            mnfld=[-1.0e-06,-1.5]
        if mxfld is None:
            mxfld=[1.0e-06,1.5]

    # plot convergence (=acceleration) rather than divergence.
    if 'div' in utrd:
        factor=[-1.0,1.0]

    if toldeg is None:
        toldeg = 0.2

    if section is None:
        # default to section along the equator in the Pacific
        section="eqpac"

    if section == "eqpac":
        lat=0.0
        lon=None
        xmin=140.0
        xmax=270
        sectitle=" along the equator"
    elif section == "135W":
        lat=None
        lon=-135.0
        xmin=-10.0
        xmax=10.0
        sectitle=" at 135W"
    elif section == "122W":
        lat=None
        lon=-122.0
        xmin=-10.0
        xmax=10.0
        sectitle=" at 122W"

    if cmap is None:
        cmap="RdYlBu_r"

    if utrd is not None:
        title=utrd+" (m/s2 - colours) and mean eastward current (m/s - lines) "+sectitle
        if 'div' in utrd:
            # in this case we are plot -1 * divergence:
            title="-"+title
        outfile=filestem+"_"+utrd+"_uo_"+section+".png"
    else:
        title="mean eastward current (m/s) "+sectitle
        outfile=filestem+"_uo_"+section+".png"

    plot_nemo_section(filenames=filenames,var_names=var_names,toldeg=toldeg,
                      lon=lon,lat=lat,xmin=xmin,xmax=xmax,depthmax=300,cmap=cmap,plot_types=["c","l"],
                      mnfld=mnfld,mxfld=mxfld,factor=factor,scientific=scientific,
                      title=title, outfile=outfile)

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-I", "--filestem",action="store",help="name stem of input files")
    parser.add_argument("-i", "--filenames",action="store",nargs="+",help="name(s) of input file(s)")
    parser.add_argument("-u", "--utrd",action="store",dest="utrd",help="name of trend field") 
    parser.add_argument("-s", "--section",action="store",dest="section",help="section to plot: eqpac, 135W, 122W")
    parser.add_argument("-c", "--cmap",action="store",dest="cmap",help="colour map (default RdYlBu_r)")
    parser.add_argument("-f", "--mnfld",action="store",dest="mnfld",type=float,nargs="+",
                        help="minimum field value(s) to plot")
    parser.add_argument("-F", "--mxfld",action="store",dest="mxfld",type=float,nargs="+",
                        help="maximum field value(s) to plot")
    parser.add_argument("-Z", "--factor", action="store",dest="factor",type=float,nargs="+",
                    help="multiplicative factor(s) to apply to plotted field(s)")
    parser.add_argument("-Q", "--toldeg", action="store",dest="toldeg",type=float,
                    help="tolerance (degrees) for finding section")

    args = parser.parse_args()

    plot_utrd_section(filestem=args.filestem,filenames=args.filenames,utrd=args.utrd,
                      section=args.section,cmap=args.cmap,
                      mnfld=args.mnfld,mxfld=args.mxfld,factor=args.factor,toldeg=args.toldeg)
