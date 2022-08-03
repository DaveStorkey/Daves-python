#! /usr/bin/env python
"""
    Script to plot a profile read in from a netcdf file.

Mar 2020 : A few updates including update to python 3. DS.
Apr 2020 : Incorporated runids keyword and style.db from one of Pierre's scripts. DS.
"""
__author__ = 'Dave Storkey (dave.storkey@metoffice.gov.uk)'
__date__ = '$Date: 2013/04/29 $'

import sys
import numpy as np
import numpy.ma as ma
import matplotlib
import matplotlib.pyplot as plt
import netCDF4

# ============================ class run =====================================
class run(object):
    def __init__(self, runid):
        # parse dbfile
        self.runid, self.name, self.line, self.color = parse_dbfile(runid)

# ============================ file parser =====================================
def parse_dbfile(runid):
    try:
        lstyle=False
        with open('style.db') as fid:
            for cline in fid:
                att=cline.split('|')
                if att[0].strip() == runid:
                    cpltrunid = att[0].strip()
                    cpltname  = att[1].strip()
                    cpltline  = att[2].strip()
                    cpltcolor = att[3].strip()
                    lstyle=True
        if not lstyle:
            raise Exception(runid+' not found in style.db')

    except Exception as e:
        print ('Issue with file : style.db')
        print (e)
        sys.exit(42)

    # return value
    return cpltrunid, cpltname, cpltline, cpltcolor

# ============================ main routine =====================================
def plot_profile(datafiles,fields,maxdepth=None,xmin=None,xmax=None,xtitle=None,legend=None,legendloc="upper right",
                 linestyles=None,colors=None,outfile=None,normfield=None,normfiles=None,runids=None,kup=None):

    if normfield is not None:
        if normfiles is not None and len(normfiles) == 1:
            # read normfield once for all if only one normfile specified
            with netCDF4.Dataset(normfiles[0], mode='r') as dataid:
                            normfld = np.squeeze(np.copy(dataid.variables[normfield]))
        else:
            # if normfile not specified try to read normfield from first file:
            with netCDF4.Dataset(datafiles[0], mode='r') as dataid:
                            normfld = np.squeeze(np.copy(dataid.variables[normfield]))
    elif normfiles is not None:
        normfield = fields[0]
        if len(fields) != 1:
            raise Exception("Error: If specifying normfile (and not normfield) then all the files must contain the same field (len(fields)==1).")
        elif len(normfiles) == 1:
            with netCDF4.Dataset(normfiles[0], mode='r') as dataid:
                            normfld = np.squeeze(np.copy(dataid.variables[normfield]))
    
    if len(datafiles) != len(fields):
        if len(datafiles) == 1:
            for i in range(len(fields)-1):
                datafiles.append(datafiles[0])
        if len(fields) == 1:
            for i in range(len(datafiles)-1):
                fields.append(fields[0])
            
    if runids is not None:
        if len(runids) == len(datafiles):
            run_lst=[]
            for runid in runids:
                 # initialise run
                 run_lst.append(run(runid))
            colors = [r.color for r in run_lst]
            linestyles = [r.line for r in run_lst]
            legend = [r.name for r in run_lst]
        else:
            raise Exception('Error: number of runids must equal the number of files.')

    else:
        # runids keyword over-rides colors, linestyles and legend keywords.
        if colors is None:
            colors=len(fields)*[None]
        else:
            if len(colors) != len(fields):
                if len(colors) == 1:
                    colors=len(fields)*colors
                else:
                    raise Exception('Error: number of elements of colors must be one or equal to number of fields.')

        if linestyles is None:
            linestyles=len(fields)*[None]
        else:
            if len(linestyles) != len(fields):
                if len(linestyles) == 1:
                    linestyles=len(fields)*linestyles
                else:
                    raise Exception('Error: number of elements of linestyles must be one or equal to number of fields.')

    count = -1
    for [filename,field,color,linestyle] in zip(datafiles,fields,colors,linestyles):

        count += 1
        dataid = netCDF4.Dataset(filename, mode='r')
        for depthname in ['depth', 'deptht', 'depthu', 'depthv', 'depthw']:
            try:
                depths = np.copy(dataid.variables[depthname])
            except KeyError:
                pass
            else:
                break
        else:
            raise KeyError

        fld = np.squeeze(np.copy(dataid.variables[field]))
        if normfiles is not None or normfield is not None:
            if len(normfiles) > 1:
                with netCDF4.Dataset(normfiles[count], mode='r') as dataid:
                            normfld = np.squeeze(np.copy(dataid.variables[normfield]))
            small = np.min(np.absolute(normfld[np.where(normfld != 0.0)]))*0.001
            print('small is : ',small)
            fld[:] = fld[:]/(normfld[:]+small)

        if maxdepth is not None:
            if kup is None:
                depths_select=depths[np.where(depths < maxdepth)]
                fld_select=fld[0:len(depths_select)]
            else:
                # if kup is set then "maxdepth" means max value of kup to plot
                kup_select=np.arange(maxdepth)
                fld_select=fld[0:len(kup_select)]
        else:
            if kup is None:
                depths_select=depths
            else:
                kup_select=np.arange(len(fld))
            fld_select=fld

        if count == 0:
            if xmin is not None:
                plt.gca().set_xlim(left=xmin)
            if xmax is not None:
                plt.gca().set_xlim(right=xmax)
            if xtitle is not None:
                plt.xlabel(xtitle)
            if kup:
                plt.ylabel('kup')
            else:
                plt.ylabel('depths (m)')

        if kup:
            plt.plot(fld_select,kup_select,color=color,linestyle=linestyle)
        else:
            plt.plot(fld_select,-depths_select,color=color,linestyle=linestyle)

    print('legendloc is ',legendloc)
    if legend is not None:
        plt.legend(legend,loc=legendloc)

    if outfile is not None:
        matplotlib.rcParams['font.size'] =8
        plt.savefig(outfile,dpi=200)
    else:
        plt.show()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--files",action="store",dest="datafiles", help="names of input files",nargs="+")
    parser.add_argument("-n","--fields",action="store",dest="fields", help="names of fields to plot",nargs="+")
    parser.add_argument("-F", "--normfiles", action="store",dest="normfiles",default=None,nargs="+",
                    help="file(s) containing profile to normalise all the others by.")
    parser.add_argument("-N", "--normfield", action="store",dest="normfield",default=None,
                    help="name of profile(s) to normalise all the others by.")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",default=None,
                    help="plot to output file - filetype by extension. If unset plot to GUI.")
    parser.add_argument("-K", "--kup", action="store_true",dest="kup",default=None,
                    help="plot against kup (model level index above bathymetry) rather than depth")
    parser.add_argument("-D", "--maxdepth", action="store",dest="maxdepth",type=int,default=None,
                    help="maximum depth to plot")
    parser.add_argument("-x", "--xmin", action="store",dest="xmin",type=float,default=None,
                    help="start of x-axis")
    parser.add_argument("-X", "--xmax", action="store",dest="xmax",type=float,default=None,
                    help="end of x-axis")
    parser.add_argument("-c", "--colors", action="store",dest="colors",default=None,nargs="+",
                    help="color(s)")
    parser.add_argument("-r", "--runids", action="store",dest="runids",default=None,nargs="+",
                    help="runid(s): if specified read descriptions, line colours and line styles from style.db")
    parser.add_argument("-s", "--linestyles", action="store",dest="linestyles",default=None,nargs="+",
                    help="linestyle(s)")
    parser.add_argument("-t", "--xtitle", action="store",dest="xtitle",default=None,
                    help="title of x-axis")
    parser.add_argument("-L", "--legend", action="store",dest="legend",default=None,nargs="+",
                    help="list of labels for legend")
    parser.add_argument("-l", "--legendloc", action="store",dest="legendloc",default="upper right",
                    help="location of legend on plot")
    args = parser.parse_args()

    plot_profile(args.datafiles,args.fields,maxdepth=args.maxdepth,xmin=args.xmin,xmax=args.xmax,xtitle=args.xtitle,
                 linestyles=args.linestyles,colors=args.colors,legend=args.legend,legendloc=args.legendloc,
                 outfile=args.outfile,normfiles=args.normfiles,normfield=args.normfield,runids=args.runids,kup=args.kup)        
