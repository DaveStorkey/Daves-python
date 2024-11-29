#! /usr/bin/env python
# 
# Routine to plot histogram of an ocean field.
#
# Aug 2017  :  Add option to do a proper volume census as opposed to just
#              counting cells. NB. hardwired to read in the T-grid cell dimensions
#              at the moment. DS.
#
# Nov 2018  :  Add option to change grid type (T,U,V,F) for calculating cell 
#              volumes for volume census. DS.
#
# June 2020 :  1. Convert to python 3
#              2. Add area_census option for 2D fields. 
#              3. Generalise so you can plot more than one dataset. 
#              4. Add legend. 
#              5. Add option for scientific format of xtick labels.
#              DS.
#
# Nov 2024  : Rename plot_histogram_multi.py -> plot_histogram_coords.py 
#             and put in a couple of fixes/tweaks. DS.

import netCDF4 as nc
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import textwrap

def plot_histogram(filenames=None, varnames=None, outfile=None, coordsfilenames=None, 
                   logy=False, nbins=None, grid=None, xmin=None, xmax=None, ymin=None, ymax=None, alpha=None,
                   rec=None, printvalues=False,fieldlabel=None,legend=None,scientific=None,title=None):

    if len(varnames) != len(filenames):
        if len(filenames) == 1:
            filenames = filenames*len(varnames)
        elif len(varnames) == 1:
            varnames = varnames*len(filenames)
        else:
            raise Exception('Error: number of varnames must equal number of filenames (or one of those must have length 1).')

    if coordsfilenames is None:
        coordsfilenames = [None]*len(filenames)
    elif len(coordsfilenames) != len(filenames):
        if len(coordsfilenames) == 1:
            coordsfilenames = coordsfilenames*len(filenames)
        else:
            raise Exception('Error: number of coordsfilenames must equal number of filenames (or one of those must have length 1).')

    if nbins is None:
        nbins=50

    for coordsfilename, filename, varname in zip(coordsfilenames,filenames,varnames):

        print(' ')
        print(coordsfilename, filename, varname)
        print(' ')

        infile=nc.Dataset(filename,'r')
        if rec is not None:
            infield=infile.variables[varname][rec]
        else:
            infield=infile.variables[varname][0]

        # exclude mask points if there is a mask
        try:
            infield_valid_points = infield[~infield.mask]
        except AttributeError:
            infield_valid_points = infield

        area_census = False
        volume_census = False
        if coordsfilename is not None:
            infield_coords=infile.variables[varname].coordinates
            if grid is None:
                print('Grid unspecified - defaulting to T-grid cell volumes')
                grid='t'
            elif grid not in 'tTuUvVfF':
                raise Exception('Must specify T, U, V or F grid.')
            else:
                grid = grid.lower()
            coordsfile = nc.Dataset(coordsfilename)
            # [0] index gets rid of degenerate time dimension:
            e1 = coordsfile.variables['e1t'][0][:]
            e2 = coordsfile.variables['e2t'][0][:]
            if 'dep' in infield_coords:
                volume_census = True
                e3coords=False
                try:
                    e3_var = infile.variables['e3'+grid]
                except KeyError:
                    try:
                        e3_var = coordsfile.variables['e3'+grid+'_0']
                    except KeyError:
                        raise Exception('Could not find e3'+grid+' field in infile or coordsfile.')
                    else:
                        e3coords=True
                if e3coords:
                    e3 = e3_var[0][:]
                elif rec is not None:
                    e3 = e3_var[rec][:]
                else:
                    e3 = e3_var[:]
                # convert from m3 to km3
                cell_volume = e1[:] * e2[:] * e3[:] * 1.e-9
            else:
                area_census = True
                # area in km2:
                cell_area = e1[:] * e2[:] * 1.e-6

        print(' ')
        print('area_census : ',area_census,'volume_census : ',volume_census)
        print(' ')

        if area_census:
            try:
                area_valid_points = cell_area[~infield.mask]
            except AttributeError:
                area_valid_points = cell_area

        if volume_census:
            try:
                volume_valid_points = cell_volume[~infield.mask]
            except AttributeError:
                volume_valid_points = cell_volume

        # use xmin and xmax (if set) to define histogram bins, not just the limits of the plot
        if xmin is not None:
            infield_valid_points = infield_valid_points[np.where(infield_valid_points > xmin)]
            if volume_census:
                volume_valid_points = volume_valid_points[np.where(infield_valid_points > xmin)]
            elif area_census:
                area_valid_points = area_valid_points[np.where(infield_valid_points > xmin)]
        else:
            # in this case set xmin to be the minimum of the *first* field
            xmin = np.min(infield_valid_points)
        if xmax is not None:
            infield_valid_points = infield_valid_points[np.where(infield_valid_points < xmax)]
            if volume_census:
                volume_valid_points = volume_valid_points[np.where(infield_valid_points < xmax)]
            elif area_census:
                area_valid_points = area_valid_points[np.where(infield_valid_points < xmax)]
        else:
            # in this case set xmax to be the maximum of the *first* field
            xmax = np.max(infield_valid_points)

        if volume_census:
            hist, bins = np.histogram(infield_valid_points, bins=nbins, range=(xmin,xmax), weights=volume_valid_points)
        elif area_census:
            hist, bins = np.histogram(infield_valid_points, bins=nbins, range=(xmin,xmax), weights=area_valid_points)
        else:
            hist, bins = np.histogram(infield_valid_points, bins=nbins, range=(xmin,xmax))
        if printvalues:
            for (b1,b2,h) in zip(bins[0:-1],bins[1:],hist):
                print(b1,b2,h)
            if area_census:
                print('total area : ',np.sum(hist),' km2')
            elif volume_census:
                print('total volume : ',np.sum(hist),' km3')
        width = 0.8 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        plt.bar(center, hist, align='center', width=width, alpha=alpha)

    if xmin is None:
        xmin = np.amin(infield_valid_points)
    if xmin is None:
        xmax = np.amax(infield_valid_points)
    if ymin is None and ymax is not None:
        ymin = 0.0
    if ymax is None and ymin is not None:
        ymax = np.amax(hist)
    plt.gca().set_xlim([xmin,xmax])
    if ymin is not None and ymax is not None:
        plt.gca().set_ylim([ymin,ymax])

    if logy:
        plt.yscale('log')
    if fieldlabel is not None:
        plt.gca().set_xlabel(fieldlabel)
    else:
        plt.gca().set_xlabel(varname)
    if area_census:
        plt.gca().set_ylabel('area (km2)')
    elif volume_census:
        plt.gca().set_ylabel('volume (km3)')
    else:
        plt.gca().set_ylabel('number of points')

    if scientific:
        plt.gca().ticklabel_format(axis='x',style='sci',scilimits=(0,0))

    if legend:
        plt.gca().legend(legend)

    if title is not None:
        title=textwrap.fill(title,70)
        plt.gcf().suptitle(title, fontsize=12, y=0.95)    

    if outfile:
        matplotlib.rcParams['font.size'] =8
        plt.savefig(outfile,dpi=200)
    else:
        plt.show()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--filenames", action="store",dest="filenames",nargs="+",
                    help="name(s) of input file(s)")
    parser.add_argument("-v", "--varnames", action="store",dest="varnames",nargs="+",
                    help="name(s) of variable(s) to plot")
    parser.add_argument("-r", "--rec", action="store",dest="rec",type=int,default=None,
                    help="record to read from file (counting from zero)")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",
                    help="name of output plot - if unset will plot to GUI")
    parser.add_argument("-l", "--fieldlabel", action="store",dest="fieldlabel",
                    help="name of field for x-axis label")
    parser.add_argument("-c", "--coordsfilenames", action="store",dest="coordsfilenames",nargs="+",
                    help="name(s) of coordinates or meshmash file(s) - if set will plot an area or volume census of input field rather than counting cells.")
    parser.add_argument("-g", "--grid", action="store",dest="grid",
                    help="type of grid: T, U, V, F - determines which cell dimensions to use for volume census.")
    parser.add_argument("-x", "--xmin", action="store",dest="xmin",type=float,default=None,
                    help="minimum value of x-axis (also used to define bins)")
    parser.add_argument("-X", "--xmax", action="store",dest="xmax",type=float,default=None,
                    help="maximum value of x-axis (also used to define bins)")
    parser.add_argument("-y", "--ymin", action="store",dest="ymin",type=float,default=None,
                    help="minimum value of y-axis")
    parser.add_argument("-Y", "--ymax", action="store",dest="ymax",type=float,default=None,
                    help="maximum value of y-axis")
    parser.add_argument("-n", "--nbins", action="store",dest="nbins",type=int,default=None,
                    help="number of bins for histogram (default 50)")
    parser.add_argument("-G", "--logy", action="store_true",dest="logy",
                    help="log scale for y-axis")
    parser.add_argument("-L", "--legend", action="store",dest="legend",nargs="+",
                    help="legend entries")
    parser.add_argument("-e", "--scientific", action="store_true",dest="scientific",
                    help="scientific format for x-axis labels")
    parser.add_argument("-p", "--print", action="store_true",dest="printvalues",
                    help="print out bins and values of the histogram as well as plotting ")
    parser.add_argument("-t", "--title", action="store",dest="title",
                    help="title for plot")
    parser.add_argument("-a", "--alpha", action="store",dest="alpha",type=float,
                    help="transparency of bars")
 
    args = parser.parse_args()

    plot_histogram(filenames=args.filenames, varnames=args.varnames, outfile=args.outfile, 
                   coordsfilenames=args.coordsfilenames, grid=args.grid,legend=args.legend,
                   logy=args.logy, rec=args.rec, nbins=args.nbins, fieldlabel=args.fieldlabel,
                   xmin=args.xmin, xmax=args.xmax, ymin=args.ymin, ymax=args.ymax, alpha=args.alpha,
                   scientific=args.scientific, printvalues=args.printvalues, title=args.title)
