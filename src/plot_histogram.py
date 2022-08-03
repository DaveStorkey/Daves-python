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
# June 2020 :  Convert to python 3. DS.
#
# Feb 2021  :  Adapt so can be passed an array from another python routine. DS.

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma

def plot_histogram(filename=None, var=None, outfile=None, coordsfilename=None, cell_area=None, cell_volume=None, 
                   logy=False, nbins=None, grid=None, xmin=None, xmax=None, ymin=None, ymax=None, rec=None, 
                   printvalues=False, fieldlabel=None, noshow=None):

    if nbins is None:
        nbins=50

    if filename is not None and isinstance(var,str):
        infile=nc.Dataset(filename,'r')
        if rec is not None:
            infield=infile.variables[var][rec]
        else:
            infield=infile.variables[var][0]
        infield_coords=infile.variables[var].coordinates
    elif isinstance(var,np.ndarray):
        infield=var
    else:
        raise Exception("Error : either specify input file and variable name or pass in numpy array")

    area_census = False
    volume_census = False
    if cell_volume is not None:
        volume_census = True
    elif cell_area is not None:
        area_census = True
    elif coordsfilename is not None:
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

    if volume_census:
        print("Doing volume census.")
    elif area_census:
        print("Doing area census")

    # exclude mask points if there is a mask
    try:
        infield_valid_points = infield[~infield.mask]
    except AttributeError:
        print("No mask for field")
        infield_valid_points = infield

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
    if xmax is not None:
        infield_valid_points = infield_valid_points[np.where(infield_valid_points < xmax)]
        if volume_census:
            volume_valid_points = volume_valid_points[np.where(infield_valid_points < xmax)]

    if volume_census:
        hist, bins = np.histogram(infield_valid_points, bins=nbins, weights=volume_valid_points)
    elif area_census:
        hist, bins = np.histogram(infield_valid_points, bins=nbins, weights=area_valid_points)
    else:
        hist, bins = np.histogram(infield_valid_points, bins=nbins)
    if printvalues:
        for (b1,b2,h) in zip(bins[0:-1],bins[1:],hist):
            print(b1,b2,h)
    width = 0.8 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=width)

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
    elif isinstance(var,str):
        plt.gca().set_xlabel(var)
    if area_census:
        plt.gca().set_ylabel('area (km2)')
    elif volume_census:
        plt.gca().set_ylabel('volume (km3)')
    else:
        plt.gca().set_ylabel('number of points')
    if outfile:
        plt.savefig(outfile,dpi=200)
    elif not noshow:
        plt.show()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--filename", action="store", dest="filename", help="name of input file")
    parser.add_argument("-v", "--varname", action="store", dest="varname", help="name of field to plot")
    parser.add_argument("-r", "--rec", action="store",dest="rec",type=int,default=None,
                    help="record to read from file (counting from zero)")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",
                    help="name of output plot - if unset will plot to GUI")
    parser.add_argument("-f", "--fieldlabel", action="store",dest="fieldlabel",
                    help="name of field for x-axis label")
    parser.add_argument("-c", "--coordsfilename", action="store",dest="coordsfilename",
                    help="name of coordinates or meshmash file - if set will plot a volume census of input field rather than counting cells.")
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
    parser.add_argument("-L", "--logy", action="store_true",dest="logy",default=False,
                    help="log scale for y-axis")
    parser.add_argument("-p", "--print", action="store_true",dest="printvalues",default=False,
                    help="print out bins and values of the histogram as well as plotting ")
 
    args = parser.parse_args()

    plot_histogram(filename=args.filename, var=args.varname, outfile=args.outfile, coordsfilename=args.coordsfilename, grid=args.grid,
                   logy=args.logy, rec=args.rec, nbins=args.nbins, fieldlabel=args.fieldlabel,
                   xmin=args.xmin, xmax=args.xmax, ymin=args.ymin, ymax=args.ymax, printvalues=args.printvalues)
