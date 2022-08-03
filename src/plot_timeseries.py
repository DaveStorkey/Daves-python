#! /usr/bin/env python2.7

'''
Created on Apr 3, 2013

@author: hadcv
'''

# =============================
# Imports
# =============================
from datetime import datetime
from numpy import argmin, invert
from numpy.ma import mask_or, getmask
import argparse
import numpy
import scipy.stats
import os
import sys
import warnings

from Iris_local_utils import errlev, cm_custom, latlon_index, date_ticks, cl_cat_time, slice_time
from NETCDF_utils import str2time, nc2dt
from iris.analysis import MEAN
from itertools import cycle
import iris

import Iris_Tim_tools as ttools
#----------------------------------------------------------


# =============================
# Argument Parsing
# =============================
args = argparse.ArgumentParser(description='Plots 1D and 2D timeseries')

# Positional
args.add_argument('file_ts1', help='Filenames containing timeseries data. '
                  'Separate timeseries should be specified as separate arguments, while files to be '
                  'concatenated into a single timeseries should be given as glob strings.',
                  metavar='timeseries', type=str, nargs='*')
args.add_argument(
    'name_var', help='Name of the variable to be plotted', metavar='variable', type=str)

# Optional
args.add_argument('-r', help='Filename containing the timeseries against which to compare',
                  metavar='timeseries_ref', type=str, dest='file_ts2', nargs=1)
args.add_argument('-l', help='Use logarithmic contouring for 2D plots',
                  dest='plot_log', action='store_true')
# Example: [ 'a = %i, b = %i' % (i,j) for (i,j) in ((1,2),(3,4),(5,6)) ]
args.add_argument('--labels', type=str, help='Legend labels for 1D plots', nargs='*', dest='plot_labels',
                  metavar='custom_labels')
args.add_argument('--title', type=str, help='Title for the diagram', nargs=1, dest='plot_title',
                  metavar='custom_title')
args.add_argument('--xaxis', type=str, help='X axis label for the diagram', nargs=1, dest='plot_xaxis',
                  metavar='x-axis_label')
args.add_argument('--xaxis_type', type=str, help='Format of the x axis ticks', nargs=1,
                  metavar='x-axis_tick_type', choices=['date', 'none'], default=['date'])
args.add_argument('--yaxis', type=str, help='Y axis label for the diagram', nargs=1, dest='plot_yaxis',
                  metavar='y-axis_label')
args.add_argument('--line_styles', type=str, help='Matplotlib line styles to use', nargs='*', metavar='line_styles')
args.add_argument('--line_widths', type=float, help='Line widths to use', nargs='*', metavar='line_widths')
args.add_argument('--line_colours', help='Line colours to use (any matplotlib specification)', nargs='*', metavar='line_colours')
args.add_argument('--figsize', help='Figure size (width, height) in inches', nargs=2, type=float, default=[8, 5])

args.add_argument('--diff', help='1D plots are produced as differences',
                  action='store_true', dest='plot_1d_diff')
args.add_argument('--diff_origin', help='2D plots are produced as differences to the first time coordinate (after bounding)',
                  action='store_true', dest='plot_diff_origin')
args.add_argument('-e', type=float, help='Fraction of data to be contoured as extrema',
                  dest='plot_ext', metavar='fraction_of_contour_extrema', default=.01)
# args.add_argument('--lat',type=float,help='Latitude range to slice',
#                                dest='coord_lat',metavar='latitude_range',nargs='*')
# args.add_argument('--lon',type=float,help='Longitude range to slice',
#                                dest='coord_lon',metavar='longitude_range',nargs='*')
args.add_argument('--loci', type=float, help='Latitude/Longitude of target location, either [lat,lon] specifying a point '
                  'or [lat_l,lon_l,lat_u,lon_u] specifying the bottom-left and top-right corners of a box',
                  dest='coord_latlon', metavar='lat_lon_loci', nargs='*')
args.add_argument('--mask', type=float, help='Latitude/Longitude of target location to mask; [lat_l,lon_l,lat_u,lon_u] specifying '
                  'the bottom-left and top-right corners of a box',
                  dest='coord_latlon_mask', metavar='lat_lon_mask', nargs=4)

args.add_argument('-z', type=float, help='Depth range to slice',
                  dest='coord_depth', metavar='depth_range', nargs='*')
args.add_argument('--z_oper', action='store', nargs=1, type=str, help='Operation with which to collapse depth coordinate',
                  dest='z_oper', choices=['max', 'mean', 'median', 'min', 'rms', 'std_dev', 'sum', 'variance'])
args.add_argument('--z_wgt', action='store', nargs=2, type=str, default=None,
                        help='The filename and variable name to be used for weighting the depth coordinate collapse operation',
                        dest='z_wgt')
args.add_argument('-n', type=int, help='Number of contours to use for 2D plots', nargs=1,
                  dest='plot_nlev', metavar='no_contours', default=[20])
args.add_argument('--from', type=str, help='Lower bound of time range to slice of form yyyy[/mm[/dd]]',
                  dest='coord_time_l', metavar='time_range_lower')
args.add_argument('--to', type=str, help='Upper bound of time range to slice of form yyyy[/mm[/dd]]',
                  dest='coord_time_u', metavar='time_range_upper')
args.add_argument('-y', type=float, help='Limits on the y axis', nargs=2,
                  dest='plot_lims', metavar='yaxis_limits', default=None)
args.add_argument('-d', type=float, help='Limits on the contours to draw for 2D plots', nargs=2,
                  dest='plot_c_lims', metavar='contour_limits', default=None)
args.add_argument(
    '-f', type=str, help='Filename to save as', dest='filename', nargs=1)

input = args.parse_args()

import matplotlib

if input.filename is not None:
    matplotlib.use('Agg')

# Matplotlib imports (after setting backend)
from matplotlib.colors import BoundaryNorm
from matplotlib.dates import DateFormatter, date2num, YearLocator, MonthLocator, DayLocator, AutoDateLocator
from matplotlib.pyplot import xlim, ylim, cm, colorbar, legend, title
from matplotlib.ticker import MaxNLocator, LogLocator, AutoMinorLocator
import matplotlib.pyplot as mplt

if input.plot_log:
    plot_locator = LogLocator(input.plot_nlev[0])
else:
    plot_locator = MaxNLocator(input.plot_nlev[0])

i_count = 0
cm_custom()

if input.file_ts2 is not None:
    file_list = input.file_ts2
    file_list.extend(input.file_ts1)
else:
    file_list = input.file_ts1


# =============================
# Loop over files in arguments
# =============================
for i in file_list:              # Produce a cube for each file

    i_count += 1

    # =============================
    # Retrieve variable cube using NEMO name
    # =============================
    cube_ts = iris.load(i, iris.Constraint(
        cube_func=lambda nemoname: nemoname.var_name == input.name_var))

    if len(cube_ts) == 0:
        raise StandardError(
            'No cubes exist in file %s with the variable %s' % (i, input.name_var))
    elif len(cube_ts) > 1:
        # BUG: The biggus operations within result in inaccessable data
        cube_ts = cl_cat_time(cube_ts)
    else:
        cube_ts = cube_ts[0]

    # =============================
    # Data masking
    # =============================
    if input.coord_latlon_mask is not None and i_count == 1:
        mask = ttools.mask_gen(cube_ts.coord('longitude').points, cube_ts.coord('latitude').points,
                               [input.coord_latlon_mask[i] for i in [2, 3, 0, 1]])
        cube_ts.data.mask = mask_or(invert(mask), getmask(cube_ts.data))

    # =============================
    # Reduce cube (multidimensional lat/lon)
    # =============================

    # REDUCE cube by depth
    if input.coord_depth is not None:
        try:
            if len(input.coord_depth) == 2:
                depth_const = iris.Constraint(depth = lambda zra: input.coord_depth[0] <= zra <= input.coord_depth[1])
            elif len(input.coord_depth) == 1:
                zra_minloc = argmin(
                    abs(cube_ts.coord('depth').points - input.coord_depth[0]))
                depth_const = iris.Constraint(depth=cube_ts.coord('depth').points[zra_minloc])

            cube_ts = cube_ts.extract(depth_const)

        except iris.exceptions.CoordinateNotFoundError:
            pass

    # REDUCE cube by time
    if (input.coord_time_l is not None) or (input.coord_time_u is not None):
        cube_ts = slice_time(cube_ts, [input.coord_time_l, input.coord_time_u])

    # Indices for target lat/lon
    if input.coord_latlon is not None:
        if len(input.coord_latlon) == 2:          	   # Point
            print 'Calculating point coordinate for cube %i' % i_count
            i_latlon = latlon_index(cube_ts, input.coord_latlon)
        elif len(input.coord_latlon) == 4:        	   # Box
            print 'Calculating bottom-left and top-right coordinates for cube %i' % i_count
            i_latlon = latlon_index(cube_ts, input.coord_latlon[0:2]) + \
                latlon_index(cube_ts, input.coord_latlon[2:4])
        else:
            raise StandardError('ERROR: invalid argument for --loci argument')
    elif cube_ts.coord('longitude').points.size == 3:  # 3x3 C1D workaround
        print 'Assuming a 3x3 C1D grid -- lat/lon indices [1,1] are used'
        i_latlon = [1, 1]
    elif cube_ts.coord('longitude').points.size == 1:  # 1x1 C1D workaround
        i_latlon = [0, 0]
    else:
        # No location specified: global domain
        i_latlon = None

    # COLLAPSE cube by lat/lon (a point or a box based on
    # bottom-left/top_right corners)
    if i_latlon is not None:
        if cube_ts.ndim == 4:   # T,Z,Y,X
            if len(i_latlon) == 2:
                cube_ts = cube_ts[:, :, i_latlon[0], i_latlon[1]]
            else:
                # Problem with wraparound indexing!
                cube_ts = cube_ts[
                    :, :, i_latlon[0]:i_latlon[2] + 1, i_latlon[1]:i_latlon[3] + 1]

        elif cube_ts.ndim == 3:  # T,Y,X
            if len(i_latlon) == 2:
                cube_ts = cube_ts[:, i_latlon[0], i_latlon[1]]
            else:
                # And here..
                cube_ts = cube_ts[
                    :, i_latlon[0]:i_latlon[2] + 1, i_latlon[1]:i_latlon[3] + 1]

    try:
        # TODO: This requires area weighting!
        cube_ts = cube_ts.collapsed(['latitude', 'longitude'], MEAN)
    except iris.exceptions.CoordinateCollapseError:
        pass

    # COLLAPSE cube by depth
    if input.z_oper is not None:
        
        collapse_fn = getattr(iris.analysis, input.z_oper[0].upper())
                    
        # Do not process single-length depth dimensions
        if 'depth' in [j.name() for j in cube_ts.dim_coords]:

            # Do not use weighting if the operator does not support it
            if 'uses_weighting' not in dir(collapse_fn):
                meshfile = None
                        
                if input.z_wgt is not None:
                    print 'Depth operator "%s" does not use weighting; --z_wgt argument is ignored' \
                    % input.z_oper[0]
                        
            # Meshmask file
            elif input.z_wgt is not None:
                meshfile = input.z_wgt[0]
            else:
                meshfile = '/data/cr1/hadcv/mesh_masks/ORCA1.nc'

            # Process weighting variable
            if meshfile is not None:

                # Load weighting variable cube
                if input.z_wgt is not None:
                    e3t = iris.load_cube(meshfile, iris.Constraint(
                        cube_func=lambda nemoname: nemoname.var_name == input.z_wgt[1]))
                else:
                    e3t = iris.load_cube(meshfile, iris.Constraint(
                        cube_func=lambda nemoname: nemoname.var_name == 'e3t_0'))

                # Remove time coordinate from weighting variable if needed
                if 'time' in [i.name() for i in e3t.dim_coords]:
                    if len(e3t.coord('time').points) > 1:
                        print 'A time-varying weighting variable was found (%s); the first time is used' \
                        % e3t.__repr__()
                                
                        e3t = e3t[sliceby(e3t,'time',0)]
                            
                # Check that at this point we have a 1D weighting variable
                if e3t.ndim != 1:
                    print 'Weighting variable has %s dimensions; currently only 1D variables are supported. ' \
                    'Weighting is skipped' % e3t.ndim
                    meshfile = None
                        
                # Extract depths to match variable being weighted if required
                if input.coord_depth is not None:     
                    if 'depth' in [i.name() for i in e3t.dim_coords]:
                        e3t = e3t.extract(depth_const)
                    else:
                        print 'No dimension coordinate with standard_name="depth" in cube %s; ' \
                        'could not extract required depths' % e3t.__repr__()

                # Broadcast the 1D weighting variable to the shape of the variable being weighted
                if e3t.shape[0] == cube_ts.shape[1]:
                    e3t_tmp = cube_ts.copy()
                    tile_shape = numpy.asarray(e3t_tmp.shape)
                    tile_shape = tile_shape[tile_shape != e3t.shape].tolist() + [1]
                            
                    try:
                        e3t_tmp.data = numpy.tile(e3t.data,tile_shape)
                        
                        # Check broadcast was done correctly
                        for j in xrange(e3t.shape[0]):
                            assert (e3t_tmp[:,j].data == e3t[j].data).all()
                        
                        e3t = e3t_tmp 

                    except:
                        print 'Could not make 1D weighting variable (%s) conform to global domain' \
                        % e3t.__repr__()

                # If weighting variable still does not conform in shape, skip weighting
                if e3t.shape != cube_ts.shape:
                    print 'Weighting cube dimensions %s do not match those of cube ' \
                          'to be collapsed %s. Weighting is skipped' % (e3t.shape, cube_ts.shape)
                    meshfile = None

                # Collapse depth dimension
                if meshfile is not None:
                    cube_ts = cube_ts.collapsed('depth', collapse_fn, weights=e3t.data)
                else:
                    cube_ts = cube_ts.collapsed('depth', collapse_fn)


    # Convert cube time units to that of the reference cube
    if (input.plot_1d_diff or cube_ts.ndim == 2 or input.xaxis_type[0] == 'none') and (i_count != 1):
        cube_ts.coord('time').convert_units(cube_out[0].coord('time').units)

    # Check whether the field has both positive and negative values
    if i_count == 1:
        if (cube_ts.data.max() > 0 and cube_ts.data.min() < 0) or \
           (input.plot_diff_origin and cube_ts.ndim == 2):
            plot_plus_minus = True
        else:
            plot_plus_minus = False

    # Difference cubes 
    if cube_ts.ndim == 1:

        if i_count == 1:
            cube_out = [cube_ts]
        else:
            cube_out += [cube_ts]

        if (i_count == len(file_list)) and input.plot_1d_diff:
            cube_out = [i - cube_out[0] for i in cube_out]

    elif cube_ts.ndim == 2:

        if i_count == 1:
            cube_out = cube_ts

            if plot_plus_minus:
                plot_colour = cm.get_cmap('idl_bgwor')
            else:
                plot_colour = cm.get_cmap('gist_ncar')

            plot_sym = False

            if input.plot_diff_origin:
                cube_out -= cube_out[0]
                break

        elif i_count == 2:
            cube_out = cube_ts - cube_out
            plot_colour = cm.get_cmap('idl_bgwor')
            plot_sym = True

        else:
            # 2D cube differencing does not currently
            # make sense of more than 2 files being
            # supplied
            break


# =============================
# Plot
# =============================

# Some figure parameters..
fig = mplt.figure(figsize=input.figsize)

if input.plot_title is None:
    if cube_ts.ndim == 1:
        plot_title = 'Timeseries plot of variable %s' % input.name_var
    else:
        plot_title = 'Hovmuller timeseries plot of variable %s' % input.name_var
else:
    plot_title = input.plot_title[0]

if input.plot_xaxis is None:
    plot_xaxis = 'Time'
else:
    plot_xaxis = input.plot_xaxis[0]

if input.plot_yaxis is None:
    if cube_ts.ndim == 1:
        plot_yaxis = '%s (%s)' % (input.name_var, cube_ts.units)
    else:
        plot_yaxis = 'Depth'
else:
    plot_yaxis = input.plot_yaxis[0]

# Plotting commands
if cube_ts.ndim == 2:
    cube_out.transpose()

    # x-axis uses dates
    if input.xaxis_type[0] == 'date':
        datetimes = [nc2dt(j) for j in cube_out.coord(
            'time').units.num2date(cube_out.coord('time').points)]

    # Otherwise xaxis uses udunits
    else:
        datetimes = cube_out.coord('time').points

    # Choice of contour levels
    if input.plot_c_lims is not None:
        plot_levs = numpy.linspace(input.plot_c_lims[0], input.plot_c_lims[1], num=input.plot_nlev[0] + 1)
    elif (input.file_ts2 == None) & (not plot_plus_minus):
        plot_levs = plot_locator.tick_values(scipy.stats.scoreatpercentile(cube_out.data.flatten(), 10), 
                                             scipy.stats.scoreatpercentile(cube_out.data.flatten(), 90))
    else:
        plot_levs = errlev(cube_out, nlev=input.plot_nlev[0], log=input.plot_log, ex=input.plot_ext)

    # Contour plotting
    clev = mplt.pcolormesh(datetimes, cube_out.coord('depth').points, cube_out.data, 
                           norm=BoundaryNorm(plot_levs, ncolors=plot_colour.N, clip=True),
                           cmap=plot_colour)

    # If not set manually, set axis limits to constrained cube times and depths
    if input.plot_lims is not None:
        ylim(input.plot_lims)
    else:
        ylim([max(cube_out.coord('depth').points), min(cube_out.coord('depth').points)])

    xlim(datetimes.min(), datetimes.max())

    colorbar(clev)

elif cube_ts.ndim == 1:

    if input.line_colours is None:
#TODO: Find the actual colormap rather than doing this by hand
#NOTE: Method superceded by axes.set_prop_cycler at matplotlib 1.5
        # Brewer qualitative palette
        colours = [(226,28,28), (55,126,184), (77,175,74), (152,78,163), (255,127,0), (255,255,51), (166,86,40), (247,129,191), (153,153,153)]
        colour_cycler = cycle( [(i/256., j/256., k/256.) for (i, j, k) in colours])
    else:
        colour_cycler = cycle(input.line_colours)
    if input.line_styles is None:
        style_cycler = cycle(['-'])
    else:
        style_cycler = cycle(input.line_styles)
    if input.line_widths is None:
        width_cycler = cycle([1.5])
    else:
        width_cycler = cycle(input.width_styles)


    for i in cube_out:

        # x-axis uses dates
        if input.xaxis_type[0] == 'date':
            datetimes = [nc2dt(j) for j in i.coord(
                'time').units.num2date(i.coord('time').points)]

        # Otherwise x-axis uses udunits
        else:
            datetimes = i.coord('time').points

        # Set x-axis limits to constrained cube times
        try:
            datetimes_min = min(min(datetimes), datetimes_min)
            datetimes_max = max(max(datetimes), datetimes_max)
        except:
            datetimes_min = min(datetimes)
            datetimes_max = max(datetimes)

        xlim(datetimes_min, datetimes_max)

        # If a reference timeseries was specified, plot with scatter
        # (i == cube_out[0] doesn't always work for some reason..?)
        if (input.file_ts2 != None) and (id(i) == id(cube_out[0])):
            mplt.plot(datetimes, i.data, 'k.')
        else:
            mplt.plot(datetimes, i.data, color=next(colour_cycler), linewidth=next(width_cycler), linestyle=next(style_cycler))

    if input.plot_lims is not None:
        ylim(input.plot_lims)

    # Add input labels, else use file name arguments to label
    if input.plot_labels == None:
        plot_labels = [os.path.basename(i) for i in file_list[0:len(cube_out)]]
    else:
        plot_labels = input.plot_labels

    plot_leg = legend(plot_labels, fancybox=True, loc=0, fontsize=11)
    plot_leg.get_frame().set_alpha(0.5)
    
    if input.plot_title is None:
        if input.z_oper is not None:
            plot_title += ' (%s over depth)' % input.z_oper[0].lower()
        elif 'depth' in [i.name() for i in cube_out[0].coords()]:
            plot_title += ', depth = ' + \
                str(cube_out[0].coord('depth').points[0]) + 'm'

    fig.gca().grid(True)

if (input.coord_latlon is not None) and (input.plot_title is None):
    if len(input.coord_latlon) == 2:
        plot_title += ', region = [%sN, %sE]' % tuple(input.coord_latlon)
            
    elif len(input.coord_latlon) == 4:
        plot_title += ', region = [%s-%sN, %s-%sE]' % tuple(input.coord_latlon[0::2]+input.coord_latlon[1::2])
        
[i.set_fontsize(12) for i in mplt.gca().get_xticklabels()]
[i.set_fontsize(12) for i in mplt.gca().get_yticklabels()]

# Format of xaxis is date formatted
if input.xaxis_type[0] == 'date':
    no_days = (max(datetimes) - min(datetimes)).days

    # Year-month
    if no_days > 180:
        fig.gca().xaxis.set_major_locator(AutoDateLocator())
        fig.gca().xaxis.set_major_formatter(DateFormatter('%Y-%m'))

        # These are very arbitrary
        if no_days <= 3600:
            fig.gca().xaxis.set_minor_locator(MonthLocator())
        elif no_days <= 18000:
            fig.gca().xaxis.set_minor_locator(YearLocator())
    # Year-month-day
    else:
        fig.gca().xaxis.set_major_locator(DayLocator(interval=10))
        fig.gca().xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))
        fig.gca().xaxis.set_minor_locator(DayLocator())

    fig.autofmt_xdate()

# Otherwise the x axis label needs to include the time units
else:
    plot_xaxis += ' (%s)' % cube_out[0].coord('time').units
    fig.gca().xaxis.set_minor_locator(AutoMinorLocator())

# Add the y minor ticks
fig.gca().yaxis.set_minor_locator(AutoMinorLocator())

# Title, axis labels
title(plot_title)
mplt.xlabel(plot_xaxis, fontsize=14)
mplt.ylabel(plot_yaxis, fontsize=14)

if input.filename is None:
    mplt.show()
else:
    mplt.savefig(input.filename[0])
