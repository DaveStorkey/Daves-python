#! /usr/bin/env python

'''
Wrapper script for plot_nemo.py to plot the bathymetry in the
region of the sills/channels identified in "Sills of the Global
Ocean" by Simon Thompson (1995) and produce a web page.

Read in details of each location from a simple database.

@author: Dave Storkey
@date: Apr 2024
'''

import os
import datetime
import distutils
import xarray as xr
import numpy as np
import plot_nemo as pn
from cutout_domain import cutout_domain

def html_pieces():
    # Define html templates for the output web page

    html_head = """
<!DOCTYPE doctype PUBLIC "-//w3c//dtd html 4.0 transitional//en"><html><head>
   <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
   <title>Sills of the Global Ocean</title>
   <meta name="author" content="Dave Storkey"></head>
<style type="text/css">
table.gridtable {
	border-width: 2px;
	border-color: black;
	border-style: solid;
}
table.gridtable th {
	border-width: 2px;
	border-color: black;
	border-style: solid;
        text-style: bold;
        text-align: center;
        font-size: 1.9em;
}
table.gridtable td {
	border-width: 2px;
	border-color: black;
	border-style: solid;
        text-align: left;
}
</style>
   <body text="#000000" bgcolor="#ffffff" link="#cc0000" vlink="#663300" alink="#ff0000">
   <font size=2>
<hr>
<center>
<h1>
<font color="#000000">Sills of the Global Ocean</font></h1>
</center>
<hr>

    """

    html_tail = """
</table>
</center>

<hr width="100%">
<font color="#000000">Dave Storkey, <a href="mailto:dave.storkey@metoffice.gov.uk">dave.storkey@metoffice.gov.uk</a></font>
<br><font color="#000000">Created: DATE_TEXT</font>
<br>&nbsp;
</body></html>
    """

    return html_head, html_tail

def lonlat_to_float(string_in):

    try:
        # allow for the case where the input is already interpretable as a float
        float_out = float(string_in)
    except ValueError:
        strings=string_in.split()
        if len(strings) == 2:
            # degrees
            float_out=float(strings[0])
        elif len(strings) == 3:
            # degrees and arcminutes
            float_out=float(strings[0])+float(strings[1])/60.0
        elif len(strings) == 4:
            # degrees, arcminutes and arcseconds
            float_out=float(strings[0])+float(strings[1])/60.0+float(strings[2])/3600.0
        else:
            raise Exception("Error: could not interpret "+string_in+" as a latitude or longitude.")
        if strings[-1] == "N" or strings[-1] == "E":
            pass
        elif strings[-1] == "S" or strings[-1] == "W":
            float_out = -float_out
        else:
            raise Exception("Error: could not interpret "+string_in+" as a latitude or longitude.")

    return float_out 

def coordstring(lat,lon):

    if lat < 0.0:
        latstring = f"{-lat:.2f}S"
    else:
        latstring = f"{lat:.2f}N"

    if lon < 0.0:
        lonstring = f"{-lon:.2f}W"
    else:
        lonstring = f"{lon:.2f}E"

    cstring = latstring+" "+lonstring

    return cstring

def plot_sills(database=None, filenames=None, vars=None, titles=None, cutout=False,
               proj=None, plotsize=None, area_infile=None, area_var=None, clobber=False,
               xsection_dir=None):

    if database is None:
        raise Exception("Error: must specify database file.")

    if filenames is None:
        raise Exception("Error: must specify at least one input file.")
    elif not isinstance(filenames,list):
        filenames=[filenames]

    if vars is None:
        vars=len(filenames)*["Bathymetry"]
    elif not isinstance(vars,list):
        vars=[vars]
    if len(vars) == 1 and len(filenames) > 1:
        vars=len(filenames)*vars

    if titles is None:
        titles=""
    if not isinstance(titles,list):
        titles=[titles]
    if len(titles) == 1 and len(filenames) > 1:
        titles=len(filenames)*titles

    if cutout is None:
        cutout=False
    if not isinstance(cutout,list):
        cutout=[cutout]
    if len(cutout) == 1 and len(filenames) > 1:
        cutout=len(filenames)*cutout

    if plotsize is None:
        plotsize = 2.0

    if not isinstance(proj,list):
        proj=[proj]
    if len(proj) == 1 and len(filenames) > 1:
        proj=len(filenames)*proj

    # parse database
    s=[]
    sections=[]
    section_text=""
    section_dir="."
    try:
        with open(database) as db:
            count=0 
            for cline in db:
                if "====" in cline:
                    # new section
                    count=0
                    section_text=cline[2:].replace("====","").strip()
                    section_dir=section_text.replace(" ","")
                    sections.append(section_text)
                    if not os.path.isdir(section_dir):
                        os.makedirs(section_dir)
                if cline[0]=="#" or "|" not in cline:
                    # allow comments in database
                    continue
                count+=1
                att=cline.split('|')
                s.append({})
                s[-1]["index"] = count
                s[-1]["section text"] = section_text
                s[-1]["section dir"] = section_dir
                s[-1]["name"] = att[0].strip()
                s[-1]["depth"] = float(att[1].strip())
                s[-1]["width"] = att[2].strip()
                s[-1]["lat"] = lonlat_to_float(att[3].strip())
                s[-1]["lon"] = lonlat_to_float(att[4].strip())
                s[-1]["coordstring"] = coordstring(s[-1]["lat"],s[-1]["lon"])
                s[-1]["ref"] = att[5].strip()
    except FileNotFoundError:
        raise Exception("Error: file "+database+" not found.")

    print(len(sections)," sections found in database "+database)
    if area_infile is None:
        area_infile=filenames[0]
    if area_var is None:
        area_var=vars[0]

    # plot areas : projection, [W,E,S,N] extent of area to plot, longitudinal offset
    #              between marker cross and label, and vertbar (True or False) for each section.
    plot_areas = { "North Atlantic"         : ["merc",-90,10,-10,70,2.0,True],
                   "South Atlantic"         : ["merc",-70,30,-65,10,2.0,True],
                   "Pacific"                : ["pc180",110,300,-50,70,3.0,False],
                   "Indian Ocean"           : ["pc0"  ,30,120,-65,30,2.0,True],
                   "Arctic"                 : ["northps",None,None,None,None,3.0,True],
                   "Indonesian Throughflow" : ["pc0"  ,110,135,-15,10,0.5,True] }

    # loop over sections and plot maps with locations of sills marked
    for sec in sections:
        sectag=sec.replace(" ","")
        print("")
        print("Plotting map for "+sec)
        if sec not in plot_areas.keys():
            print("Can't find specs for "+sec+". Skipping.")
            continue
        else:
           spec=plot_areas[sec]
        secsills = [sill for sill in s if sill["section text"]==sec]
        print(len(secsills)," sills found in "+sec)
        print("")
        outfile=os.path.join(sectag,sectag+"_map.png")
        if len(secsills) > 0 and (clobber or not os.path.exists(outfile)):
            plot_text=[]
            plot_points=[]
            for sill in secsills:
                plot_text=plot_text+[sill["lon"]+spec[5],sill["lat"],str(sill["index"])]
                plot_points=plot_points+[sill["lon"]-0.1,sill["lat"]-0.1,sill["lon"]+0.1,sill["lat"]+0.1,
                                         sill["lon"]-0.1,sill["lat"]+0.1,sill["lon"]+0.1,sill["lat"]-0.1]
            (cslines, cscolor, csarrows) = pn.plot_nemo(filenames=area_infile,sca_names=area_var,proj=spec[0],
                                           plot_types="b",cmap="cividis_r",mnfld=0.0,mxfld=5500.0,nlevs=21,
                                           west=spec[1],east=spec[2],south=spec[3],north=spec[4],vertbar=spec[6],
                                           facecolor="lightgrey",text=plot_text,textbgcolor="white",textsize="xx-small",
                                           draw_points=plot_points,outfile=outfile)

    print(len(s)," sills found in database "+database)

    # loop over entries in database and bathymetry files and plot 
    for sill in s:
        print("")
        print("Plotting "+sill["name"])
        print("")
        draw_points = [sill["lon"]-0.1,sill["lat"]-0.1,sill["lon"]-0.1,sill["lat"]+0.1,
                       sill["lon"]-0.1,sill["lat"]+0.1,sill["lon"]+0.1,sill["lat"]+0.1,
                       sill["lon"]+0.1,sill["lat"]+0.1,sill["lon"]+0.1,sill["lat"]-0.1,
                       sill["lon"]+0.1,sill["lat"]-0.1,sill["lon"]-0.1,sill["lat"]-0.1]
        draw_fmt    = 4 * ["k-"]
        depmin = max(0.0   ,sill["depth"]-1000.0)
        depmax = min(5500.0,sill["depth"]+1000.0)
        compressed_name = sill["name"].replace(" ","").replace("/","")
        if xsection_dir is not None:
            xfilename = os.path.join(xsection_dir,"section_XTRAC_"+compressed_name+".dat")
            if os.path.exists(xfilename):
                with open(xfilename) as xfilein:
                    xfile = xfilein.read().splitlines()
                    npoints = int(xfile[1])
                    for point in range(npoints-1):                       
                        start_lon,start_lat = xfile[point+2].split()
                        end_lon  ,end_lat   = xfile[point+3].split()
                        draw_points = draw_points + [float(start_lon),float(start_lat),float(end_lon),float(end_lat)]
                        draw_fmt    = draw_fmt    + ["r-"]

        for filename,var,title,precut,prj in zip(filenames,vars,titles,cutout,proj):
            filestem=filename.replace(".nc","")
            outfile = os.path.join(sill["section dir"],title+"_"+compressed_name+".png")

            if clobber or not os.path.exists(outfile):
                south = sill["lat"]-plotsize
                north = sill["lat"]+plotsize
                west = sill["lon"]-plotsize
                east = sill["lon"]+plotsize
                if filename == filenames[0]:
                    print("south, north, west, east : ",south,north,west,east)

                print("")
                print("precut : ",precut)
                print("")
                if precut:
                    cutout_domain(infile=filename,invar=var,outfile=filestem+"_"+compressed_name+".nc",
                                 south=south,north=north,west=west,east=east)
                    filename=filestem+"_"+compressed_name+".nc"
                    south = None
                    north = None
                    west  = None
                    east  = None
              
                title=[title+": "+sill["name"],sill["coordstring"]+" "+str(sill["depth"])+"m"]
                (cslines, cscolor, csarrows) = pn.plot_nemo(filenames=filename,sca_names=var,
                                               plot_types="b",cmap="tab20b_r",proj=prj,
                                               mnfld=depmin,mxfld=depmax,clip=True,nlevs=21,
                                               west=west,east=east,south=south,north=north,
                                               vertbar=True,outfile=outfile,title=title,
                                               facecolor="white",draw_points=draw_points,
                                               draw_fmt=draw_fmt)

        if precut:
            try:
                os.remove(filename)
            except FileNotFoundError:
                pass

    # Create web page
    print("")
    print("Creating web page.")
    print("")
    html_head, html_tail = html_pieces()
    html_tail = html_tail.replace("DATE_TEXT","{:%d %b %Y %H:%M}".format(datetime.datetime.now()))
    ncols=str(len(titles)).strip()

    if os.path.exists("ocean_sills.html"):
        os.rename("ocean_sills.html","ocean_sills_OLD.html")
    with open("ocean_sills.html","x") as webfile:
        webfile.write(html_head)
        webfile.write("<hr><center><h2>Contents</h2></center><hr>\n")
        for sec in sections:
            sectag=sec.replace(" ","")
            webfile.write("<a href='#"+sectag+"'>"+sec+"</a><br>\n")
        for sec in sections:
            sectag=sec.replace(" ","")
            if sec != sections[0]:
                webfile.write("</table>")
            webfile.write("<hr><center><h2 id='"+sectag+"'>"+sec+"</h2></center><hr>\n")
            imagefile=os.path.join(sill["section dir"],sectag+"_map.png")
            webfile.write("<center><table COLS='2' WIDTH='90%' style='font-size:90%' NOSAVE >")
            webfile.write("<tr>\n")
            imagefile=os.path.join(sectag,sectag+"_map.png")
            webfile.write("<td><center><a href="+imagefile+"><img src="+imagefile+" alt='blah' height='800'></center></td>\n")
            webfile.write("<td>\n")
            for sill in [sill for sill in s if sill["section text"]==sec]:
                webfile.write("<a href='#"+sill["name"].replace(" ","")+"'>"+str(sill["index"])+". "+sill["name"]+"</a><br>\n")            
            webfile.write("</td></tr></table>\n")
            webfile.write("<center><table BORDER=1 COLS='"+ncols+"' WIDTH='90%' style='font-size:90%' NOSAVE >\n")
            for sill in [sill for sill in s if sill["section text"]==sec]:
                webfile.write("<tr>\n")
                webfile.write("<th colspan='"+ncols+"'><a id='"+sill["name"].replace(" ","")+"'>"+str(sill["index"])+". "+sill["name"]+"</a><br>"+sill["coordstring"]+" "+str(sill["depth"])+"m</th>\n")
                webfile.write("</tr>\n")
                webfile.write("<tr>\n")
                for title in titles:
                    imagefile=os.path.join(sill["section dir"],title+"_"+sill["name"].replace(" ","")+".png")
                    webfile.write("<td><center><a href="+imagefile+"><img src="+imagefile+" alt='blah' height='350'></a></center></td>\n")
                webfile.write("</tr>\n")
        webfile.write(html_tail)


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--database", action="store",dest="database",
        help="name of database file with details of locations to plot.")
    parser.add_argument("-i", "--filenames", action="store",dest="filenames",nargs="+",
        help="name(s) of database file with details of locations to plot.")
    parser.add_argument("-v", "--vars", action="store",dest="vars",nargs="+",
        help="name(s) of variables to plot")
    parser.add_argument("-t", "--titles", action="store",dest="titles",nargs="+",
        help="descriptions of the data from each input file eg. GEBCO2021 or ORCA025")
    parser.add_argument("-C", "--cutout", action="store",dest="cutout",nargs="+",
        type=lambda x:bool(distutils.util.strtobool(x)),
        help="if True cutout the data prior to plotting (for large data sets)")
    parser.add_argument("-P", "--proj", action="store",dest="proj",nargs="+",
        help="projection (default PlateCarree")
    parser.add_argument("-Z", "--plotsize", action="store",dest="plotsize",type=float,
        help="half size of plots in degrees (default = 2.0)")
    parser.add_argument("-S", "--xsection_dir", action="store",dest="xsection_dir",
        help="directory with definitions of Marine_Val sections to be plotted where applicable.")
    parser.add_argument("-A", "--area_infile", action="store",dest="area_infile",
        help="Which file to use for plotting section areas (default first input file).")
    parser.add_argument("--area_var", action="store",dest="area_var",
                        help="Which variable to use for plotting section areas (default first input variable).")
    parser.add_argument("-X", "--clobber", action="store_true",dest="clobber",
        help="if true, overwrite existing plots (default false)")
    args = parser.parse_args()
 
    plot_sills(database=args.database,filenames=args.filenames,vars=args.vars,
               titles=args.titles,cutout=args.cutout,proj=args.proj,plotsize=args.plotsize,
               area_infile=args.area_infile,area_var=args.area_var,clobber=args.clobber,
               xsection_dir=args.xsection_dir)
