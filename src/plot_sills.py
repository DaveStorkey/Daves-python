#! /usr/bin/env python

'''
Wrapper script for plot_nemo.py to plot the bathymetry in the
region of the sills/channels identified in "Sills of the Global
Ocean" by Simon Thompson (1995).

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

def plot_sills(database=None, filenames=None, vars=None, titles=None, cutout=False,
               proj=None, clobber=False):

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

    if not isinstance(proj,list):
        proj=[proj]
    if len(proj) == 1 and len(filenames) > 1:
        proj=len(filenames)*proj

    # parse database
    s=[]
    section_text=""
    section_dir="."
    try:
        with open(database) as db:
            for cline in db:
                if "====" in cline:
                    # new section
                    section_text=cline[2:].replace("====","").strip()
                    section_dir=section_text.replace(" ","")
                    if not os.path.isdir(section_dir):
                        os.makedirs(section_dir)
                if cline[0]=="#" or "|" not in cline:
                    # allow comments in database
                    continue
                att=cline.split('|')
                s.append({})
                s[-1]["section text"] = section_text
                s[-1]["section dir"] = section_dir
                s[-1]["name"] = att[0].strip()
                s[-1]["depth"] = float(att[1].strip())
                s[-1]["width"] = att[2].strip()
                s[-1]["lat"] = lonlat_to_float(att[3].strip())
                s[-1]["lon"] = lonlat_to_float(att[4].strip())
                s[-1]["ref"] = att[5].strip()
    except FileNotFoundError:
        raise Exception("Error: file "+database+" not found.")

    print(len(s)," entries found in database "+database)

    # loop over entries in database and bathymetry files and plot 
    for sill in s:
        print("")
        print("Plotting "+sill["name"])
        print("")
        draw_points = [sill["lon"]-0.1,sill["lat"]-0.1,sill["lon"]-0.1,sill["lat"]+0.1,
                       sill["lon"]-0.1,sill["lat"]+0.1,sill["lon"]+0.1,sill["lat"]+0.1,
                       sill["lon"]+0.1,sill["lat"]+0.1,sill["lon"]+0.1,sill["lat"]-0.1,
                       sill["lon"]+0.1,sill["lat"]-0.1,sill["lon"]-0.1,sill["lat"]-0.1]
        depmin = max(0.0   ,sill["depth"]-1000.0)
        depmax = min(5500.0,sill["depth"]+1000.0)
        compressed_name = sill["name"].replace(" ","")

        for filename,var,title,precut,prj in zip(filenames,vars,titles,cutout,proj):
            filestem=filename.replace(".nc","")
            outfile = os.path.join(sill["section dir"],title+"_"+compressed_name+".png")

            if clobber or not os.path.exists(outfile):
                south = sill["lat"]-2.5
                north = sill["lat"]+2.5
                west = sill["lon"]-2.5
                east = sill["lon"]+2.5
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

                (cslines, cscolor, csarrows) = pn.plot_nemo(filenames=filename,sca_names=var,
                                               plot_types="b",cmap="tab20b_r",proj=prj,
                                               mnfld=depmin,mxfld=depmax,west=west,east=east, 
                                               south=south,north=north,vertbar=True,nlevs=21,
                                               outfile=outfile,title=title+": "+sill["name"],
                                               facecolor="white",draw_points=draw_points)

#        if precut:
#            try:
#                os.remove(filename)
#            except FileNotFoundError:
#                pass

    # Create web page
    print("")
    print("Creating web page.")
    print("")
    html_head, html_tail = html_pieces()
    html_tail = html_tail.replace("DATE_TEXT","{:%d %b %Y %H:%M}".format(datetime.datetime.now()))
    ncols=str(len(s)).strip()

    if os.path.exists("ocean_sills.html"):
        os.rename("ocean_sills.html","ocean_sills_OLD.html")
    with open("ocean_sills.html","x") as webfile:
        webfile.write(html_head)
        section_text=""
        for sill in s:
            if sill["section text"] != section_text:
                if len(section_text) > 0:
                    webfile.write("</table>")
                section_text=sill["section text"]
                webfile.write("<hr><center><h2>"+section_text+"</h2></center><hr>")
                webfile.write("<center><table BORDER=1 COLS='"+ncols+"' WIDTH='90%' style='font-size:90%' NOSAVE >")
            webfile.write("<tr>\n")
            webfile.write("<th>"+sill["name"]+"</th>\n")
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
    parser.add_argument("-X", "--clobber", action="store_true",dest="clobber",
        help="if true, overwrite existing plots (default false)")
    args = parser.parse_args()
 
    plot_sills(database=args.database,filenames=args.filenames,vars=args.vars,
               titles=args.titles,cutout=args.cutout,proj=args.proj,clobber=args.clobber)
