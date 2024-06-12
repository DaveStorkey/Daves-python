#! /usr/bin/env python

'''
Wrapper script for plot_nemo.py to plot NEMO momentum
trends diagnostics and produce a web page. Mainly for 
testing the code. 

@author: Dave Storkey
@date: June 2024
'''

import os
import xarray as xr
import numpy as np
import plot_nemo as pn
import datetime

def html_pieces():

    html_head = '''
<!DOCTYPE doctype PUBLIC "-//w3c//dtd html 4.0 transitional//en"><html><head>
   <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
   <title>Momentum trend diagnostics in NEMO: eORCA025 GO8p6 configuration</title>
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
}
table.gridtable td {
	border-width: 2px;
	border-color: black;
	border-style: solid;
        text-align: left;
}
</style>
   <body text="#000000" bgcolor="#ffffff" link="#cc0000" vlink="#663300" alink="#ff0000">
   <font size=3>
<hr>
<center>
<h1>
<font color="#000000">Momentum trend diagnostics in NEMO:<br>TITLE</font></h1>
</center>
<hr>

<h2>Index</h2>

<a href="#intro">1. Introduction</a><br>
<a href="#u3d_100m">2. Balance of 3D "U" momentum trends at 100m</a><br>
<a href="#v3d_100m">3. Balance of 3D "V" momentum trends at 100m</a><br>
<a href="#u2d">4. Balance of depth-integrated "U" momentum trends</a><br>
<a href="#v2d">5. Balance of depth-integrated "V" momentum trends</a><br>

<br>
<!--------------------------------------------------------------------------------------------->
<hr>
<a name="intro"></a><h2>1. Introduction</h2>
<p>
Put your own text here.
</p>
<br>
    '''

    html_tail = '''
<hr width="100%">
<font color="#000000">Dave Storkey, <a href="mailto:dave.storkey@metoffice.gov.uk">dave.storkey@metoffice.gov.uk</a></font>
<br><font color="#000000">Last update: DATE_TIME</font>
<br>&nbsp;
</body></html>
    '''

    table_header = '''
<!--------------------------------------------------------------------------------------------->
<hr>
<a name="TABLE_NAME"></a><h2>SECTION_NUMBER. TABLE_TITLE</h2>
<p>
Note the different colour scales: the plots on each row of the table are plotted to the same scale and 
the first row of the table contains the largest trends. 
</p>
<center><table BORDER=1 COLS=NCOLS WIDTH="90%" NOSAVE >
    '''

    return html_head, html_tail, table_header

        
def process_fields(uvtype, trdtype, infile):
    '''
    Process input fields if necessary
    '''
    indata = xr.open_dataset(infile)    

    if trdtype == "adv":
        if uvtype+"_adv" in indata.data_vars:
            # nothing to do
            return
        # ADV = RVO + KEG + ZAD (for vector invariant form)
        for intype in "rvo", "keg", "zad":
            if uvtype+"_"+intype not in indata.data_vars:
                raise Exception("Could not find "+uvtype+"_"+intype+" in file "+infile)
        trd_adv = indata.data_vars[uvtype+"_rvo"] + indata.data_vars[uvtype+"_keg"] + indata.data_vars[uvtype+"_zad"]         
        trd_adv.to_netcdf(infile,mode="a")

    elif trdtype == "res":
        # TBD!
        return
    
def plot_mmtm_trends(infiles=None, tag=None, title=None, level=None, clobber=False):

    # Tuples to specify what is plotted, contour limits and the layout of tables
    # on the web page. Change according to your purposes. 
    #
    # These default layout and contour limits are appropriate for the first month 
    # of spin up from rest of ORCA025. 

#    trd3d_table = ( ( 1.0e-05 , "hpg", "pvo", "zdf" ),
#                    ( 1.0e-07 , "adv", "ldf", "tot" ),
#                    ( 1.0e-11 , "res" ), 
#                  )

    trd3d_table = ( ( 1.0e-05 , "hpg", "pvo", "zdf" ),
                    ( 1.0e-07 , "ldf", "tot", "adv" ),
                  )

    trd2d_table = ( ( 5.0e-02 , "spg2d_huv", "pvo2d_huv", "frc2d_huv" ),
                    ( 5.0e-04 , ""        , "tau2d_huv", ""         ), 
                    ( 1.0e-04 , "tfr2d_huv", "tot2d_huv", ""         ),
                    ( 5.0e-06 , "bfr2d_huv", ""        , ""         ), 
                  )

    trdtitle = { "hpg" : "horizontal pressure gradient"          , 
                 "pvo" : "planetary vorticity"                   ,
                 "zdf" : "vertical diffusion"                    ,
                 "adv" : "advection"                             ,
                 "ldf" : "lateral diffusion"                     ,
                 "tot" : "total trend excluding time filter"     ,
                 "res" : "total minus sum of trends"             ,
                 "spg2d_huv" : "surface pressure gradient"        ,
                 "pvo2d_huv" : "planetary vorticity"              ,
                 "frc2d_huv" : "forcing from baroclinic solution" ,
                 "tau2d_huv" : "trend due to wind stress"         ,
                 "tfr2d_huv" : "trend due to ice-ocean stress"    , 
                 "tot2d_huv" : "total depth-integrated trend"     ,
                 "bfr2d_huv" : "trend due to bottom friction"     ,
                 "res2d_huv" : "total minus sum of trends (2D)"   ,
               }    

    if level is None:
        # default to something close to 100m (for the Drakkar 75 level set).
        level = 24

    # make plots
    for trdtable in (trd3d_table,trd2d_table):
        for uvtype, infile in zip(("utrd","vtrd"),infiles):
            for row in trdtable:
                mnfld=-1.0*row[0]
                mxfld=row[0]
                for trdtype in row[1:]:
                    trdfld = trdtype.replace("huv","h"+uvtype[0])

                    # generate derived trends (if not already in file)
                    for trd in ("adv", "res"):
                        if trd in trdfld:
                            process_fields(uvtype, trdfld, infile)

                    if "2d" in trdtype:
                        lev_to_plot = 0
                    else:
                        lev_to_plot = level

                    # plot trend
                    outfile = uvtype+"_"+trdfld+".png"
                    if len(trdtype) > 0 and ( clobber or not os.path.exists(outfile) ):
                        print(">>> plotting "+uvtype+" "+trdtype)
                        (cslines, cscolor, csarrows) = \
                                pn.plot_nemo(filenames=infile, sca_names=uvtype+"_"+trdfld, 
                                plot_types="c", cmap="RdBu_r", nlevs=21, level=lev_to_plot, 
                                mnfld=mnfld, mxfld=mxfld, scientific=True, 
                                title=uvtype.upper()+" "+trdtype.upper(), outfile=outfile)

    # generate web page
    html_head, html_tail, table_header = html_pieces()       
    if title is not None:
        title="<br>"+title
    else:
        title=""
    html_head = html_head.replace("TITLE",title)
    html_tail = html_tail.replace("DATE_TIME","{:%d %b %Y %H:%M}".format(datetime.datetime.now()))
    if tag is not None or tag == "":
        tag=tag+"_"
    else:
        tag=""    

    sec_count=1
    with open(tag+"mmtm_trends.html","w") as webfile:
        webfile.write(html_head)
        for ndim, trdtable in zip(("3D","2D"),(trd3d_table,trd2d_table)):
            ncols_table = max([len(row) for row in trdtable])
            for uv, uvtype in zip(("U","V"),("utrd","vtrd")):
                sec_count += 1
                webfile.write(table_header.replace("TABLE_NAME",uvtype+ndim)\
                                          .replace("SECTION_NUMBER",str(sec_count))\
                                          .replace("TABLE_TITLE",ndim+" "+uv+" momentum trends at 100m")\
                                          .replace("NCOLS",str(ncols_table)) )
                for row in trdtable:
                    mxfld = str(row[0]) 
                    webfile.write("<tr>\n") 
                    for trd in row[1:]:
                         trdfld = trd.replace("huv","h"+uv.lower())
                         if trdfld == "":
                             webfile.write("<th></th>\n")
                         else:
                             webfile.write("<th>\n "+trdtitle[trd]+" ("+trdfld.upper()+")<br>\n")
                             webfile.write("-"+mxfld+" m2/s2 -- +"+mxfld+" m2/s2\n</th>")
                    webfile.write("</tr>\n") 
                    webfile.write("<tr>\n") 
                    for trd in row[1:]:
                         trdfld = trd.replace("huv","h"+uv.lower())
                         if trd == "":
                             webfile.write("<td></td>\n")
                         else:
                             webfile.write("<td>\n<center><a href='"+uvtype+"_"+trdfld+".png'>")
                             webfile.write("<img src='"+uvtype+"_"+trdfld+".png' alt='"+uvtype+"_"+trdfld+"' height='350'></a></center>\n</td>")
                    webfile.write("</tr>\n") 
                webfile.write("</table></center>\n</p>\n")
        webfile.write(html_tail)

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles",action="store",dest="infiles",nargs=2,
        help="names of two input files (U and V).")
    parser.add_argument("-t", "--tag",action="store",dest="tag",
        help="tag for this experiment")
    parser.add_argument("-T", "--title",action="store",dest="title",
        help="title for the web page")
    parser.add_argument("-L", "--level",action="store",type=int, dest="level",
        help="model level number to plot (for 3D fields)")
    parser.add_argument("-X", "--clobber", action="store_true",dest="clobber",
        help="if true, overwrite existing plots (default false)")
    args = parser.parse_args()
 
    plot_mmtm_trends(infiles=args.infiles,tag=args.tag,title=args.title,level=args.level,clobber=args.clobber)

