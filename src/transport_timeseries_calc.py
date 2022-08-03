#! /usr/bin/env python
'''
Routine to calculate timeseries of transports through sections.
Based on code by Tim Graham. 

Created April 2017

@author: Dave Storkey 
'''
import pickle
import transports
import numpy as np
import glob
from def_sections_hadgem3 import def_sections_hadgem3
        
def transport_timeseries_calc(outfile=None, sections_list=None):

    # List of dictionaries to define data sets (user defined)
    runs=[
           {'jobid':'u-ai599',
           'description' : 'GC3.1',
           'resolution'  : '025e',
           'meshfile':'/data/cr1/hadtd/CMA/HadGEM3/ORCA025/mesh_mask_GO6.nc',
           'datadir':'/scratch/frsy/transports/u-ai599_annual_means/',
           'Tstring':'/nemo_ai599*T.nc',
           'Ustring':'/nemo_ai599*U.nc',
           'Vstring':'/nemo_ai599*V.nc',
           'data'   : {} },
           {'jobid':'u-au939',
           'description' : 'GC4',
           'resolution'  : '025e',
           'meshfile':'/data/cr1/hadtd/CMA/HadGEM3/ORCA025/mesh_mask_GO6.nc',
           'datadir':'/scratch/frsy/transports/u-au939_annual_means/',
           'Tstring':'/nemo_au939*T.nc',
           'Ustring':'/nemo_au939*U.nc',
           'Vstring':'/nemo_au939*V.nc',
           'data'   : {} },
#          {'jobid':'u-ak817',
#           'description' : 'GO6-GSI8.1',
#           'resolution'  : '025e',
#           'meshfile':'/data/cr1/hadtd/CMA/HadGEM3/ORCA025/mesh_mask_GO6.nc',
#           'datadir':'/scratch/frsy/ocean_transports/u-ak817_5day_means/',
#           'Tstring':'/*5d*grid_T*.nc',
#           'Ustring':'/*5d*grid_U*.nc',
#           'Vstring':'/*5d*grid_V*.nc',
#           'data'   : {} }
           ]
   
    # Dictionary to define sections (leave this alone unless you're defining a new section)
    # Each section defined by a list : [ name of section in definitions file, density range if applicable ]
    all_sections={'dms_net'      : ['dms', None ],
                  'dms_overflow' : ['dms',(27.8,60)],
                  'if_net'       : ['if', None ],
                  'if_overflow'  : ['if',(27.8,60)],
                  'fs_net'       : ['fs', None ],
                  'fs_overflow'  : ['fs',(27.8,60)],
                  'bering'       : ['bering', None ],
                  'davis'        : ['davis_strait', None ],
                  'itf'          : ['itf', None ],
                  'acc'          : ['acc', None ],
                  'lab_ar7w'     : ['lab_ar7w', None ],
                  's_greenland'  : ['s_greenland', None ],
                  'gibraltar'    : ['gibraltar', None ]}

    if outfile is None:
        outfile = "transport_timeseries"
        for run in runs:
            outfile=outfile+'_'+run['jobid']
        outfile=outfile+'.pkl'

    if sections_list is None:
        sections_list = all_sections.keys()

    #Loop through the model runs and add transports to the dictionaries 
    for run in runs:
        Tfiles=glob.glob(run['datadir']+run['Tstring'])
        Ufiles=glob.glob(run['datadir']+run['Ustring'])
        Vfiles=glob.glob(run['datadir']+run['Vstring'])
        Tfiles.sort()
        Ufiles.sort()
        Vfiles.sort()
        
        #Define the sections
        for section in sections_list:
            if section not in all_sections.keys():
               raise Exception('Error: requested cross section '+section+' not defined.')

            (sec,sw_ne) = def_sections_hadgem3(all_sections[section][0],run['resolution'])
    
            coords=[coord for coord in zip(sec[0],sec[1])]
            if all_sections[section][1] is not None:
                dens_range = np.array(all_sections[section][1])
                run['data'][section]=transports.calc_trans_rho(Tfiles,Tfiles,Ufiles,Vfiles,
                                          coords,sec[2],sec[3],sec[4],sec[5],
                                          dens_range,run['meshfile'],sw_ne) 
            else:
                run['data'][section]=transports.calc_trans_depth(Tfiles,Tfiles,Ufiles,Vfiles,
                                          coords,sec[2],sec[3],sec[4],sec[5],
                                          [0,-1],run['meshfile'],sw_ne,tracers=False) 
        
    fid = open(outfile,'wb')
    pickle.dump(runs,fid)
    fid.close()
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",default=None,
                    help="name of output pickle file")
    parser.add_argument("-s", "--sections", action="store",dest="sections_list",nargs="+",default=None,
                    help="list of sections to calculate. (If unset calculate all).")

    args = parser.parse_args()

    transport_timeseries_calc(outfile=args.outfile, sections_list=args.sections_list)        
