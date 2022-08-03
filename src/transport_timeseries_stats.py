#! /usr/bin/env python

'''
Routine to calculate statistics of timeseries of transports through straits
previously calculated and pickled by calc_transport_timeseries.py

@author: Dave Storkey
@date: April 2017
'''

import pickle
import numpy as np
from collections import OrderedDict

def plot_transport_timeseries(pickle_files, outfile=None, table=False):

    # Ordered dictionary to define sections (leave this alone unless you're defining a new section)
    # Each section defined by a list : [ name of section in pickle file, positive negative or net flow, description ]
    all_sections=OrderedDict([('dms_net'      , ['dms_net'     , 'net', 'Denmark Strait net'        ]),
                              ('dms_overflow' , ['dms_overflow', 'net', 'Denmark Strait overflow'   ]),
                              ('if_net'       , ['if_net'      , 'net', 'Iceland-Faeroes net'       ]),
                              ('if_overflow'  , ['if_overflow' , 'net', 'Iceland-Faeroes overflow'  ]),
                              ('fs_net'       , ['fs_net'      , 'net', 'Faeroes-Scotland net'      ]),
                              ('fs_overflow'  , ['fs_overflow' , 'net', 'Faeroes-Scotland overflow' ]),
                              ('bering_net'   , ['bering'      , 'net', 'Bering Strait net'         ]),
                              ('davis_net'    , ['davis'       , 'net', 'Davis Strait net'          ]),
                              ('davis_northwd', ['davis'       , 'pos', 'Davis Strait northward'    ]),
                              ('davis_southwd', ['davis'       , 'neg', 'Davis Strait southward'    ]),
                              ('itf'          , ['itf'         , 'net', 'Indonesian Throughflow'    ]),
                              ('acc_net'      , ['acc'         , 'net', 'Drake Passage net'         ]),
                              ('acc_eastwd'   , ['acc'         , 'pos', 'Drake Passage eastward'    ]),
                              ('acc_westwd'   , ['acc'         , 'neg', 'Drake Passage westward'    ]),
                              ('lab_ar7w'     , ['lab_ar7w'    , 'net', 'Labrador Sea AR7 section'  ]),
                              ('s_greenland'  , ['s_greenland' , 'net', 'South Greenland'           ]),
                              ('gib_inward'   , ['gibraltar'   , 'pos', 'Gibraltar Strait inward'   ]),
                              ('gib_outward'  , ['gibraltar'   , 'neg', 'Gibraltar Strait outward'  ]) ] )

    runs=[]
    for infile in pickle_files:
        with open(infile,'rb') as file_id:
            # runs is a list of dictionaries, one for each run.
            data_in = pickle.load(file_id)
        runs = runs + data_in

    if table:
        print '<html>'
        print '<table width="90%" border="1" cellpadding="5" cellspacing="0">'
        print '<tr>'
        print '<td>'
        runs_all_keys=[]
        for run in runs:
            print '<td><b> '+run['jobid']+' </b></td>'
            runs_all_keys = runs_all_keys+run['data'].keys()
        print '</tr>'
        for sec_to_plot in all_sections.keys():
            sec_in_data = all_sections[sec_to_plot][0]
            sec_type = all_sections[sec_to_plot][1]
            if sec_in_data not in runs_all_keys:
                pass
            else:
                print '<tr>'
                print '<td><b> '+all_sections[sec_to_plot][2]+' </b></td>'
                for run in runs:
                    data=run['data']
                    if sec_in_data not in data.keys():
                        print '<td></td>'
                    elif sec_type == 'pos':
                        print '<td> '+str(np.around(np.mean((data[sec_in_data][0][:,0])*1e-06),decimals=2))+' &plusmn; '+ \
                                      str(np.around(np.std((data[sec_in_data][0][:,0])*1e-06),decimals=2)),' Sv</td>'
                    elif sec_type == 'neg':
                        print '<td> '+str(np.around(np.mean((data[sec_in_data][0][:,1])*1e-06),decimals=2))+' &plusmn; '+ \
                                      str(np.around(np.std((data[sec_in_data][0][:,1])*1e-06),decimals=2)),' Sv</td>'
                    elif sec_type == 'net':
                        print '<td> '+str(np.around(np.mean((data[sec_in_data][0][:,0]+data[sec_in_data][0][:,1])*1e-06),decimals=2))+' &plusmn; '+ \
                                      str(np.around(np.std((data[sec_in_data][0][:,0]+data[sec_in_data][0][:,1])*1e-06),decimals=2)),' Sv</td>'
                print '</tr>'
        print '</table>'
        print '</html>'
    else:
        for run in runs:
            print 'Stats for run : ',run['jobid']
            data=run['data']
            for sec_to_plot in all_sections.keys():
                sec_in_data = all_sections[sec_to_plot][0]
                sec_type = all_sections[sec_to_plot][1]
                if sec_in_data not in data.keys():
                    pass
                elif sec_type == 'pos':
                    print '  ',sec_to_plot,':'
                    print '     mean                : ',np.mean((data[sec_in_data][0][:,0])*1e-06),' Sv'
                    print '     standard deviation  : ',np.std((data[sec_in_data][0][:,0])*1e-06),' Sv'
                    print ' '
                elif sec_type == 'neg':
                    print '  ',sec_to_plot,':'
                    print '     mean                : ',np.mean((data[sec_in_data][0][:,1])*1e-06),' Sv'
                    print '     standard deviation  : ',np.std((data[sec_in_data][0][:,1])*1e-06),' Sv'
                    print ' '
                elif sec_type == 'net':
                    print '  ',sec_to_plot,':'
                    print '     mean                : ',np.mean((data[sec_in_data][0][:,0]+data[sec_in_data][0][:,1])*1e-06),' Sv'    
                    print '     standard deviation  : ',np.std((data[sec_in_data][0][:,0]+data[sec_in_data][0][:,1])*1e-06),' Sv'    
                    print ' '

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("pickle_files", nargs="+", help="name(s) of input pickle file(s)")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",default=None,
                    help="plot to output file - filetype by extension. If unset plot to GUI.")
    parser.add_argument("-t", "--table", action="store_true",dest="table",default=False,
                    help="output results as an HTML table.")

    args = parser.parse_args()

    plot_transport_timeseries( args.pickle_files, outfile=args.outfile, table=args.table)
