#! /usr/bin/env python

'''
Routine to plot timeseries of transports through straits
previously calculated and pickled by calc_transport_timeseries.py

@author: Dave Storkey
@date: April 2017
'''

import pickle
from collections import OrderedDict
import matplotlib
import matplotlib.pyplot as plt

def plot_transport_timeseries(pickle_file, outfile=None, section=None):

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

    file_id = open(pickle_file,'rb')
    # runs is a list of dictionaries, one for each run.
    runs = pickle.load(file_id)

    if section is None:
        raise Exception('Error: must specify section to plot')

    labels_for_legend=[]
    for run in runs:
        secpkl = all_sections[section]
        if secpkl[0] not in run['data'].keys():
            raise Exception('Error: requested section not in data for '+run['jobid'])
        if secpkl[1] == 'pos':
            plt.plot((run['data'][secpkl[0]][0][:,0])*1e-06)
        elif secpkl[1] == 'neg':
            plt.plot((run['data'][secpkl[0]][0][:,1])*1e-06)
        elif secpkl[1] == 'net':
            plt.plot((run['data'][secpkl[0]][0][:,0]+run['data'][secpkl[0]][0][:,1])*1e-06)
        labels_for_legend.append(run['jobid'])
        
    plt.legend(labels_for_legend,loc=3,fontsize='medium')
    plt.gcf().suptitle(secpkl[2], fontsize=12, y=0.95)
    plt.gca().set_ylabel('transport (Sv)')

    if outfile is not None:
        matplotlib.rcParams['font.size'] =8
        plt.savefig(outfile,dpi=200)
    else:
        plt.show()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("pickle_file", help="name of input picklefile")
    parser.add_argument("-o", "--outfile", action="store",dest="outfile",default=None,
                    help="plot to output file - filetype by extension. If unset plot to GUI.")
    parser.add_argument("-s", "--section", action="store",dest="section",default=None,
                    help="section to plot. If unset plot first section in file.")

    args = parser.parse_args()

    plot_transport_timeseries( args.pickle_file, outfile=args.outfile, section=args.section )
