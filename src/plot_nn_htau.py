#! /usr/bin/env python

'''
Plot nn_htau latitudinal dependence.

Created Jan 2022

@author: Dave Storkey
'''

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def plot_nn_htau(profiles=None,outfilename=None):

    if profiles is None:
        profiles=[1,4]

    # latitude in radians:
    phi = (np.arange(201)-100) * 0.005 * np.pi
    htau = np.zeros(phi.shape)

    for profile in profiles:

        print('profile : ',profile)

        if profile == 1:
            htau = np.maximum( 0.5, np.minimum( 10., 45. * np.abs(np.sin(phi[:])) ) )
        elif profile == 4:
            htau = np.where(phi < 0.0, np.maximum( 0.5, np.minimum( 30., 45. * np.abs(np.sin(phi[:])) ) ),
                                       np.maximum( 0.5, np.minimum( 10., 45. * np.abs(np.sin(phi[:])) ) ) )  
        elif profile == 5:
            htau = np.maximum( 0.5, np.minimum( 10., 45. * np.abs(np.sin(phi[:])) ) )
            htau = htau + np.where(phi > -0.698, 0.0,
                                   np.maximum( 0.5, np.minimum( 20., 135. * np.abs(np.sin(phi[:]+0.698)) ) ) )

        plt.plot(phi*180./np.pi,htau)

    plt.gca().set_xlabel("latitude (deg)")
    plt.gca().set_ylabel("htau (m)")

    #plt.legend(legend,loc=1,fontsize='medium')

    if outfilename is not None:
        matplotlib.rcParams['font.size'] =8
        plt.savefig(outfilename,dpi=200)
    else:
        plt.show()
        
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--profiles", action="store", dest="profiles", nargs="+", type=int,
                    help="profiles to plot")
    parser.add_argument("-o", "--outfile", action="store", dest="outfilename", 
                    help="name of output file")

    args = parser.parse_args()

    plot_nn_htau(profiles=args.profiles,outfilename=args.outfilename)
