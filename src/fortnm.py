#!/usr/bin/env python
"""
    Manipulates Fortran Namelists using the ROSE namelist parser
    Requires python 2.7 for the OrderedDict class. 
"""

__author__ = 'Dave Storkey (dave.storkey@metoffice.gov.uk)'

import sys
import namelist

from collections import OrderedDict


def cmp(x, y):
    if x < y:
        return -1
    elif x > y:
        return 1
    else:
        return 0

class NamelistFileDict(OrderedDict):
    """
    Uses the ROSE namelist parser to parse a file full of namelists and return the result
    as a nested Ordered Dictionary.
    """

    def __init__(self,filename=None):
        OrderedDict.__init__(self)
        if filename is not None:
            self.update(self.parse(filename))
   
    def parse(self,filename):
        """Parse namelist file into NamelistFileDict object using ROSE namelist parser"""
        mynmlfile  = OrderedDict()
        groups = namelist.parse([filename])
        for group in groups:
            mynmlfile[group.name] = OrderedDict()
            for object in group.objects:
                mynmlfile[group.name][object.lhs] = object.get_rhs_as_string()
        return mynmlfile

    def display(self):

        for nml in self.keys():
            print('!==========================')
            print('&'+nml)
            print('!==========================')
            for par in self[nml].keys():
                print(par," = ",self[nml][par])
            print('/')
            print()

def difference_namelists(nmlfile1, nmlfile2):
    """Displays differences between two namelist files"""

    # First find the namelists in common
    nml_commonlist = [item for item in nmlfile1.keys() if item in nmlfile2.keys()]
    nml_exlist1 = [item for item in nmlfile1.keys() if item not in nmlfile2.keys()]
    nml_exlist2 = [item for item in nmlfile2.keys() if item not in nmlfile1.keys()]
    if len(nml_exlist1) > 0:
        print()
        print('==============================================')
        print('The following namelists only exist in file 1 :')
        for item in nml_exlist1:    
            print(item,',',)
        print()
        print('==============================================')
    if len(nml_exlist2) > 0:
        print()
        print('==============================================')
        print('The following namelists only exist in file 2 :')
        for item in nml_exlist2:    
            print(item,',',)
        print()
        print('==============================================')
    print()
    for nml in nml_commonlist:
        print()
        print('====================')
        print('Namelist',nml,':')
        print('====================')
        par1 = nmlfile1[nml]
        par2 = nmlfile2[nml]
        par_commonlist = [item for item in par1.keys() if item in par2.keys()]
        par_exlist1 = [item for item in par1.keys() if item not in par2.keys()]
        par_exlist2 = [item for item in par2.keys() if item not in par1.keys()]
        for item in par_commonlist:
            value1 = nmlfile1[nml][item]
            value2 = nmlfile2[nml][item]
            if cmp(value1,value2) != 0:
                print(item,'=',value1,'<=>',value2)
        for item in par_exlist1:
            value1 = nmlfile1[nml][item]
            print(item,'=',value1,'<=>','UNSPECIFIED')
        for item in par_exlist2:
            value2 = nmlfile2[nml][item]
            print(item,'=','UNSPECIFIED','<=>',value2)
    return None

def merge_namelists(nmlfile1, nmlfile2):
    """Merges two namelist files"""

    # 1. NB. This is just another pointer. We are modifying nmlfile1 in place.
    #    Did try doing nmlfile1.copy() but it didn't like it!
    nml_merged = nmlfile1

    # 2. Add in any namelists that are in the second file but not the first file 
    nml_exlist2 = [item for item in nmlfile2.keys() if item not in nmlfile1.keys()]
    for nml in nml_exlist2:
        nml_merged.update(item=nmlfile2[nml])

    # 3. Merge the namelists that the two files have in common using the python dict.update
    #    method. The second namelist will take precedence where the two files have parameters in common. 
    nml_commonlist = [item for item in nmlfile1.keys() if item in nmlfile2.keys()]
    for nml in nml_commonlist:
        nml_merged[nml].update(nmlfile2[nml])

    # 4. Write out the merged namelist
    nml_merged.display()    


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("filename1", help="name of first input file")
    parser.add_argument("filename2", nargs="?", default=None, help="name of second input file (optional)")
    parser.add_argument("-d", "--difference", action="store_const",dest="action",const="difference",default="difference",
                    help="difference namelists")
    parser.add_argument("-m", "--merge", action="store_const",dest="action",const="merge",default="difference",
                    help="merge namelists")
    args = parser.parse_args()

    nmlfile1 = NamelistFileDict(args.filename1)
    if args.filename2 is not None:
        nmlfile2 = NamelistFileDict(args.filename2)
        if args.action == "difference":
            difference_namelists(nmlfile1,nmlfile2)
        if args.action == "merge":
            merge_namelists(nmlfile1,nmlfile2)
    else:
        nmlfile1.display()
