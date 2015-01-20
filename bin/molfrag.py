#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
# molfrag.py - separate each different molecule from an xyz file
# Copyright (C) yyyy  name of author

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

# Author: Riccardo Petraglia <riccardo.petraglia@gmail.com>

"""Separate each molecule in a xyz file.

Each molecule present in a single xyz file will be separated from 
all the others with only the information included in the file.

:Author: `Riccardo Petraglia <riccardo.petraglia@gmail.com>

Requirements
------------
* numpy
* elements (distribuited with this software)

Examples
--------
Go in the $root/example directory and run
$python molfrag.py lastframe.xyz

"""

import os
import sys

def main():

    root = os.path.dirname(os.path.realpath(__file__))
    root = os.path.split(root)[0]
    sys.path.append(root)

    import lib.kernel  as krn
    
    args = __parser()

    krn.fragKern.factor = args.coeff
    
    atomnos, atomcoords, cstring = krn.io.readxyz(args.xyzfile)

    voisins = krn.fragKern.neighbours(atomnos, atomcoords)
    monomers = krn.fragKern.molsplit(voisins)
    
    i = 0
    for m in monomers.keys():
        monomer = []
        for atom in monomers[m]:
            monomer.append([atomnos[atom]]+atomcoords[atom])
#        fname = '%s_mono%i.xyz'%(args.out_prefix, m+1)
        #        writexyz(fname, monomer, cstring)
        print 'Mon %i molecular mass: %12.6f' % tuple([i,krn.structparams.molecular_mass(monomers[i],atomnos)])
        i += 1

    # print "Mon 1 molecular mass: ",
    # print "Mon 2 molecular mass: ",frg.molecular_mass(monomers[1],atomnos)

    # print "Mon 1: ", sorted(monomers[0])
    # print "Mon 2: ", sorted(monomers[1])

    for a in monomers.keys():
        outfn = args.out_prefix+str(a)+'.xyz'
        xyz_w = '%5d\n\n' % len(monomers[a])
        for b in sorted(monomers[a]):
            xyz_w += '%3s  %12.6f  %12.6f  %12.6f\n' % tuple([atomnos[b]]+atomcoords[b])
            
        open(outfn,'w').write(xyz_w)



    
#def __writexyz()
    


def __parser():
    import argparse
    parser = argparse.ArgumentParser(version='%prog 1.0', 
                                     description='Separe different molecules from the same xyz file.')
    
    parser.add_argument('-o','--out-prefix',
                        action = 'store',
                        dest = 'out_prefix',
                        type = str,
                        metavar = '<STRING>',
                        default = 'frag_',
                        help = 'prefix of the output files (default: frag_)')

    parser.add_argument('-c','--coefficient',
                        action = 'store',
                        dest = 'coeff',
                        type = float,
                        metavar = '<FLOAT>',
                        default = '1.05',
                        help = 'select a multiplicating factor to for the covalent radii')

    parser.add_argument('xyzfile',
                        action = 'store',
                        type = str,
                        metavar = '<XYZFILE>',
                        help = 'file containing the starting structure')



    return parser.parse_args()



if __name__ == '__main__':
    main()
