#!/usr/bin/python

import os, sys
import numpy as np
from elements import ELEMENTS

class io(object):
    @classmethod
    def readxyz(self, fname):

        atomnos = []
        atomcoords = []
        linecount = 0
        for line in open(fname).readlines():
            linecount += 1
            words = line.split()
            if linecount == 2 : cstring = line
            if len(words) == 4 and linecount > 2:
                atomnos.append(words[0].upper())
                atomcoords.append(map(float, words[1:]))
        return atomnos, atomcoords,cstring

    def writexyz(self, fname, nxyzs, cstring):

        f = open(fname, 'w')

        f.write('%i\n  %s\n'%(len(nxyzs), cstring.rstrip()))
        for nxyz in nxyzs:
            symb = ELEMENTS[nxyz[0]].symbol
            f.write("%2s%12.6f%12.6f%12.6f\n"%(symb,nxyz[1],nxyz[2],nxyz[3]))
        f.close()

        return

class aux(object):
    
    @classmethod
    def intersect(self, list1, list2):
        result=[]
        for x in list1:
            if x in list2:
                result.append(x)

        return result

    @classmethod
    def union(self,list1, list2):
        result = list1
        for x in list2:
            if x not in result:
                result.append(x)

        return result

class fragKern(object):
    factor = 1.05
    
    @classmethod
    def neighbours(self, atomnos, atomcoords):

        neig = {}
        for at1 in xrange(len(atomnos)):
            neig[at1] = []
            xyz1 = np.array(atomcoords[at1])
            rad1 = ELEMENTS[atomnos[at1]].covrad
            for at2 in range(len(atomnos)):
                if at2 == at1 :
                    continue

                xyz2 = np.array(atomcoords[at2])
                rad2 = ELEMENTS[atomnos[at2]].covrad
                d = np.sqrt(np.dot(xyz1 - xyz2, xyz1 - xyz2))
                if  d < self.factor * (rad1 + rad2):
                    neig[at1].append(at2)

        return neig

    @classmethod
    def molsplit(self, neighbours):
        cluster={}
        cl = 0
        keys = neighbours.keys()

        while len(keys) > 0:
            at1 = keys[0]
            if at1 in neighbours.keys():
                cluster[cl] = [at1]
                aux.union(cluster[cl],neighbours[at1])
                for at2 in cluster[cl]:
                    aux.union(cluster[cl],neighbours[at2])
                    del neighbours[at2]
            keys = neighbours.keys()
            cl += 1

        return cluster


class structparams(object):

    @classmethod
    def molecular_mass(self, atom_indexes,atomnos):
        mass = 0.
        for at in atom_indexes:
            mass += ELEMENTS[atomnos[at]].mass
        return mass

    @classmethod
    def baricenter(self, atom_indexes,atomnos,atomcoords):
        atomcoords_a = []
        for at in atomcoords:
            atomcoords_a.append(np.array(at))

        b = np.zeros(3)
        for at in atom_indexes:
            b += atomcoords_a[at]*ELEMENTS[atomnos[at]].mass
        return b/self.molecular_mass(atom_indexes,atomnos)


if __name__ == '__main__':

    fname = sys.argv[1]
    bname = os.path.splitext(fname)[0]

    atomnos, atomcoords, cstring = io.readxyz(sys.argv[1])

    neighbours = fragKern.neighbours(atomnos, atomcoords)
#    print voisins

    monomers = fragKern.molsplit(neighbours)
    
    for m in monomers.keys():
        monomer = []
        for atom in monomers[m]:
            monomer.append([atomnos[atom]]+atomcoords[atom])
        fname = '%s_mono%i.xyz'%(bname, m+1)
#        writexyz(fname, monomer, cstring)

    print "Mon 1 molecular mass: ",structparams.molecular_mass(monomers[0],atomnos)
    print "Mon 2 molecular mass: ",structparams.molecular_mass(monomers[1],atomnos)

    print "Mon 1: ", sorted(monomers[0])
    print "Mon 2: ", sorted(monomers[1])
    
    mon1_b = structparams.baricenter(monomers[0],atomnos,atomcoords)
    mon2_b = structparams.baricenter(monomers[1],atomnos,atomcoords)
    print 'Monomer Distance: %12.4f' % np.linalg.norm(mon1_b - mon2_b)
    #print type(baricenter(range(len(atomnos)),atomnos,atomcoords))
    
    
