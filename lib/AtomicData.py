#!/usr/bin/env python2.7

import re

class AtomicData(object):

    def __init__(self, properties):
        self._atomic_numb = int(properties[0])
        self.atomic_mass = float(properties[1])
        self.atomic_name = str(properties[2])
        self.atomic_symb = str(properties[3])
        self.atomic_radius= float(properties[4])
        
    @property
    def atomic_numb(self):
        return self._atomic_numb


class AtomicDB(object):
    
    atomic_mass = []
    atomic_name = []
    atomic_symb = []
    atomic_numb = []
    atomic_radius = []
    
    init = False

    def __init__(self):

        atomicMassRe = re.compile('\s*^\s*\[?(\d+\.?\d*)\]?\(?\d?\)?\s*$')

        if not AtomicDB.init:
            AtomicDB.init = True
            self.atomic_radius=[0.35,0.28,1.28,0.96,0.84,0.76,0.71,0.66,
                         0.57,0.58,1.66,1.41,1.21,1.11,1.07,1.05,
                         1.02,1.06,2.03,1.76,1.70,1.60,1.53,1.39,
                         1.61,1.52,1.50,1.24,1.32,1.22,1.22,1.20,
                         1.19,1.20,1.20,1.16,2.20,1.95,1.90,1.75,
                         1.64,1.54,1.47,1.46,1.42,1.39,1.45,1.44,
                         1.42,1.39,1.39,1.38,1.39,1.40,2.44,2.15,
                         2.07,2.04,2.03,2.01,1.99,1.98,1.98,1.96,
                         1.94,1.92,1.92,1.89,1.90,1.87,1.87,1.75,
                         1.70,1.62,1.51,1.44,1.41,1.36,1.36,1.32,
                         1.45,1.46,1.48,1.40,1.50,1.50,2.60,2.21,
                         2.15,2.06,2.00,1.96,1.90,1.87,1.80,1.69]

            for line in open('pt-data1.csv', 'r'):
                fields = line.split(',')
                if fields[3]: 
                    self.atomic_mass.append(float(atomicMassRe.match(fields[3]).group(1)))
                else:
                    self.atomic_mass.append(float(-10))
                self.atomic_name.append(str(fields[2]).strip())
                self.atomic_symb.append(str(fields[1]).strip())
                self.atomic_numb.append(int(fields[0]))


    def find(self,arg):
        properties = [self.atomic_numb, self.atomic_mass,  self.atomic_name, self.atomic_symb, self.atomic_radius]

        if type(arg) is int:
            ind = arg - 1
            print map(lambda x: x[ind], properties)
#            print 'asdasd', [x for x in properties]



class AtFind(object):
    
    def __init__(self, arg):

        self.DB = AtomicDB()
        self.properties = [self.DB.atomic_numb, self.DB.atomic_mass,  self.DB.atomic_name, self.DB.atomic_symb, self.DB.atomic_radius]

        if type(arg) is int:
            self._selectByZ(arg)
            return 

        elif type(arg) is str:
            print 'isString'
        else:
            print arg+' unknown'

    def _selectByZ(self,z):
        ind = z -1
        return AtomicData(map(lambda x: x[ind], self.properties))
#        return AtomicData(self.DB.atomic_numb[ind], self.DB.atomic_mass[ind], self.DB.atomic_name[ind], self.DB.atomic_symb[ind], self.DB.atomic_radius[ind])


        





def findAtNumb(atnumb):
    index = atnumb - 1
    return ([atomic_name[index], atomic_mass[index], atomic_symb[index]])
