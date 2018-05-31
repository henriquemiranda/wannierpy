#Copyright (c) 2015, Henrique Miranda
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
#    * Neither the name of the wannierpy project nor the
#      names of its contributors may be used to endorse or promote products
#      derived from this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
#DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
#ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from itertools import product, repeat
import numpy as np
from cmath import exp,pi
import pickle
import os
from matplotlib import pyplot as plt
from math import ceil
I = complex(0,1)

def red_car(red,lat):
    """
    Convert reduced coordinates to cartesian
    """
    return np.array(map( lambda coord: coord[0]*lat[0]+coord[1]*lat[1]+coord[2]*lat[2], red))

def car_red(car,lat):
    """
    Convert cartesian coordinates to reduced
    """
    return np.array(map( lambda coord: np.linalg.solve(np.array(lat).T,coord), car))

def plot_matrix(mat):
    from pygnuplot import Pygnuplot
    nx,ny = np.array(mat).shape
    p = Pygnuplot(xrange="[-.5:%lf]"%(nx-0.5),yrange='[-.5:%lf]'%(ny-0.5))
    p.writeline('unset key')
    p.addplot(mat,pwith='image',cmd='array=(%d,%d) flip=y'%(nx,ny),nx=nx,ny=ny)
    p.draw()

class Wannier_win():
    def __init__(self,seedname):
        f = open(seedname+".win",'r')

        #read cartesian vectors
        self.lat = []
        for line in f:
            if "begin unit_cell_cart" in line:
                for line in f:
                    if "end unit_cell_cart" in line:
                        break
                    if "bohr" in line or "angstroem" in line:
                        continue
                    else:
                        self.lat.append( [float(x) for x in line.strip().split()] )

        #read atomic postions
        cart = False #flag to know if the coordinates of the atoms are cartesian 
        f.seek(0)
        self.atoms = []
        self.atypes = []
        for line in f:
            if "begin atoms_cart" in line or  "begin atoms_frac" in line:
                if 'cart' in line: cart = True
                for line in f:
                    if "end" in line:
                        break
                    if "bohr" in line or "angstroem" in line:
                        continue
                    else:
                        line_strip = line.strip().split()
                        atom = [line_strip[0], [ float(x) for x in line_strip[1:] ]]
                        if cart:
                            atom = car_red([atom],self.lat)
                        self.atoms.append( atom )
        f.close()
        self.natoms = len(self.atoms)

    def get_atoms(self):
        """ Get the atoms from the wannier file
        """
        sym = [a[0] for a in self.atoms]
        pos = [a[1] for a in self.atoms]
        return self.lat, sym, pos

    def __str__(self):
        s = ""
        s += "lattice vectors:\n"
        s += "\n".join(["%12.8lf "*3%tuple(x) for x in self.lat])+"\n"
        s += "atom positions:\n"
        s += "\n".join([ ("%4s "+"%12.8lf "*3)%(x[0],x[1][0],x[1][1],x[1][2]) for x in self.atoms])
        return s


