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
from copy import copy
from matplotlib import pyplot as plt
from math import ceil
I = complex(0,1)

def red_car(red,lat):
    """
    Convert reduced coordinates to cartesian
    """
    return np.array([coord[0]*lat[0]+coord[1]*lat[1]+coord[2]*lat[2] for coord in red])

def car_red(car,lat):
    """
    Convert cartesian coordinates to reduced
    """
    return np.array([np.linalg.solve(np.array(lat).T,coord) for coord in car])

def plot_matrix(mat):
    from pygnuplot import Pygnuplot
    nx,ny = np.array(mat).shape
    p = Pygnuplot(xrange="[-.5:%lf]"%(nx-0.5),yrange='[-.5:%lf]'%(ny-0.5))
    p.writeline('unset key')
    p.addplot(mat,pwith='image',cmd='array=(%d,%d) flip=y'%(nx,ny),nx=nx,ny=ny)
    p.draw()

class Wannier_hr():
    def __init__(self,seedname,dopickle=True):
        """ Read the hamiltonian from wannier interpolation
        """
        self.kpoints_interpolate_hash = None
        self.kpoints_calc_hash = None
        self.seedname = seedname
        self.supercell = (1,1,1)

        #check for the existence of a pickle
        if os.path.isfile(seedname+'.npy') and dopickle:
            print("reading from", seedname+'.npy')
            self.read_npy()
        else:
            print("reading from", seedname+'_hr.dat')
            self.read_hr(seedname)

            #make a pickle
            if dopickle: self.write_npy()
                
    def read_npy(self,filename=None):
        """
        read npy from the 
        """
        if filename is None: filename = self.seedname+".npy"
        f = open(filename,'rb')
        self.nwann,self.npoints,self.degeneracy,self.ham,self.points = pickle.load(f) 
        f.close()
        
    def write_npy(self,filename=None):
        """
        write npy to the hard drive
        """
        if filename is None: filename = self.seedname+".npy"
        f = open(filename,'wb')
        pickle.dump([self.nwann,self.npoints,self.degeneracy,self.ham,self.points],f)
        f.close()

    def read_hr(self,seedname):
        """
        Read the hamiltonian in real space
        """
        f = open(seedname+"_hr.dat",'r')
        f.readline()
        self.nwann = int(f.readline())
        self.npoints = int(f.readline())
        self.degeneracy = []
        for i in range(int(ceil(self.npoints/15.0))):
            self.degeneracy += [int(i) for i in f.readline().strip().split()]
        self.ham = np.zeros([self.npoints,self.nwann,self.nwann],dtype=complex)
        self.points = np.zeros([self.npoints,3],dtype=int)
        for n in range(self.npoints):
            for i,j in product(range(self.nwann),repeat=2):
                line = f.readline().strip().split()
                self.points[n] = np.array([int(k) for k in line[:3]])
                self.ham[n,i,j] = complex(float(line[-2]),float(line[-1]))
        f.close()

    def write_hr(self):
        """
        write the hamiltonian in real space
        """
        print('Writting _hr.dat file not implemented yet')
        raise NotImplementedError
        f = open("%s_hr.dat"%seedname,'w')
        f.close()
 
    def get_ham(self):
        """
        Get the hamiltonian as a big matrix
        """
        xmin = abs(min(self.points[:,0]))
        ymin = abs(min(self.points[:,1]))
        zmin = abs(min(self.points[:,2]))

        xmax = abs(max(self.points[:,0]))
        ymax = abs(max(self.points[:,1]))
        zmax = abs(max(self.points[:,2]))

        xiter = range(-xmin,xmax)
        yiter = range(-ymin,ymax)
        ziter = range(-zmin,zmax)

        xdim = len(xiter)+1
        ydim = len(yiter)+1
        zdim = len(ziter)+1
        
        print("xrange: %d %d"%(xmin,xmax))
        print("yrange: %d %d"%(ymin,ymax))
        print("zrange: %d %d"%(zmin,zmax))

        nw = self.nwann
        ham = np.zeros([xdim,nw,ydim,nw,zdim],dtype=complex)

        for n,point in enumerate(self.points):
            nx,ny,nz = point+np.array([xmin,ymin,zmin])
            ham[nx,:,ny,:,nz] = self.ham[n,:,:]

        return ham.reshape([xdim*nw,ydim*nw,zdim])

    def plot_ham(self):
        """ generate the hamiltonian as a big matrix and plot it using pygnuplot
        """
        ham = self.get_ham()
        n = len(ham)
        plot_matrix(ham.real)

    def set_supercell(self,nx,ny,nz):
        self.supercell = [nx,ny,nz]

    def get_ham_kpoint(self,kpoint=(0,0,0)):
        nx,ny,nz = self.supercell
        nw = self.nwann

        #new ham
        hamk = np.zeros([nx,ny,nz,nw,nx,ny,nz,nw],dtype=complex)
        for xi,yi,zi in product(range(nx),range(ny),range(nz)):
            for n,point in enumerate(self.points):
                i,j,k = point

                x,y,z = xi+i, yi+j, zi+k
                atomi,atomj,atomk = x%nx,y%ny,z%nz       
                celli,cellj,cellk = x/nx,y/ny,z/nz

                phase = exp(I*2*pi*np.dot(kpoint,[celli,cellj,cellk]))
                hamk[xi,yi,zi,:,atomi,atomj,atomk,:] += self.ham[n]*phase/self.degeneracy[n]

        ham = hamk.reshape([nx*ny*nz*nw,nx*ny*nz*nw])

        return ham

    def interpolate_bands(self,kmesh_calc,wigner_mesh,kpoints_interpolate,backup=True,filename=False):
        """
        Function to interpolate a quantity in the bandstructure

        wannier     regular mesh        wannier     any kmesh
         mesh          kpoints          r-mesh
                  
           Hr    ->       Hk       ->     Hr      ->   Hk
                 1.                2.             3.
                 kmesh_calc        wigner_mesh    kpoints_interpolate

        Arguments:
        wigner_mesh   -> (points,degeneracy)
        bands_calc    -> (kpoints,eig)
        kpoints_interpolate -> kpoints
        backup        -> Do not change the hamiltonian of this instance
        filename      -> Save the new hamiltonian in a file
        """
        # 1. calculate the eigenvectors at the kpoints
        # Hr -> Hk
        kpoints_calc, eig_calc = kmesh_calc
        kpoints_calc = np.array(kpoints_calc)
        kpoints_interpolate = np.array(kpoints_interpolate)

        #hash kpoints_calc and kpoints_interpolate
        kpoints_interpolate_hash = hash(str(kpoints_interpolate))
        kpoints_calc_hash        = hash(str(kpoints_calc))

        #avoid recalculating things        
        if not self.kpoints_interpolate_hash == kpoints_interpolate_hash:
            h = self.get_ham_kpoints(kpoints_interpolate)
            _, self.eiv_interpolate = list(zip(*[ np.linalg.eigh(hk) for hk in h ]))
            self.kpoints_interpolate_hash = kpoints_interpolate_hash 
        if not self.kpoints_calc_hash == kpoints_calc_hash:
            h = self.get_ham_kpoints(kpoints_calc)
            _, self.eiv_calc = list(zip(*[ np.linalg.eigh(hk) for hk in h ]))
            self.kpoints_calc_hash = kpoints_calc_hash 

        # 2. calculate the Hr hamiltonian
        # Hk -> Hr 

        #backup
        if backup:
            bk_points = copy(self.points)
            bk_degeneracy = copy(self.degeneracy)
            bk_ham = copy(self.ham)
    
        #increase the number of real space components for the interpolation
        points, degeneracy = wigner_mesh
        self.points     = points
        self.degeneracy = degeneracy
        self.npoints    = len(points)
        nbands = self.nwann

        #replace the eigenvalues in the hamiltonian
        h = [ np.dot(np.dot(u,np.eye(nbands)*e),np.conjugate(u.T)) for e,u in zip(eig_calc,self.eiv_calc) ]

        #calculate Hr
        self.ham = self.get_ham_real(h,kpoints_calc)

        #save this hamiltonian in a file?
        if filename:
            self.write_npy(filename=filename)

        # 3. interpolate the hamiltonian at the different kpoints
        # Hr -> Hk
        h          = np.array([ self.get_ham_kpoint(kpt) for kpt in kpoints_interpolate ])
        eig_interp = np.array([ np.dot(np.dot(np.conjugate(u.T),hk),u).diagonal().real for hk,u in zip(h,self.eiv_interpolate) ])

        #restore
        if backup:
            self.points = bk_points
            self.degeneracy = bk_degeneracy
            self.ham = bk_ham
            self.npoints = len(bk_points)

        return eig_interp

    def get_ham_real(self,hamk,kpoints):
        """
        Get the hamiltonian in real space using a list of hamiltonians
        on a regular mesh of k-points
        TODO: This can be done using a FFT library
        """

        nkpoints = len(hamk)
        ham = np.zeros([self.npoints,self.nwann,self.nwann],dtype=complex)
        for nk,kpoint in enumerate(kpoints):
            for n in range(self.npoints):
                ham[n] += hamk[nk]*exp(-I*2*pi*np.dot(kpoint,self.points[n]))*self.degeneracy[n]
        return ham/nkpoints

    def get_ham_kpoint(self,kpoint=(0,0,0)):
        """
        Get the hamiltonian at a certain k-point
        TODO: This can be done using a FFT library
        """
        hamk = np.zeros([self.nwann,self.nwann],dtype=complex)
        for n in range(self.npoints):
            hamk += self.ham[n]*exp(I*2*pi*np.dot(kpoint,self.points[n]))/self.degeneracy[n]
        return hamk

    def get_hopping_matrices(self):
        """
        Get the interaction matrices in real space
        """
        xij_int = [ tuple(x) for x in self.points ]

        #create a dicitonary with the different positions and empty matrices
        return dict(list(zip(xij_int,[self.ham[i] for i in range(self.npoints)])))

    def get_ham_kpoints(self,kpoints):
        """
        Get the hamiltonian for a list of k-points
        """
        nkpoints = len(kpoints)
        nbands   = self.nwann
        h = np.zeros([nkpoints,nbands,nbands],dtype=np.complex64)
        for nk,k in enumerate(kpoints):
            h[nk] = self.get_ham_kpoint(kpoint=k)
        return h

    def get_eigvalsh(self,kpoint=(0,0,0)):
        return np.linalg.eigvalsh(self.get_ham_kpoint(kpoint))

    def get_eigvals(self,kpoint=(0,0,0)):
        return np.linalg.eigvals(self.get_ham_kpoint(kpoint))

    def get_eigh(self,kpoint=(0,0,0)):
        return np.linalg.eigh(self.get_ham_kpoint(kpoint))
 
    def __str__(self):
        s = ""
        s += "nwann: %d\n"%self.nwann
        s += "npoints: %d\n"%self.npoints
        s += "degeneracy: \n%s\n"%str(self.degeneracy)
        return s
