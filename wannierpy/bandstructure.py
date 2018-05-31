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

class Bandstructure():
    def __init__(self,kpath,divisions):
        self.points,self.kpath = list(zip(*kpath))
        self.kpath = np.array(self.kpath)
        self.divisions = np.array(divisions)
        self.fermi = 0
        self.gen_kpath()
 
    def gen_kpath(self):
        """ Generate the path in reduced coordinates
        """
        kpoints = []
        for i in range(len(self.divisions)):
            for j in range(self.divisions[i]):
                kpoints.append(self.kpath[i] + float(j)/self.divisions[i]*(self.kpath[i+1]-self.kpath[i]))
        kpoints.append(self.kpath[-1])
        self.kpoints = kpoints
        return kpoints

    def set_fermi(self,fermi):
        self.fermi = fermi
   
    def estimate_band_connection(self, prev_eigvecs, eigvecs, prev_band_order):
        """ this function is used to estimate the bands connection
        taken from phonopy: http://atztogo.github.io/phonopy/ 
        """
        metric = np.abs(np.dot(prev_eigvecs.conjugate().T, eigvecs))
        connection_order = []
        indices = list(range(len(metric)))
        indices.reverse()
        for overlaps in metric:
            maxval = 0
            for i in indices:
                val = overlaps[i]
                if i in connection_order:
                    continue
                if val > maxval:
                    maxval = val
                    maxindex = i
            connection_order.append(maxindex)

        band_order = [connection_order[x] for x in prev_band_order]

        return band_order

    def cut_bands(self,band_min=None,band_max=None):
        """ A function to change the number of bands
        """
        if not band_min: band_min = 0
        if not band_max: band_max = len(self.eigenvalues[0,:])
        band_max += 1 #include the last band
        self.nbands = band_max-band_min
        self.eigenvalues = self.eigenvalues[:,band_min:band_max]

    def order_bands(self):
        """ with this routine we order the bands according to their caracter
        """
        order = list(range(self.nbands))
        for k in range(1,self.nkpoints):
            order = self.estimate_band_connection(self.eigenvectors[k-1],self.eigenvectors[k],order)
            eig = [self.eigenvalues[k,i] for i in order]
            self.eigenvalues[k,:] = eig

    def write_netcdf(self,filename,atoms):
        """ Write the bandstructure in a netcdf format
            filename -> output filename
            lat -> lattice in bohr
            red -> atom positions in reduced coordinates
            atypes -> atom species
        """
        try:
            from netCDF4 import Dataset
        except:
            print("Error importing netcdf")
            exit()

        natoms = len(atoms)
        ncfile = Dataset(filename, 'w')
        ncfile.createDimension('complex', 2)
        ncfile.createDimension('number_of_cartesian_dimensions', 3)
        ncfile.createDimension('number_of_reduced_dimensions', 3)
        ncfile.createDimension('number_of_kpoints', self.nkpoints)
        ncfile.createDimension('number_of_atoms', natoms)
        ncfile.createDimension('number_of_bands', self.nbands)

        primvecs_nc  = ncfile.createVariable('primitive_vectors','f4',('number_of_cartesian_dimensions','number_of_cartesian_dimensions'))
        kpoints_nc   = ncfile.createVariable('kpoints','f8',('number_of_kpoints','number_of_reduced_dimensions'))
        eig_nc       = ncfile.createVariable('eigenvalues','f8',('number_of_kpoints','number_of_bands'))
        atoms_pos_nc = ncfile.createVariable('reduced_atom_positions','f4',('number_of_atoms','number_of_cartesian_dimensions'))
        atypes_nc    = ncfile.createVariable('atom_species','i4',('number_of_atoms'))
        weights_nc   = ncfile.createVariable('weights','f8',('number_of_kpoints'))

        #ctructural variables
        primvecs_nc[:] = atoms.get_cell()
        atoms_pos_nc[:]   = atoms.get_scaled_positions()
        atypes_nc[:] = atoms.get_atomic_numbers()

        #bandstrcture variables
        kpoints_nc[:]  = self.kpoints
        eig_nc[:] = self.eigenvalues-self.fermi
        weights_nc[:] = [1.0] * self.nkpoints

        ncfile.close()

    def plot_bandstructure(self,ylim=None,title=None,filename=None,show=True):
        """ plot the bandstructure using python 
        """
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',serif="Computer Modern Roman",size=20)
        x = list(range(self.nkpoints))
        for n in range(self.nbands): 
            plt.plot(x,self.eigenvalues[:,n]-self.fermi)
        plt.ylabel("Energy [eV]")

        #plot the labels
        x = 0
        pos = []
        for l,xd in zip(self.points,self.divisions):
            pos.append(x)
            x += xd
            plt.axvline(x,color='k')
        pos.append(x)       
 
        plt.xticks(pos,self.points)
        plt.xlim([0,self.nkpoints])
        if ylim:  plt.ylim(ylim)
        if title: plt.title(title)
        if filename: plt.savefig(filename,transparent=True)
        if show: plt.show()
        return plt

    def unfold(self,indexes,vectors):
        """ Unfold a bandstructure of a supercell into a unit cell
        """
        assert self.nbands == len(indexes) == len(vectors), "The number of indexes does not correspond to the number of orbitals"
       
        #get eigenvectors
        f = open('unfold.dat','w')
        norbitals = max(indexes)+1
        weigths = np.zeros([norbitals],dtype=complex)
        for nk,k in enumerate(self.kpoints):
            eiv = self.eigenvectors[nk]
            for n in range(self.nbands):
                weigths[:] = 0
                for o,v in enumerate(vectors):
                    weigths[indexes[o]] += eiv[n,o]*exp(I*2*pi*np.dot(k,v))
                val = abs(np.dot(weigths.conj(),weigths))*4
                f.write("%d %lf %lf\n"%(nk,self.eigenvalues[nk,n],val))
        f.close()

    def get_bandstructure_eivecs(self,eig):
        """ get the bnsdstructure and eigenvectors
        """        
        kpoints = self.kpoints
        eigenvalues = []
        eigenvectors = []
        for k in kpoints:
            eigs, eivs = eig(k)
            eigenvalues.append(  eigs )
            eigenvectors.append( eivs )
        self.eigenvalues  = np.array(eigenvalues)
        self.eigenvectors = np.array(eigenvectors)
        self.nkpoints, self.nbands = self.eigenvalues.shape
        return self.eigenvalues, self.eigenvectors

    def get_bandstructure(self,eig):
        """ get the bandstructure in an array.
        eig is a fuction that takes as input the kpoint and returns the eigenvalues
        """
        kpoints = self.kpoints
        eigenvalues = []
        for k in kpoints:
            eigenvalues.append(eig(k))
        self.eigenvalues = np.array(eigenvalues)
        self.nkpoints, self.nbands = self.eigenvalues.shape
        return self.eigenvalues

    def write_bandstructure(self,filename):
        """ write the bandstructure in a data file to be used by gnuplot
        """
        f = open(filename,'w')
        for n in range(self.nbands):
            for k in range(self.nkpoints):
                f.write("%3d %12.8lf\n"%(k,self.eigenvalues[k,n]))
            f.write("\n\n")
        f.close()


