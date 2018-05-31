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

from qepy import *
from itertools import product
import argparse

def boolean(boolean):
    return [".false.",".true."][boolean]

class Pw2wann():
    _pw2wann = 'pw2wannier90.x'

    def __init__(self,prefix,outdir='./',amn=True,unk=True,mmn=True,filename=None):
        self.outdir = outdir
        self.prefix = prefix
        self.seedname = prefix
        self.write_amn  = amn
        self.write_mmn  = mmn
        self.write_unk  = unk
        if filename == None:
            self.filename = "%s.pw2wann"%self.prefix
        else:
            self.filename = filename

    def write(self):
        f = open(self.filename,'w')
        f.write(str(self))
        f.close()

    def __str__(self):
        s = "&inputpp\n"
        s += "%10s = '%s'\n"%("outdir",self.outdir)
        s += "%10s = '%s'\n"%("prefix",self.prefix)
        s += "%10s = '%s'\n"%("seedname",self.seedname)
        s += "%10s = %s\n"  %("write_amn",boolean(self.write_amn))
        s += "%10s = %s\n"  %("write_mmn",boolean(self.write_mmn))
        s += "%10s = %s\n"  %("write_unk",boolean(self.write_unk))
        s += "/\n"
        return s
    
    def run(self):
        self.write()
        command = '%s < %s > %s.log'%(self._pw2wann,self.filename,self.prefix)
        print(command)
        os.system(command)
 
class Wannier():
    _scf = 'scf'
    _nscf = 'nscf'
    _wannier90 = 'wannier90.x'

    def __init__(self,prefix,num_bands,num_wann,scf_kpoints,nscf_kpoints,projections,options=None,dry=False,filename=None):
        self.dry = dry
        self.projections = projections
        self.prefix = prefix
        self.qe_input = qein(prefix+'.scf')
        self.cell = self.qe_input.cell_parameters
        self.atoms = self.qe_input.atoms
        self.scf_kpoints = scf_kpoints
        self.nscf_kpoints = nscf_kpoints
        self.num_bands = num_bands
        self.num_wann = num_wann
        self.options = options
        if filename == None:
            self.filename = "%s.win"%self.prefix
        else:
            self.filename = filename
        self.gen_kpoints()

    def run_gs(self):
        os.system('mkdir -p %s'%self._scf)
        self.qe_input.kpoints = self.scf_kpoints
        filename = "%s/%s.scf"%(self._scf,self.prefix)
        if self.dry:
            self.qe_input.write(filename)
        else:
            self.qe_input.run(filename)


    def gen_kpoints(self):
        """ create a list of kpoints in the full k-mesh
        """
        kpts = self.nscf_kpoints
        nkpts = kpts[0]*kpts[1]*kpts[2]
        self.klist = []
        for i,j,k in product(list(range(kpts[0])),list(range(kpts[1])),list(range(kpts[2]))):
            self.klist.append([float(i)/kpts[0],float(j)/kpts[1],float(k)/kpts[2],1.0/nkpts])
    
    def run_nscf(self):
        #check if teh scf calculation is present
        scf_save  = "%s/%s.save"%(self._scf,self.prefix)
        nscf_save = "%s/%s.save"%(self._nscf,self.prefix)
        if os.path.isdir(scf_save):
            os.system('mkdir -p %s'%self._nscf)
            os.system('cp -r %s %s'%(scf_save,nscf_save))
        else:
            print("scf calculation not found!")
            exit()
    
        self.qe_input.control['calculation'] = "'nscf'"
        self.qe_input.system['force_symmorphic'] = '.true.'
        self.qe_input.system['nbnd'] = self.num_bands
        self.qe_input.ktype = 'crystal'
        self.qe_input.klist = self.klist
        filename = "%s/%s.nscf"%(self._nscf,self.prefix)
        #create the input file
        if self.dry:
            self.qe_input.write(filename)
        else:
            self.qe_input.run(filename)

    def write(self):
        f = open(self.filename,"w")
        f.write(str(self))
        f.close()

    def __str__(self):
        s = ""
        s += "%10s = %5d\n"%("num_wann",self.num_wann)
        s += "%10s = %5d\n"%("num_bands",self.num_bands)
        s += ("%10s = "+"%4d "*3+"\n")%tuple(["mp_grid"]+self.nscf_kpoints)
        #print options
        if self.options:
            for option,value in self.options:
                s += "%10s = %10s\n"%(option,value)
        s += "\nbegin unit_cell_cart\n"
        s += "bohr\n"
        #write the cell parameters
        for vec in self.cell:
            s += ("%12.8lf"*3+"\n")%tuple(vec)
        s += "end unit_cell_cart\n"
        #write the projections
        s += "begin projections\n"
        for atom,proj in list(self.projections.items()):
            s+= "%s: %s\n"%(atom,proj)
        s += "end projections\n"
        #write the postiions of the atoms
        s += "\nbegin atoms_cart\n"
        for atom in self.atoms:
            s += ("%2s "+"%12.8lf"*3+"\n")%tuple([atom[0]]+atom[1])
        s += "end atoms_cart\n"
        #write the klist
        s += "\nbegin kpoints\n"
        for k in self.klist:
            s += ("%12.8lf"*3+"\n")%tuple(k[:3])
        s += "end kpoints\n"
        return s

    def run_wannierpp(self):
        self.write()
        command = "%s -pp %s"%(self._wannier90,self.filename)
        print(command)
        if self.dry: print(command)
        else:   os.system(command)

    def run_pw2wann(self):
        self.write()
        pw2wann = Pw2wann(self.prefix,outdir="./%s"%self._nscf,unk=False)
        if self.dry: pw2wann.write()
        else:        pw2wann.run()

    def run_wannier(self):
        self.write()
        command = "%s %s"%(self._wannier90,self.prefix)
        print(command)
        os.system(command)
