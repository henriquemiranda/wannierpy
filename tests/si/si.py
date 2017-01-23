#Copyright (c) 2016, Henrique Miranda
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
import argparse

parser = argparse.ArgumentParser(description='Calculate wannier bands.')
parser.add_argument('-s',  '--scf',       action="store_true", help='run scf calculation')
parser.add_argument('-n',  '--nscf',      action="store_true", help='run nscf calculation')
parser.add_argument('-wp', '--wannierpp', action="store_true", help='run wannier preprocessing')
parser.add_argument('-pw', '--pw2wann',    action="store_true", help='run pw2wann')
parser.add_argument('-w',  '--wannier',   action="store_true", help='run wannier')
parser.set_defaults(run_all=True)
args = parser.parse_args()

w = Wannier('si',12,8,[4,4,4],[6,6,6],{'Si':"s,px,py,pz"},
            [('hr_plot','true'),
             ('num_iter','1000'),
             ('num_print_cycles','200'),
             ('dis_num_iter','500'),
             ('dis_froz_max','6.0'),
             ('dis_win_max','17.5')])

if args.scf:       w.run_gs()
if args.nscf:      w.run_nscf()
if args.wannierpp: w.run_wannierpp()
if args.pw2wann:   w.run_pw2wann()
if args.wannier:   w.run_wannier()
