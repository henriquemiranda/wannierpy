from wannierpy import *

kmesh = [12,12,1]
kmesh_dos = [50,50,1]
debug = False
broad=2.0

#kpath
G = [ 0.0,  0.0,  0.0] # Gamma point
M = [ 0.5,  0.0,  0.0] # M point
K = [2./3,-1./3,  0.0] # K point
gri = [500,200,500]

def plot_bandstructure(filename,w,real=True):
    kw = np.array([G,M,K,G])

    #generate reduced coordinates
    k_red = []
    for i in range(len(gri)):
        for j in range(gri[i]):
            k_red.append(kw[i] + float(j)/gri[i]*(kw[i+1]-kw[i]))

    f = open(filename,"w")
    if real:
        eig = np.array([w.get_eigvalsh(k) for k in k_red])
    else:
        eig = np.array([w.get_eigvals(k) for k in k_red])
    nkpoints, nbands = eig.shape
    for n in range(nbands):
        for nk in range(nkpoints):
            e = eig[nk,n]
            f.write(("%d %12.8lf %12.8lf\n")%(nk,e.real,e.imag))
        f.write("\n\n")
    f.close()

#read wannier hamiltonian on real space 6x6 mesh
w = Wannier_hr('mos2')

#plot bandstructure
plot_bandstructure('read.dat',w)
