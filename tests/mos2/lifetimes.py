from wannierpy import *

np.set_printoptions(precision=3,linewidth=150,suppress=True)

kmesh = [12,12,1]
kmesh_dos = [50,50,1]
debug = False
broad=2.0

#kpath
G = [ 0.0,  0.0,  0.0] # Gamma point
M = [ 0.5,  0.0,  0.0] # M point
K = [2./3,-1./3,  0.0] # K point
gri = [100,50,100]

def regular_mesh(nx,ny,nz):
    kmesh = []
    for i,j,k in product(xrange(nx),xrange(ny),xrange(nz)):
        kmesh.append([float(i)/nx,float(j)/ny,float(k)/nz])
    #kmesh = np.array(kmesh) - np.array([0.5,0.5,0.5])
    return np.array(kmesh)

def histogram_eig(eiv,emin=-5.0,emax=5.0,step=0.01,sigma=0.05):
    """
    Histogram of eigenvalues using lorentzians
    """
    eiv = np.array(eiv)
    #sigma = 0.005
    x = np.arange(emin,emax,step,dtype=np.float32)
    y = np.zeros([len(x)],dtype=np.float32)

    #lorentzian stuff
    s2 = (.5*sigma)**2
    c = (.5*sigma)

    eiv     = eiv.flatten()
    eiv     = eiv[emin < eiv]
    eiv     = eiv[eiv < emax]
    
    for e in eiv:
        x1 = (x-e)**2
        y += c/(x1+s2)
    return x, y

def plot_bandstructure(filename,w):
    kw = np.array([G,M,K,G])
    gri = [100,50,100]

    #generate reduced coordinates
    k_red = []
    for i in xrange(len(gri)):
        for j in xrange(gri[i]):
            k_red.append(kw[i] + float(j)/gri[i]*(kw[i+1]-kw[i]))

    f = open(filename,"w")
    for nk,eigk in enumerate([w.get_eigvals(k) for k in k_red]):
        for eig in eigk:
            f.write(("%d %12.8lf %12.8lf\n")%(nk,eig.real,eig.imag))
        f.write("\n\n")
    f.close()

#read wannier hamiltonian on real space 6x6 mesh
w = Wannier_hr('mos2')


#plot bandstructure
plot_bandstructure('read.dat',w)

#calculate eigenvalues on a regular mesh 6x6 kpoints
full_kmesh = regular_mesh(*kmesh)
nkpoints = len(full_kmesh)
nbands = w.nwann
eig = np.zeros([nkpoints,nbands])
eiv = np.zeros([nkpoints,nbands,nbands],dtype=np.complex64)
h   = np.zeros([nkpoints,nbands,nbands],dtype=np.complex64)
for nk,k in enumerate(full_kmesh):
    h[nk] = w.get_ham_kpoint(kpoint=k)
    eig[nk], eiv[nk] = np.linalg.eigh(h[nk])

#calculate dos
kmesh_dos = regular_mesh(*kmesh_dos)
nkpoints_dos = len(kmesh_dos)
eig_dos = np.zeros([nkpoints_dos,nbands])
for nk,k in enumerate(kmesh_dos):
    hi = w.get_ham_kpoint(kpoint=k)
    eig_dos[nk] = np.linalg.eigvalsh(hi) 

d = 0.1
emin = np.min(eig_dos)-d
emax = np.max(eig_dos)+d
energies, dos = histogram_eig(eig_dos,emin=emin,emax=emax)

#normalize dos to 1
dos = dos/np.max(dos)*broad

if debug:
    import matplotlib.pyplot as plt
    plt.plot(energies,dos)
    plt.show()

#calculate lifetimes for the points on the regular mesh
from scipy.interpolate import interp1d
f = interp1d(energies, dos)
#lifetimes = np.array([[f(eig) for eig in eigk ] for eigk in eig])

#get ordered eigenvalues
eig_ordered = [ np.diagonal(np.dot(np.dot(np.conjugate(u.T),hi),u)) for hi,u in zip(h,eiv) ]
eig_ordered = np.array([ [e+f(e)*broad*complex(0,1) for e in eigk] for eigk in eig_ordered ])
#eig_ordered = np.array([ [e for e in eigk] for eigk in eig_ordered ])

#rebuild the hamiltonian with the different eigenvalues
h = [ np.dot(np.dot(u,np.eye(nbands)*e),np.conjugate(u.T)) for e,u in zip(eig_ordered,eiv) ] 

#replace eigenvalues on a regular mesh by eigenvalues +lifetimes
w.ham = w.get_ham_real(h,full_kmesh)

#plot band structure
plot_bandstructure('final.dat',w)

#fold back to real space and write in a new *_hr.dat file
#interpolate to a 24x24 mesh
#plot bandstructure with lifetimes
