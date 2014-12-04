import h5py
import numpy as np
import scipy.stats as stats
import pylab as pl

# These are the "Tableau 20" colors as RGB.
#tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
#             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
#             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
#             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
#             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)];

tableau20 = [ '#1F77B4', '#AEC7E8', '#FF7F0E', '#FFBB78', '#2CA02C', '#98DF8A' ];

# Set some physical constants
gs   = -0.5 * -1.0;
muB  = 9.27400915e-24;
mRb  = 1.443160648e-25;
pi   = 3.14159265;
kB   = 1.3806503e-23;
hbar = 1.05457148e-34;
T    = 20.e-6;
B0 = 0.01;
dBdx = 20.;
d2Bdx2 = 40000;

tres = 26;
nAtoms = 1e6;
ntrials = 1e4;
nCells = 100**3 + 1;
dt = 1e-6;

Fn = nAtoms / ntrials;

time = np.zeros((tres));
pos = np.zeros((ntrials,3,tres));
vel = np.zeros((ntrials,3,tres));
psiUp = np.zeros((ntrials,2,tres));
psiDn = np.zeros((ntrials,2,tres));
atomID = np.zeros((ntrials,1,tres));
N = np.zeros((tres));
collisionCount = np.zeros((nCells,1,tres));

f = h5py.File('outputData.h5');

dset = f.require_dataset('atomData/simuatedTime',(1,1,tres),False,False);
dset.read_direct(time);

dset = f.require_dataset('atomData/positions',(ntrials,3,tres),False,False);
dset.read_direct(pos);

dset = f.require_dataset('atomData/velocities',(ntrials,3,tres),False,False);
dset.read_direct(vel);

dset = f.require_dataset('atomData/psiUp',(ntrials,2,tres),False,False);
dset.read_direct(psiUp);

dset = f.require_dataset('atomData/psiDn',(ntrials,2,tres),False,False);
dset.read_direct(psiDn);

dset = f.require_dataset('atomData/atomID',(ntrials,1,tres),False,False);
dset.read_direct(atomID);

dset = f.require_dataset('atomData/atomNumber',(1,1,tres),False,False);
dset.read_direct(N);

dset = f.require_dataset('atomData/collisionCount',(nCells,1,tres),False,False);
dset.read_direct(collisionCount);

f.close()

totalColl = np.sum( collisionCount, 0 )[0,:];
collRate = totalColl / ( time[1]-time[0] ) / N ;
tau = np.mean(collRate[10:-1])

time = time * tau;

Ek = np.zeros((N.size,))
Ep = np.zeros((N.size,))
Et = np.zeros((N.size,))
Temp = np.zeros((N.size,))

Bx = dBdx*pos[:,0,:] - 0.5*d2Bdx2*pos[:,0,:]*pos[:,2,:];
By =-dBdx*pos[:,1,:] - 0.5*d2Bdx2*pos[:,1,:]*pos[:,2,:];
Bz = B0 + 0.5*d2Bdx2*( pos[:,2,:]*pos[:,2,:] - 0.5*( pos[:,0,:]*pos[:,0,:] + pos[:,1,:]*pos[:,1,:] ) );
B  = np.sqrt( Bx**2 + By**2 + Bz**2 );

Bnx = Bx / B;
Bny = By / B;
Bnz = Bz / B;

for i in range(0,N.size):
    kinetic = 0.5 * mRb * np.sum(vel[:,:,i]**2, 1)
    n = atomID[0:N[i],0,i].astype(int)
    Ek[i] = np.sum( kinetic[n], 0 ) / N[i] / kB * 1.e6
    proj = 2. * Bnx[n,i] * ( psiUp[n,0,i]*psiDn[n,0,i] + psiUp[n,1,i]*psiDn[n,1,i] ) + \
           2. * Bny[n,i] * ( psiUp[n,0,i]*psiDn[n,1,i] - psiUp[n,1,i]*psiDn[n,0,i] ) + \
           2. * Bnz[n,i] * ( psiUp[n,0,i]**2 + psiUp[n,1,i]**2 - 0.5 )
    Ep[i] = np.sum( 0.5*gs*muB*(B[n,i]*proj - B0), 0 ) / N[i] / kB * 1.e6
    Et[i] = Ek[i] + Ep[i]
    
    Temp[i] = 2./3. * np.sum( kinetic[n], 0) / N[i] / kB * 1.e6

pl.figure(1)
pl.subplot(3,1,1)
pl.plot(time,(Et-Et[0])/Et[0]*100.,
        lw=2.5, color=tableau20[3])
pl.ylabel(r'$\%\Delta E_T$')
ax = pl.gca()
ax.locator_params('y',tight=True, nbins=4)
pl.grid()

norm = psiUp[:,0,:]**2 + psiUp[:,1,:]**2 + psiDn[:,0,:]**2 + psiDn[:,1,:]**2
dNorm = np.mean(norm - 1., 0)

pl.subplot(3,1,2)
pl.plot(time,dNorm*100.,
        lw=2.5, color=tableau20[4]);
pl.ylabel(r'$\%\Delta \vert\Psi\vert^2$')
ax = pl.gca()
ax.locator_params('y',tight=True, nbins=4)
pl.grid()

pl.subplot(3,1,3)
pl.plot(time,N*Fn,
        lw=2.5, color=tableau20[5]);
pl.ylabel(r'$N$')
ax = pl.gca()
ax.locator_params('y',tight=True, nbins=1)
pl.xlabel(r'$t \times \tau_{c}$')
pl.grid()

pl.savefig("/Users/miMac/Documents/versionControlledFiles/miThesis/gfx/Ehrenfest/ehrenfestIPConserve.eps")

fig = pl.figure(2)
n = atomID[0:N[0],0,0].astype(int)
ri = np.sqrt( 0.5*pos[n,0,0]**2 + 0.5*pos[n,1,0]**2 + pos[n,2,0]**2 );

nri, binsri, patches = pl.hist(ri,100)
nri = np.append([0], nri , axis=0)

n = atomID[0:N[-1],0,-1].astype(int)
rf = np.sqrt( 0.5*pos[n,0,-1]**2 + 0.5*pos[n,1,-1]**2 + pos[n,2,-1]**2 );

nrf, binsrf, patches = pl.hist(rf,100)
nrf = np.append([0], nrf , axis=0)

maxwell = binsri**2 * np.exp(-0.25*gs*muB*d2Bdx2*binsri**2 / (kB*Temp[0]*1.e-6))
maxwell = maxwell / maxwell.max() * nri.max()

pl.close(fig)

pl.figure(3)
pl.plot( binsri, nri*Fn, label=r'$n(\mathbf{r},t_\mathrm{inital})$',
        lw=2.5, color=tableau20[0] )
pl.plot( binsrf, nrf*Fn, label=r'$n(\mathbf{r},t_\mathrm{final})$',
        lw=2.5, color=tableau20[1] )
pl.plot( binsri, maxwell*Fn, label=r'$n_\mathrm{IP}(\mathbf{r})$',
        lw=2.5, color=tableau20[2] )
pl.ylabel(r'$n(\mathbf{r})$')
pl.xlabel(r'$\vert\mathbf{r}\vert$')
ax = pl.gca()
ax.locator_params('x',tight=True, nbins=6)
pl.legend(loc='best')
pl.grid()

pl.savefig("/Users/miMac/Documents/versionControlledFiles/miThesis/gfx/Ehrenfest/ehrenfestIPDist.eps")

pl.show()