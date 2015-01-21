import h5py
import numpy as np
import scipy.stats as stats
import pylab as pl

# Set some physical constants
gs   = -0.5 * -1.0;
muB  = 9.27400915e-24;
mRb  = 1.443160648e-25;
pi   = 3.14159265;
kB   = 1.3806503e-23;
hbar = 1.05457148e-34;

nAtoms = 1e6;

B0 = 0.01;
dBdx = 20.;
d2Bdx2 = 40000;

Fn = nAtoms / ntrials;

Ek = np.zeros((tres))
Ep = np.zeros((tres))
Et = np.zeros((tres))
Temp = np.zeros((N.size,))

Bx = dBdx*pos[:,0,:] - 0.5*d2Bdx2*pos[:,0,:]*pos[:,2,:];
By =-dBdx*pos[:,1,:] - 0.5*d2Bdx2*pos[:,1,:]*pos[:,2,:];
Bz = B0 + 0.5*d2Bdx2*( pos[:,2,:]*pos[:,2,:] - 0.5*( pos[:,0,:]*pos[:,0,:] + pos[:,1,:]*pos[:,1,:] ) );
B  = np.sqrt( Bx**2 + By**2 + Bz**2 );

Bnx = Bx / B;
Bny = By / B;
Bnz = Bz / B;

for i in range(0,tres):
    kinetic = 0.5 * mRb * np.sum(vel[:,:,i]**2, 1)
    n = atomID[0:N[i],0,i].astype(int)
    Ek[i] = np.sum( kinetic[n], 0 ) / N[i] / kB * 1.e6
    proj = 2. * Bnx[n,i] * ( psiUp[n,0,i]*psiDn[n,0,i] + psiUp[n,1,i]*psiDn[n,1,i] ) + \
           2. * Bny[n,i] * ( psiUp[n,0,i]*psiDn[n,1,i] - psiUp[n,1,i]*psiDn[n,0,i] ) + \
           2. * Bnz[n,i] * ( psiUp[n,0,i]**2 + psiUp[n,1,i]**2 - 0.5 )
    Ep[i] = np.sum( 0.5*gs*muB*(B[n,i]*proj - B0), 0 ) / N[i] / kB * 1.e6
    Et[i] = Ek[i] + Ep[i]

    Temp[i] = 2./3. * np.sum( kinetic[n], 0) / N[i] / kB * 1.e6

norm = psiUp[:,0,:]**2 + psiUp[:,1,:]**2 + psiDn[:,0,:]**2 + psiDn[:,1,:]**2
dNorm = np.mean(norm - 1., 0)

fig = pl.figure(1)

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