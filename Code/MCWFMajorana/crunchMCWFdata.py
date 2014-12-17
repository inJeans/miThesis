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

dBdz = 2.5;

Ek = np.zeros((tres))
Ep = np.zeros((tres))
Et = np.zeros((tres))

spinUp = np.zeros((tres));
spinDn = np.zeros((tres));

xmcwf = np.zeros((tres));
vmcwf = np.zeros((tres));

xup = np.zeros((tres));
vup = np.zeros((tres));
xdn = np.zeros((tres));
vdn = np.zeros((tres));

Bx = 1.e-6;
By = 0.0;
Bz =-1.0 * dBdz * pos[:,2,:];
B  = np.sqrt( Bx**2 + By**2 + Bz**2 );

Bnx = Bx / B;
Bny = By / B;
Bnz = Bz / B;

for i in range(0,tres):
    kinetic = 0.5 * mRb * np.sum(vel[:,:,i]**2, 1)
    n = atomID[0:ntrials,0,i].astype(int)
    Ek[i] = np.mean( kinetic[n], 0 ) / kB * 1.e6
    Ep[i] = np.mean( (2.*atomIsSpinUp[n,0,i]-1.)*0.5*gs*muB*B[n,i], 0 ) / kB * 1.e6
    Et[i] = Ek[i] + Ep[i]

    spinUp[i] = np.mean( atomIsSpinUp[n,0,i], 0 )
    spinDn[i] = 1 - spinUp[i]

    xmcwf[i] = np.mean( pos[n,2,i] )*1.e6
    vmcwf[i] = np.mean( vel[n,2,i] )*1.e3

    if time[i] < 1.25:
        nup = np.where( atomIsSpinUp[:,0,i] )

    xup[i] = np.mean( pos[nup,2,i] )*1.e6
    vup[i] = np.mean( vel[nup,2,i] )*1.e3

    if time[i] < 1.25:
        ndn = np.where( 1 - atomIsSpinUp[:,0,i] )
    
    xdn[i] = np.mean( pos[ndn,2,i] )*1.e6
    vdn[i] = np.mean( vel[ndn,2,i] )*1.e3
                       
Eht = Ehk + Ehp