import numpy as np
import scipy.stats as stats
import pylab as pl
import os
import time

start = time.clock()

gs = 0.5;
mRb  = 1.443160648e-25;
muB = 9.27401e-24;
hbar = 1.05457e-34;
kB   = 1.3806503e-23;

dt = 1.e-7;
nsteps = np.ceil(0.003 / dt ).astype(int);

t    = np.zeros( (nsteps) );
t[0] = 0.;

x    = np.zeros( (nsteps) );
x[0] = -5.e-6;

v    = np.zeros( (nsteps) );
v[0] = 0.;

Ek    = np.zeros( (nsteps) );
Ek[0] = 0.5 * mRb * v[0]**2;

Bx = 1.0e-6;
dBdz = 2.5;

Bz =-dBdz*x[0];
magB = np.sqrt(Bx**2 + Bz**2);
nx = Bx / magB;
nz = Bz / magB;

a =-0.5 * gs * muB * dBdz * np.array( [ [-1.0, 0.],
                                        [ 0., 1.0] ] ) / mRb;
acc = 0.;

psi      = np.zeros( (2,nsteps), dtype=complex );
psi[:,0] = 0.5 * np.array( [ [1 + nx + nz],
                             [1 + nx - nz] ] )[:,0] / np.sqrt(1+nx);

Ep    = np.zeros( (nsteps) );
Ep[0] = (0.5 * gs * muB * ( Bx * (psi[0,0]*psi[1,0].conj() + psi[1,0]*psi[0,0].conj()) +\
                            Bz * (2.*psi[0,0]*psi[0,0].conj() - 1.) )).real;

Pup = 0.5*np.array( [ [1+nz, nx], [ nx, 1-nz ] ] );
Pdn = 0.5*np.array( [ [1-nz,-nx], [-nx, 1+nz ] ] );

pup    = np.zeros( (nsteps) );
pup[0] = np.dot( psi[:,0].conj().transpose(), np.dot( Pup, psi[:,0] ) ).real;
pdn    = np.zeros( (nsteps) );
pdn[0] = np.dot( psi[:,0].conj().transpose(), np.dot( Pdn, psi[:,0] ) ).real;

fliptime = 0

for i in range(1,nsteps):
    
    t[i] = t[i-1] + dt;
    
    acc = np.dot( psi[:,i-1].conj().transpose(), np.dot( a, psi[:,i-1] ) ).real;
    
    Bz =-dBdz*x[i-1];
    magB = np.sqrt(Bx**2 + Bz**2);
    nx = Bx / magB;
    nz = Bz / magB;
    
    theta = gs*muB*magB*dt / 2. / hbar;
    
    cosx = np.cos(theta);
    sinx = np.sin(theta);
    
    U = np.array( [ [cosx - 1j * nz * sinx, -1j * nx * sinx],
                    [-1j * nx * sinx, cosx + 1j * nz * sinx] ] );
        
    psi[:,i] = np.dot( U, psi[:,i-1] );
         
    Pup = 0.5*np.array( [ [1+nz, nx], [ nx, 1-nz ] ] );
    Pdn = 0.5*np.array( [ [1-nz,-nx], [-nx, 1+nz ] ] );
         
    pup[i] = np.dot( psi[:,i].conj().transpose(), np.dot( Pup, psi[:,i] ) ).real;
    pdn[i] = np.dot( psi[:,i].conj().transpose(), np.dot( Pdn, psi[:,i] ) ).real;

    v[i] = v[i-1] + acc  * dt;
    x[i] = x[i-1] + v[i] * dt;

    Bz =-dBdz*x[i];

    Ek[i] = 0.5 * mRb * v[i]**2;
    Ep[i] = (0.5 * gs * muB * ( Bx * (psi[0,i]*psi[1,i].conj() + psi[1,i]*psi[0,i].conj()) + \
                                Bz * (2.*psi[0,i]*psi[0,i].conj() - 1.) )).real;

    if (np.sign(x[i]) != np.sign(x[i-1]) ):
        fliptime = i-1
                                     
pops = (psi.conj() * psi).real;
                                     
k = gs*muB*Bx**2 / (hbar*v[fliptime]*dBdz);
                                     
print "\nProbability of spin flip is %f\n" % np.exp(-0.5*k*np.pi);
print "This simulation predicts %f\n" % pops[0,-1];

majProb = np.array( [ [np.exp(-0.5*k*np.pi)], [1.-np.exp(-0.5*k*np.pi)] ] )
majProbRot = (0.5 + -dBdz*x[fliptime]*(-0.5+majProb[0]) / np.sqrt(Bx**2 + (-dBdz*x[fliptime])**2)) * pup[fliptime];

#pl.figure(1)
#pl.plot(t*1e3, pops[0,:], label=r'$\langle\psi_{\uparrow}\vert\psi_{\uparrow}\rangle$')
#pl.plot(t*1e3, pops[1,:], label=r'$\langle\psi_{\downarrow}\vert\psi_{\downarrow}\rangle$')
#pl.plot(t*1e3, pops[0,:]+pops[1,:], label=r'$\langle\Psi\vert\Psi\rangle$')
#pl.plot([t[0]*1e3, t[-1]*1e3], [majProb[0], majProb[0]], label=r'$W\left(-\frac{1}{2},\frac{1}{2}\right)$')
#pl.ylabel(r'$\langle\psi_{\uparrow,\downarrow}(t)\vert\psi_{\uparrow,\downarrow}(t)\rangle$')
#pl.xlabel('time (ms)')
#pl.axis( [0, 1, 0, 1] )
#pl.grid()
#ax = pl.gca()
#ax.set_aspect(5)
#pl.legend(loc='best')
#pl.savefig("/Users/miMac/Documents/versionControlledFiles/miThesis/gfx/Ehrenfest/labframeNoFlip.eps")

pl.figure(1)
pl.plot(t*1e3, pup, label=r'$\langle\phi_{\uparrow}\vert\phi_{\uparrow}\rangle$')
pl.plot(t*1e3, pdn, label=r'$\langle\phi_{\downarrow}\vert\phi_{\downarrow}\rangle$')
pl.plot(t*1e3, pup+pdn, label=r'$\langle\Phi\vert\Phi\rangle$')
pl.plot([t[0]*1e3, t[-1]*1e3], [majProbRot, majProbRot], label=r'$W\left(-\frac{1}{2},\frac{1}{2}\right)$')
pl.ylabel(r'$\langle\phi_{\uparrow,\downarrow}(t)\vert\phi_{\uparrow,\downarrow}(t)\rangle$')
pl.xlabel('time (ms)')
#pl.axis( [0, 1, 0, 1] )
pl.grid()
ax = pl.gca()
#ax.set_aspect(5)
pl.legend(loc='best')
#pl.savefig("/Users/miMac/Documents/versionControlledFiles/miThesis/gfx/Ehrenfest/ehrenfestSpin.eps")
np.savez('localPops.npz',t=t,pup=pup,pdn=pdn)

pl.figure(2)
lns1 = pl.plot(t*1e3, x*1.e6, label=r'$x$')
ax = pl.gca();
ax2 = ax.twinx();
lns2 = ax2.plot(t*1e3, v*1.e3, 'g', label=r'$v$')
ax.set_ylabel( r'$x\, (\mu\mathrm{m})$')
ax2.set_ylabel(r'$v\, (\mathrm{mm}\,\mathrm{s}^{-1})$')
ax2.set_xlabel('time (ms)')
lns = lns1+lns2
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc='best')
ax.grid()
#pl.axis( [-5, 5, 0, 1] )
#ax = pl.gca()
#ax.set_aspect(5)
#pl.savefig("/Users/miMac/Documents/versionControlledFiles/miThesis/gfx/Ehrenfest/ehrenfestPos.eps")
np.savez('trajectory.npz',t=t,x=x,v=v)

pl.figure(3)
pl.plot(t*1e3, Ek/kB*1.e6, label=r'$E_k$')
pl.plot(t*1e3, Ep/kB*1.e6, label=r'$E_p$')
pl.plot(t*1e3, (Ek+Ep)/kB*1.e6, label=r'$E_T$')
pl.ylabel(r'$E\, (\mu\mathrm{K})$')
pl.xlabel('time (ms)')
pl.grid()
#pl.axis( [-5, 5, 0, 1] )
#ax = pl.gca()
#ax.set_aspect(5)
pl.legend(loc='best')
#pl.savefig("/Users/miMac/Documents/versionControlledFiles/miThesis/gfx/Ehrenfest/ehrenfestEnergy.eps")
np.savez('energy.npz',t=t,Ehk=Ek,Ehp=Ep)

pl.figure(1)
pl.plot(t*1e3, pops[0,:], label=r'$\langle\phi_{\uparrow}\vert\phi_{\uparrow}\rangle$')
pl.plot(t*1e3, pops[1,:], label=r'$\langle\phi_{\downarrow}\vert\phi_{\downarrow}\rangle$')
pl.plot(t*1e3, pops[0,:]+pops[1,:], label=r'$\langle\Phi\vert\Phi\rangle$')
pl.plot([t[0]*1e3, t[-1]*1e3], [majProbRot, majProbRot], label=r'$W\left(-\frac{1}{2},\frac{1}{2}\right)$')
pl.ylabel(r'$\langle\phi_{\uparrow,\downarrow}(t)\vert\phi_{\uparrow,\downarrow}(t)\rangle$')
pl.xlabel('time (ms)')
#pl.axis( [0, 1, 0, 1] )
pl.grid()
ax = pl.gca()
#ax.set_aspect(5)
pl.legend(loc='best')
#pl.savefig("/Users/miMac/Documents/versionControlledFiles/miThesis/gfx/Ehrenfest/ehrenfestSpin.eps")

end = time.clock()
print end - start

pl.show()