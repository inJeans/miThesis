import numpy as np
import scipy.stats as stats
import pylab as pl
import os
import time
import seaborn as sns

# sns.set(style='ticks', palette='Set2')
sns.set_palette('Set2')
set2_colours = sns.color_palette(palette=None)
# sns.set(color_codes=True)

start = time.clock()

gs = 0.5;
muB = 9.27401e-24;
hbar = 1.05457e-34;

dt = 1.e-6;
nsteps = np.ceil(0.015 / dt * 2).astype(int);

t    = np.zeros( (nsteps) );
t[0] = -0.015;

Bx = 2.0e-7;
C = 0.0025;

Bz =-C*t[0];
magB = np.sqrt(Bx**2 + Bz**2);
nx = Bx / magB;
nz = Bz / magB;

psi      = np.zeros( (2,nsteps), dtype=complex );
psi[:,0] = 0.5 * np.array( [ [1 + nx + nz],
                             [1 + nx - nz] ] )[:,0] / np.sqrt(1+nx);

Pup = 0.5*np.array( [ [1+nz, nx], [ nx, 1-nz ] ] );
Pdn = 0.5*np.array( [ [1-nz,-nx], [-nx, 1+nz ] ] );

pup    = np.zeros( (nsteps) );
pup[0] = np.dot( psi[:,0].conj().transpose(), np.dot( Pup, psi[:,0] ) ).real;
pdn    = np.zeros( (nsteps) );
pdn[0] = np.dot( psi[:,0].conj().transpose(), np.dot( Pdn, psi[:,0] ) ).real;

for i in range(1,nsteps):
    
    t[i] = t[i-1] + dt;
    
    Bz =-C*t[i];
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
                                     
pops = (psi.conj() * psi).real;
                                     
k = gs*muB*Bx**2 / (hbar*C);
                                     
print "\nProbability of spin flip is %f\n" % np.exp(-0.5*k*np.pi);
print "This simulation predicts %f\n" % pops[0,-1];

majProb = np.array( [ [np.exp(-0.5*k*np.pi)], [1.-np.exp(-0.5*k*np.pi)] ] )
majProbRot = 0.5 + nz*(-0.5+majProb[0]);

pl.figure(1)
pl.plot(t*1e3, pops[0,:], label=r'$\langle\psi_{\uparrow}\vert\psi_{\uparrow}\rangle$')
pl.plot(t*1e3, pops[1,:], label=r'$\langle\psi_{\downarrow}\vert\psi_{\downarrow}\rangle$')
pl.plot(t*1e3, pops[0,:]+pops[1,:], label=r'$\langle\Psi\vert\Psi\rangle$')
pl.plot([t[0]*1e3, t[-1]*1e3], [majProb[0], majProb[0]], label=r'$W\left(-\frac{1}{2},\frac{1}{2}\right)$')
pl.ylabel(r'$\langle\psi_{\uparrow,\downarrow}(t)\vert\psi_{\uparrow,\downarrow}(t)\rangle$')
pl.xlabel('time (ms)')
pl.axis( [-5, 5, 0, 1] )
ax = pl.gca()
ax.set_aspect(5)
pl.legend(loc='best')
pl.savefig("/Users/miMac/Documents/versionControlledFiles/miThesis/gfx/Ehrenfest/labframeNoFlip.eps")
                                     
pl.figure(2)
pl.plot(t*1e3, pup, label=r'$\langle\phi_{\uparrow}\vert\phi_{\uparrow}\rangle$')
pl.plot(t*1e3, pdn, label=r'$\langle\phi_{\downarrow}\vert\phi_{\downarrow}\rangle$')
pl.plot(t*1e3, pup+pdn, label=r'$\langle\Phi\vert\Phi\rangle$')
pl.plot([t[0]*1e3, t[-1]*1e3], [majProbRot, majProbRot], label=r'$W\left(-\frac{1}{2},\frac{1}{2}\right)$')
pl.ylabel(r'$\langle\phi_{\uparrow,\downarrow}(t)\vert\phi_{\uparrow,\downarrow}(t)\rangle$')
pl.xlabel('time (ms)')
pl.axis( [-5, 5, 0, 1] )
ax = pl.gca()
ax.set_aspect(5)
pl.legend(loc='best')
pl.savefig("/Users/miMac/Documents/versionControlledFiles/miThesis/gfx/Ehrenfest/rotframeNoFlip.eps")

pl.figure(3)
# plops = np.row_stack((pops[0,:],pops[1,:]))
# fig = pl.stackplot( t*1e3, plops )
# sns.tsplot(pops[0,:])
pl.fill_between(t*1e3,0,pops[1,:]+pops[0,:],linewidth=0.0,facecolor=set2_colours[0])
pl.fill_between(t*1e3,0,pops[0,:],linewidth=0.0,facecolor=set2_colours[1])
pl.ylabel(r'$\langle\psi_{\uparrow,\downarrow}(t)\vert\psi_{\uparrow,\downarrow}(t)\rangle$')
pl.xlabel('time (ms)')
pl.axis( [-5, 5, 0, 1] )
# pl.setp(fig, edgecolor=set2_colours[0])
# pl.setp(fig, linewidth=0.0)
ax = pl.gca()
ax.set_aspect(5)
pl.savefig("/Users/miMac/Desktop/test.eps")
pl.show()

end = time.clock()
print end - start

pl.show()