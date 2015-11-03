from __future__ import print_function

import numpy as np
import scipy.stats as stats
import os
import time

def main():
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
                                     
    print( "Probability of spin flip is {0}".format( np.exp(-0.5*k*np.pi) ) );
    print( "This simulation predicts    {0}\n".format( pops[0,-1] ) );

    end = time.clock()
    print( "This simulation took {0}s.".format( end - start ) );

    f = open('majorana_data.npz', 'w')
    np.savez( f, k=k, nz=nz, t=t, pops=pops, pup=pup, pdn=pdn )
    f.close() 

if __name__ == "__main__":
    main()