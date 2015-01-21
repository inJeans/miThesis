import h5py
import numpy as np
import os

# Set some physical constants
gs   = -0.5 * -1.0;
muB  = 9.27400915e-24;
mRb  = 1.443160648e-25;
pi   = 3.14159265;
kB   = 1.3806503e-23;
hbar = 1.05457148e-34;
T    = 20.e-6;
dBdr = 8746.;

tres = 101;
ntrials = 1e6;
dt = 1e-6;

ts  = ['7', '8', '9', '11', '12', '13'];
ran = [0.15, 0.13, 0.3, 0.26, 0.22, 0.2];

#for t in range(0,6):
for t in [2]:

    time = np.zeros((tres));
    vel = np.zeros((ntrials,3,tres));
    isPerturb = np.zeros((ntrials,1,tres));
    N = np.zeros((tres));
    collisionCount = np.zeros((10**3+1,1,tres));

    f = h5py.File('/Users/miMac/Documents/versionControlledFiles/miThesis/Code/HomoWalraven/' + ts[t] + 'uK.h5')

    dset = f.require_dataset('atomData/simuatedTime',(1,1,tres),False,False);
    dset.read_direct(time);

    dset = f.require_dataset('atomData/velocities',(ntrials,3,tres),False,False);
    dset.read_direct(vel);

    dset = f.require_dataset('atomData/isPerturb',(ntrials,1,tres),False,False);
    dset.read_direct(isPerturb);

    dset = f.require_dataset('atomData/atomNumber',(1,1,tres),False,False);
    dset.read_direct(N);

    dset = f.require_dataset('atomData/collisionCount',(10**3+1,1,tres),False,False);
    dset.read_direct(collisionCount);

    f.close()

    totalColl = np.sum( collisionCount, 0 )[0,:];
    collRate = totalColl / ( time[1]-time[0] ) / N ;

    time = time * collRate;

    Temp = np.zeros((N.size,))
    Tperturb = np.zeros((N.size,))

    for i in range(0,N.size):
        kinetic = 0.5 * mRb * np.sum(vel[0:N[i],:,i]**2, 1)
        n = np.where( np.isfinite(kinetic) )
        Temp[i] = 2./3. * np.sum( kinetic[n], 0) / N[i] / kB * 1.e6
        kineticPerturb = isPerturb[0:N[i],0,i] * 0.5 * mRb * np.sum(vel[0:N[i],:,i]**2, 1)
        Tperturb[i] = 2./3. * np.sum( kineticPerturb[n], 0) / (0.01*ntrials) / kB * 1.e6

    if t== 2 or t == 3 or t == 5:
        vel = np.zeros((ntrials,3,tres));
        isPerturb = np.zeros((ntrials,1,tres));
        N = np.zeros((tres));
        collisionCount = np.zeros((10**3+1,1,tres));

        f = h5py.File('/Users/miMac/Documents/versionControlledFiles/miThesis/Code/HomoWalraven/' + ts[t] + 'buK.h5')
    
        dset = f.require_dataset('atomData/velocities',(ntrials,3,tres),False,False);
        dset.read_direct(vel);
    
        dset = f.require_dataset('atomData/isPerturb',(ntrials,1,tres),False,False);
        dset.read_direct(isPerturb);
    
        dset = f.require_dataset('atomData/atomNumber',(1,1,tres),False,False);
        dset.read_direct(N);
    
        f.close()
    
        Temp = np.zeros((N.size,))
        Tperturb2 = np.zeros((N.size,))
    
        for i in range(0,N.size):
            kinetic = 0.5 * mRb * np.sum(vel[0:N[i],:,i]**2, 1)
            n = np.where( np.isfinite(kinetic) )
            Temp[i] = 2./3. * np.sum( kinetic[n], 0) / N[i] / kB * 1.e6
            kineticPerturb = isPerturb[0:N[i],0,i] * 0.5 * mRb * np.sum(vel[0:N[i],:,i]**2, 1)
            Tperturb2[i] = 2./3. * np.sum( kineticPerturb[n], 0) / (0.01*ntrials) / kB * 1.e6

        Tperturb = 0.5 * (Tperturb + Tperturb2);

    z, cov = np.polyfit( time[0:ran[t]*tres], np.log(np.abs(Temp[-1] - Tperturb[0:ran[t]*tres])), 1, cov=True )
    
    print "The temperature is " + ts[t]
    print "The thermalisation time is", -z[0]
    print "Thermalisation in %f collisions", (1./-z[0])*4

    np.save('Te' + ts[t] + 'uK.npy', Temp)
    np.save('Tp' + ts[t] + 'uK.npy', Tperturb)
    np.save('ti' + ts[t] + 'uK.npy', time)

