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
T    = 20.e-6;
B0 = 0.01;
dBdx = 20.;
d2Bdx2 = 40000;

tres = 501;
ntrials = 1e5;
nCells = 100**3 + 1;
dt = 1e-6;

time = np.zeros((tres));
pos = np.zeros((ntrials,3,tres));
vel = np.zeros((ntrials,3,tres));
psiUp = np.zeros((ntrials,2,tres));
psiDn = np.zeros((ntrials,2,tres));
atomID = np.zeros((ntrials,1,tres));
atomIsSpinUp = np.zeros((ntrials,1,tres));

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

dset = f.require_dataset('atomData/atomIsSpinUp',(ntrials,1,tres),False,False);
dset.read_direct(atomIsSpinUp);

f.close()

time = time*1.e3

npzf = np.load('energy.npz')
teh = npzf['t']
Ehk = npzf['Ehk']/kB*1.e6
Ehp = npzf['Ehp']/kB*1.e6

npzf = np.load('localPops.npz')
pupEh = npzf['pup']
pdnEh = npzf['pdn']

npzf = np.load('trajectory.npz')
xEh = npzf['x']*1.e6
vEh = npzf['v']*1.e3

teh = teh*1.e3