import h5py
import numpy as np
import pylab as pl
import os

# These are the "Colour Blind 10" colors as RGB.
tableau20 = [(0, 107, 164), (255, 128, 14), (171, 171, 171), (89, 89, 89),
             (95, 158, 209), (200, 82, 0), (162, 200, 236), (137, 137, 137),
             (255, 188, 121), (207, 207, 207)]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

tres = 101;
ntrials = 1e6;
dt = 1e-6;

ts  = ['7', '8', '9', '11', '12', '13'];
ran = [0.15, 0.13, 0.15, 0.28, 0.22, 0.06];

markers = ['^', '<', 's', 'd', '>', 'v']

pl.figure(1)
pl.xlabel(r'$t / \tau_\mathrm{c}$', fontsize=14)
pl.ylabel(r'$T/T_\mathrm{bulk}(0)$', fontsize=14)

plotIndices = np.linspace(0,tres-1,20).astype(int);

time = np.load('ti11uK.npy')*4;
temp = np.load('Te11uK.npy');
pl.plot( time[plotIndices], temp[plotIndices]/temp[0], marker='o', markersize=9, markeredgecolor='none', color=tableau20[6], lw=0)

for t in range(0,6):
    time = np.load('ti'  + ts[t] + 'uK.npy')*4;
    temp = np.load('Te'  + ts[t] + 'uK.npy');
    tper = np.load('Tp' + ts[t] + 'uK.npy');

    z, cov = np.polyfit( time[0:ran[t]*tres], np.log(np.abs(temp[-1] - tper[0:ran[t]*tres])), 1, cov=True )

    print "The thermalisation time is", -z[0]
    print "Thermalisation in %f collisions", (1./-z[0])
    print cov[1,1]

    name = r'$T_{p} = $' + ts[t] + r'$\mu\mathrm{K},\, \tau/\tau_\mathrm{c} =$' + '%.2f' % (1./-z[0]) + r'$\pm$' + '%.3f' % (cov[1,1])

    pl.plot( time, (temp[-1] + np.sign(tper[0] - temp[0]) * np.exp(z[1] + z[0]*time))/temp[0], color=tableau20[t], lw=1 )
    pl.plot( time[plotIndices], tper[plotIndices]/temp[0], marker=markers[t], markersize=9, markeredgecolor='none', color=tableau20[t], lw=0,label=name )

pl.ylim(0.65, 1.4)
pl.xlim(0, 15)
pl.legend(loc='best',fontsize=12)

pl.savefig("/Users/miMac/Documents/versionControlledFiles/miThesis/gfx/Thermalisation/walravenHomo.eps")

pl.show()
