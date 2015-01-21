import h5py
import numpy as np
import scipy.stats as stats
import pylab as pl

# These are the "Colour Blind 10" colors as RGB.
tableau20 = [(0, 107, 164), (255, 128, 14), (171, 171, 171), (89, 89, 89),
             (95, 158, 209), (200, 28, 0), (137, 137, 137), (162, 200, 236),
             (255, 188, 121), (207, 207, 207)]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

#tableau20 = [ '#1F77B4', '#AEC7E8', '#FF7F0E', '#FFBB78', '#2CA02C', '#98DF8A' ];

pl.figure(1)
pl.plot(time,Ek, label=r'$\langle E_k\rangle$', color=tableau20[0], lw=2)
pl.plot(time,Ep, label=r'$\langle E_p\rangle$', color=tableau20[1], lw=2)
pl.plot(time,Et, label=r'$\langle E_T\rangle$', color=tableau20[3], lw=2)
pl.plot(teh,Ehk, alpha=0.5, color=tableau20[0])
pl.plot(teh,Ehp, alpha=0.5, color=tableau20[1])
pl.xlabel('time (ms)', fontsize=14)
pl.ylabel(r'$E (\mu\mathrm{K})$', fontsize=14)
pl.axis([0, 3, -5, 5])
pl.legend(loc='best')

# Remove the plot frame lines. They are unnecessary chartjunk.
ax = pl.gca()
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)

# Ensure that the axis ticks only show up on the bottom and left of the plot.
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

# Limit the range of the plot to only where the data is.
# Avoid unnecessary whitespace.
pl.ylim(-3, 5)
pl.xlim(0, 3)

# Make sure your axis ticks are large enough to be easily read.
# You don't want your viewers squinting to read your plot.
pl.yticks(range(-3, 6, 1), [str(x) for x in range(-3, 6, 1)], fontsize=14)
pl.xticks(fontsize=14)

# Remove the tick marks; they are unnecessary with the tick lines we just plotted.
pl.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")

pl.grid()

pl.savefig("/Users/miMac/Documents/versionControlledFiles/miThesis/gfx/MCWF/mcwfMajoranaEnergy.eps")

fig = pl.figure(2)
pl.plot(time, spinUp, label=r'$\langle\phi_{\uparrow}\vert\phi_{\uparrow}\rangle$',
        color=tableau20[0], lw=2)
pl.plot(time, spinDn, label=r'$\langle\phi_{\downarrow}\vert\phi_{\downarrow}\rangle$',
        color=tableau20[1], lw=2)
pl.plot(time, spinUp+spinDn, label=r'$\langle\Phi\vert\Phi\rangle$',
        color=tableau20[3], lw=2)
pl.plot(teh, pupEh, alpha=0.5, color=tableau20[0])
pl.plot(teh, pdnEh, alpha=0.5, color=tableau20[1])
pl.xlabel('time (ms)')
pl.ylabel(r'$\langle\phi_{\uparrow,\downarrow}(t)\vert\phi_{\uparrow,\downarrow}(t)\rangle$')
pl.legend(loc='best')

# Remove the plot frame lines. They are unnecessary chartjunk.
ax = pl.gca()
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)

# Ensure that the axis ticks only show up on the bottom and left of the plot.
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

# Limit the range of the plot to only where the data is.
# Avoid unnecessary whitespace.
pl.ylim(0, 1)
pl.xlim(0, 3)

# Make sure your axis ticks are large enough to be easily read.
# You don't want your viewers squinting to read your plot.
#pl.yticks(range(0, 1, 0.2), [str(x) for x in range(0, 1, 0.2)], fontsize=14)
pl.yticks(fontsize=14)
pl.xticks(fontsize=14)

# Remove the tick marks; they are unnecessary with the tick lines we just plotted.
pl.tick_params(axis="both", which="both", bottom="off", top="off",
               labelbottom="on", left="off", right="off", labelleft="on")

pl.grid()

pl.savefig("/Users/miMac/Documents/versionControlledFiles/miThesis/gfx/MCWF/mcwfMajoranaSpin.eps")

pl.figure(3)
lns1 = pl.plot(time,xmcwf, color=tableau20[0], label=r'$z$', lw=2)
pl.plot(teh,xEh, alpha=0.5, color=tableau20[0])
pl.plot(time,xup,'--', color=tableau20[0])
pl.plot(time,xdn,'--', color=tableau20[0])
ax = pl.gca();
ax2 = ax.twinx();
lns2 = ax2.plot(time,vmcwf, color=tableau20[1], label=r'$v$', lw=2)
ax2.plot(teh,vEh, alpha=0.5, color=tableau20[1])
ax2.plot(time,vup,'--', color=tableau20[1])
ax2.plot(time,vdn,'--', color=tableau20[1])
ax.set_ylabel( r'$x\, (\mu\mathrm{m})$')
ax2.set_ylabel(r'$v\, (\mathrm{mm}\,\mathrm{s}^{-1})$')
ax2.set_xlabel('time (ms)')
lns = lns1+lns2
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc='best')

# Remove the plot frame lines. They are unnecessary chartjunk.
#ax.spines["top"].set_visible(False)
#ax.spines["bottom"].set_visible(False)
#ax.spines["right"].set_visible(False)
#ax.spines["left"].set_visible(False)

# Ensure that the axis ticks only show up on the bottom and left of the plot.
#ax.get_xaxis().tick_bottom()
#ax.get_yaxis().tick_left()

# Limit the range of the plot to only where the data is.
# Avoid unnecessary whitespace.
#pl.ylim(0, 1)
#pl.xlim(0, 3)

# Make sure your axis ticks are large enough to be easily read.
# You don't want your viewers squinting to read your plot.
#pl.yticks(range(0, 1, 0.2), [str(x) for x in range(0, 1, 0.2)], fontsize=14)
#pl.yticks(fontsize=14)
#pl.xticks(fontsize=14)

# Remove the tick marks; they are unnecessary with the tick lines we just plotted.
#pl.tick_params(axis="both", which="both", bottom="off", top="off",
#               labelbottom="on", left="off", right="off", labelleft="on")

ax.grid()


pl.savefig("/Users/miMac/Documents/versionControlledFiles/miThesis/gfx/MCWF/mcwfMajoranaTrajectory.eps")

pl.show()