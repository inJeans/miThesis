import h5py
import numpy as np
import scipy.stats as stats
import pylab as pl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

# These are the "Color Blind 10" colors as RGB.
tableau20 = [(0, 107, 164), (255, 128, 14), (171, 171, 171), (89, 89, 89),
             (95, 158, 209), (200, 28, 0), (137, 137, 137), (162, 200, 236),
             (255, 188, 121), (207, 207, 207)]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

pl.figure(1)
pl.subplot(3,1,1)
pl.plot(time,(Et-Et[0])/Et[0]*100.,
        lw=2.5, color=tableau20[3])
pl.ylabel(r'$\Delta E_T\,(\%)$', fontsize=14)
ax = pl.gca()
ax.locator_params('y',tight=True, nbins=4)
ax.locator_params('x',tight=True, nbins=6)
pl.yticks(fontsize=14)
# make these tick labels invisible
pl.setp( ax.get_xticklabels(), visible=False)
ax.grid()

pl.subplot(3,1,2)
pl.plot(time,dNorm*100.,
        lw=2.5, color=tableau20[4]);
pl.ylabel(r'$\Delta \vert\Psi\vert^2\,(\%)$', fontsize=14)
ax = pl.gca()
ax.locator_params('y',tight=True, nbins=4)
ax.locator_params('x',tight=True, nbins=6)
pl.yticks(fontsize=14)
# make these tick labels invisible
pl.setp( ax.get_xticklabels(), visible=False)
ax.grid()

pl.subplot(3,1,3)
pl.plot(time,N*Fn,
        lw=2.5, color=tableau20[5]);
pl.ylabel(r'$N$', fontsize=14)
ax = pl.gca()
ax.locator_params('y',tight=True, nbins=3)
ax.locator_params('x',tight=True, nbins=6)
ax.set_yticklabels([' ', r'$10^6$', ''])
pl.yticks(fontsize=14)
pl.xticks(fontsize=14)
pl.xlabel(r'$t \times \tau_{c}$', fontsize=14)
ax.grid()

# Remove the plot frame lines. They are unnecessary chartjunk.
#ax = pl.gca()
#ax.spines["top"].set_visible(False)
#ax.spines["bottom"].set_visible(False)
#ax.spines["right"].set_visible(False)
#ax.spines["left"].set_visible(False)

# Ensure that the axis ticks only show up on the bottom and left of the plot.
#ax.get_xaxis().tick_bottom()
#ax.get_yaxis().tick_left()

# Limit the range of the plot to only where the data is.
# Avoid unnecessary whitespace.
#pl.ylim(-3, 5)
#pl.xlim(0, 3)

# Make sure your axis ticks are large enough to be easily read.
# You don't want your viewers squinting to read your plot.
#pl.yticks(range(-3, 6, 1), [str(x) for x in range(-3, 6, 1)], fontsize=14)
#pl.xticks(fontsize=14)

# Remove the tick marks; they are unnecessary with the tick lines we just plotted.
#pl.tick_params(axis="both", which="both", bottom="off", top="off",
#                labelbottom="on", left="off", right="off", labelleft="on")

#pl.grid()

pl.savefig("/Users/miMac/Documents/versionControlledFiles/miThesis/gfx/MCWF/mcwfIPConserve.eps")

fig = pl.figure(2)
pl.plot( binsri, nri*Fn, label=r'$n(\mathbf{r},t_\mathrm{inital})$',
        lw=2.5, color=tableau20[0] )
pl.plot( binsrf, nrf*Fn, label=r'$n(\mathbf{r},t_\mathrm{final})$',
        lw=2.5, color=tableau20[1] )
pl.plot( binsri, maxwell*Fn, label=r'$n_\mathrm{IP}(\mathbf{r})$',
        lw=2.5, color=tableau20[2] )
pl.ylabel(r'$n(\mathbf{r})$', fontsize=14)
pl.xlabel(r'$\vert\mathbf{r}\vert$', fontsize=14)
ax = pl.gca()
ax.locator_params('x',tight=True, nbins=6)
ax.locator_params('y',tight=True, nbins=6)
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
#pl.ylim(0, 1)
#pl.xlim(0, 3)

# Make sure your axis ticks are large enough to be easily read.
# You don't want your viewers squinting to read your plot.
#pl.yticks(range(0, 1, 0.2), [str(x) for x in range(0, 1, 0.2)], fontsize=14)
pl.yticks(fontsize=14)
pl.xticks(fontsize=14)

# Remove the tick marks; they are unnecessary with the tick lines we just plotted.
pl.tick_params(axis="both", which="both", bottom="off", top="off",
               labelbottom="on", left="off", right="off", labelleft="on")

pl.grid()

pl.savefig("/Users/miMac/Documents/versionControlledFiles/miThesis/gfx/MCWF/mcwfIPDist.eps")

pl.figure(3)

pl.plot(time,Ek,time,Ep)

pl.show()