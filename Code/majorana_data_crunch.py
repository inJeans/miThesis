import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

def load_npz( filename="majorana_data.npz" ):
	f = open(filename, 'r')
	npz_file = np.load(f)
	k = npz_file['k']
	nz = npz_file['nz']
	t = npz_file['t']
	pops = npz_file['pops']
	pup = npz_file['pup']
	pdn = npz_file['pdn']
	f.close()
	return k, nz, t, pops, pup, pdn

def main():
	k, nz, t, pops, pup, pdn = load_npz()
	
	majProb = np.array( [ [np.exp(-0.5*k*np.pi)], 
		                  [1.-np.exp(-0.5*k*np.pi)] ] )
	majProbRot = 0.5 + nz*(-0.5+majProb[0]);

	plt.figure(1)
	sns.set_palette('Set2')
	set2_colours = sns.color_palette(palette=None)
	plt.fill_between(t*1e3, pops[0,:]+pops[1,:], linewidth=0.0, facecolor=set2_colours[0], 
		             label=r'$\langle\psi_{\uparrow}\vert\psi_{\uparrow}\rangle$')
	plt.fill_between(t*1e3, pops[0,:], linewidth=0.0, facecolor=set2_colours[1], zorder=1,
		             label=r'$\langle\psi_{\downarrow}\vert\psi_{\downarrow}\rangle$')
	plt.plot([t[0]*1e3, t[-1]*1e3], [majProb[0], majProb[0]], color=set2_colours[2], 
		     label=r'$W\left(-\frac{1}{2},\frac{1}{2}\right)$')
	plt.ylabel(r'$\langle\psi_{\uparrow,\downarrow}(t)\vert\psi_{\uparrow,\downarrow}(t)\rangle$')
	plt.xlabel('time (ms)')
	plt.axis( [-2, 2, 0, 1] )
	ax = plt.gca()
	ax.set_aspect(3)
	ax.grid(b=True,which='major',color='white',alpha=0.8,linestyle='-')
	ax.set_axisbelow(False)
	plt.legend(loc='best')
	plt.savefig("/Users/miMac/Documents/versionControlledFiles/miThesis/gfx/Ehrenfest/labframeNoFlip.eps")
                                     
	plt.figure(2)
	sns.set_palette('Set3')
	set3_colours = sns.color_palette(palette=None)
	plt.fill_between(t*1e3, pup+pdn, linewidth=0.0, facecolor=set3_colours[4], 
		             label=r'$\langle\phi_{\uparrow}\vert\phi_{\uparrow}\rangle$')
	plt.fill_between(t*1e3, pup, linewidth=0.0, facecolor=set3_colours[1], 
		             label=r'$\langle\phi_{\downarrow}\vert\phi_{\downarrow}\rangle$')
	plt.plot([t[0]*1e3, t[-1]*1e3], [majProbRot, majProbRot], color=set3_colours[2], 
		     label=r'$W\left(-\frac{1}{2},\frac{1}{2}\right)$')
	plt.ylabel(r'$\langle\phi_{\uparrow,\downarrow}(t)\vert\phi_{\uparrow,\downarrow}(t)\rangle$')
	plt.xlabel('time (ms)')
	plt.axis( [-2, 2, 0, 1] )
	ax = plt.gca()
	ax.set_aspect(3)
	plt.legend(loc='best')
	plt.savefig("/Users/miMac/Documents/versionControlledFiles/miThesis/gfx/Ehrenfest/rotframeNoFlip.eps")

	plt.show()

if __name__ == "__main__":
	main()