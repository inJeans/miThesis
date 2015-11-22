""" Short File Description

    Longer file description
"""
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import font_manager as font_manager
from matplotlib import patches as patches
import plotly.plotly as py
import seaborn as sns

def load_npz(filename='majorana_data.npz'):
    """ Read in the wavefunction data from the Schrodinger simulation """

    f = open(filename, 'r')
    npz_file = np.load(f)
    k = npz_file['k']
    nz = npz_file['nz']
    t = npz_file['t']
    pops = npz_file['pops']
    pup = npz_file['prob_up']
    pdn = npz_file['prob_dn']
    f.close()
    return k, nz, t, pops, pup, pdn


def main():
    # Sign in to plotly
    py.sign_in('inJeans', 'hl1vdk1zv2')

    # Set some common figure properties
    path = '/Users/miMac/Library/Fonts/AvenirNext-Regular.ttf'
    prop = font_manager.FontProperties(fname=path, size=10)
    figure_dim = (3.8, 2.25)
    figure_dir = '../../gfx/Chapter-4/'

    # Load wavefunction data
    k, nz, t, pops, pup, pdn = load_npz('no_flip_data.npz')
    tms = t * 1e3  # Time in miliseconds

    # Create plots in the lab frame
    maj_prob = np.array([[np.exp(-0.5*k*np.pi)],
                         [1.-np.exp(-0.5*k*np.pi)]])

    lab_fig = plt.figure(num=1, figsize=figure_dim)
    sns.set_palette('Set2')
    set2_colours = sns.color_palette(palette=None)

    # Plot spin up
    plt.fill_between(tms, pops[0, :]+pops[1, :],
                     linewidth=0.0, facecolor=set2_colours[0])
    # Plot spin down
    plt.fill_between(tms, pops[0, :],
                     linewidth=0.0, facecolor=set2_colours[1])
    # Plot Majorana prediction
    plt.plot([tms[0], tms[-1]], [maj_prob[0], maj_prob[0]],
             color=set2_colours[2], linewidth=0.75)

    plt.ylabel(r'$\langle\psi_{\uparrow,\downarrow}(t)\vert\psi_{\uparrow,\downarrow}(t)\rangle$',
               fontproperties=prop)
    plt.xlabel('time (ms)', fontproperties=prop)

    plt.axis([-2, 2, 0, 1])
    ax = plt.gca()
    ax.set_xticklabels(ax.get_xticks(), fontproperties=prop)
    ax.set_yticklabels(ax.get_yticks(), fontproperties=prop)
    ax.grid(b=True, which='major', color='white',
            linestyle='-', linewidth=0.5)
    ax.set_axisbelow(False)
    lab_fig.savefig(figure_dir+'labframeNoFlip.eps', bbox_inches='tight')

    # Plot lab frame legend
    lab_leg = plt.figure(num=2, figsize={figure_dim[0], 0.15})
    ax = plt.subplot()
    ax.set_axis_off()
    lines = {}
    # Dummy plot for spin-up
    plt.plot([], [], color=set2_colours[0], linewidth=10)
    up_label = r'$\langle\psi_{\uparrow}\vert\psi_{\uparrow}\rangle$'
    # Dummy plot for legend
    plt.plot([], [], color=set2_colours[1], linewidth=10)
    dn_label = r'$\langle\psi_{\downarrow}\vert\psi_{\downarrow}\rangle$'
    # Dummy plot for majorana prediction
    plt.plot([], [], color=set2_colours[2], linewidth=0.75)
    maj_label = r'$W\left(-\frac{1}{2},\frac{1}{2}\right)$'

    lgd = ax.legend([up_label, dn_label, maj_label], ncol=3)
    lab_leg.savefig(figure_dir+'labframe_legend.eps', bbox_inches='tight')

    # Create plots in the co-rotating frame
    maj_prob_rot = 0.5 + nz*(-0.5+maj_prob[0])

    rot_fig = plt.figure(num=3, figsize=figure_dim)
    sns.set_palette('Set3')
    set3_colours = sns.color_palette(palette=None)
    # Plot spin up
    plt.fill_between(tms, pup+pdn, linewidth=0.0, facecolor=set3_colours[4])
    # Plot spin down
    plt.fill_between(tms, pup, linewidth=0.0, facecolor=set3_colours[3])
    # Plot majorana prediction
    plt.plot([tms[0], tms[-1]], [maj_prob_rot, maj_prob_rot],
             color=set2_colours[4], linewidth=0.75)

    plt.ylabel(r'$\langle\phi_{\uparrow,\downarrow}(t)\vert\phi_{\uparrow,\downarrow}(t)\rangle$',
               fontproperties=prop)
    plt.xlabel('time (ms)', fontproperties=prop)

    plt.axis([-2, 2, 0, 1])
    ax = plt.gca()
    ax.set_xticklabels(ax.get_xticks(), fontproperties=prop)
    ax.set_yticklabels(ax.get_yticks(), fontproperties=prop)
    ax.grid(b=True, which='major', color='white',
            linestyle='-', linewidth=0.5)
    ax.set_axisbelow(False)
    # lgd = plt.legend([pup_label, pdn_label, pmaj_label],
    #                  loc='center left', bbox_to_anchor=(1, 0.5))
    # plt.savefig(figure_dir+'rotframeFlip.eps',
    #             bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.savefig(figure_dir+'rotframeNoFlip.eps', bbox_inches='tight')

    # Plot rot frame legend
    lab_leg = plt.figure(num=4, figsize={figure_dim[0], 0.15})
    ax = plt.subplot()
    ax.set_axis_off()
    lines = {}
    # Dummy plot for spin-up
    plt.plot([], [], color=set3_colours[4], linewidth=10)
    pup_label = r'$\langle\phi_{\uparrow}\vert\phi_{\uparrow}\rangle$'
    # Dummy plot for legend
    plt.plot([], [], color=set3_colours[3], linewidth=10)
    pdn_label = r'$\langle\phi_{\downarrow}\vert\phi_{\downarrow}\rangle$'
    # Dummy plot for majorana prediction
    plt.plot([], [], color=set2_colours[4], linewidth=0.75)
    pmaj_label = r'$W\left(-\frac{1}{2},\frac{1}{2}\right)$'

    lgd = ax.legend([pup_label, pdn_label, pmaj_label], ncol=3)
    lab_leg.savefig(figure_dir+'rotframe_legend.eps', bbox_inches='tight')

    plt.show()

if __name__ == "__main__":
    main()
