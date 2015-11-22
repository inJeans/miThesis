""" Short File Description

    Longer file description
"""
from __future__ import print_function

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
    figure_dim = (2.75, 2.)
    figure_dir = '../../gfx/Chapter-4/'

    # Load wavefunction data
    k, nz, t, pops, pup, pdn = load_npz('flip_data.npz')
    tms = t * 1e3  # Time in miliseconds

    # Create plots in the lab frame
    maj_prob = np.array([[np.exp(-0.5*k*np.pi)],
                         [1.-np.exp(-0.5*k*np.pi)]])

    lab_fig, lab_ax = plt.subplots()
    sns.set_palette('Set2')
    set2_colours = sns.color_palette(palette=None)
    # Plot spin up
    lab_ax.plot(tms, pops[0, :]+pops[1, :],
                color=set2_colours[0],
                label=r'$\langle\psi_{\uparrow}\vert\psi_{\uparrow}\rangle$')
    # Plot spin down
    lab_ax.plot(tms, pops[0, :],
                color=set2_colours[1],
                label=r'$\langle\psi_{\downarrow}\vert\psi_{\downarrow}\rangle$')
    # Plot Majorana prediction
    plt.plot([tms[0], tms[-1]], [maj_prob[0], maj_prob[0]],
             color=set2_colours[2],
             label=r'$W\left(-\frac{1}{2},\frac{1}{2}\right)$')

    # lab_ax.ylabel(r'$\langle\psi_{\uparrow,\downarrow}(t)\vert\psi_{\uparrow,\downarrow}(t)\rangle$',
    #                fontproperties=prop)
    # lab_ax.xlabel('time (ms)', fontproperties=prop)

    lab_ax.axis([-2, 2, 0, 1])
    update = {'data': [{'fill': 'tozeroy'}]}
    print("Uploading data...")
    plot_url = py.plot_mpl(lab_fig, update=update, filename='lab_frame_flip')

    # # Create plots in the co-rotating frame
    # maj_prob_rot = 0.5 + nz*(-0.5+maj_prob[0])

    # rot_fig = plt.figure(num=2, figsize=figure_dim)
    # sns.set_palette('Set3')
    # set3_colours = sns.color_palette(palette=None)
    # # Plot spin up
    # plt.fill_between(tms, pup+pdn, linewidth=0.0, facecolor=set3_colours[4])
    # # Dummy plot for legend
    # plt.plot([], [], color=set3_colours[4], linewidth=10)
    # pup_label = r'$\langle\phi_{\uparrow}\vert\phi_{\uparrow}\rangle$'
    # # Plot spin down
    # plt.fill_between(tms, pup, linewidth=0.0, facecolor=set3_colours[1])
    # # Dummy plot for legend
    # plt.plot([], [], color=set3_colours[1], linewidth=10)
    # pdn_label = r'$\langle\phi_{\downarrow}\vert\phi_{\downarrow}\rangle$'
    # # Plot majorana prediction
    # plt.plot([tms[0], tms[-1]], [maj_prob_rot, maj_prob_rot],
    #          color=set3_colours[2])
    # pmaj_label = r'$W\left(-\frac{1}{2},\frac{1}{2}\right)$'

    # plt.ylabel(r'$\langle\phi_{\uparrow,\downarrow}(t)\vert\phi_{\uparrow,\downarrow}(t)\rangle$',
    #            fontproperties=prop)
    # plt.xlabel('time (ms)', fontproperties=prop)

    # plt.axis([-2, 2, 0, 1])
    # ax = plt.gca()
    # ax.set_xticklabels(ax.get_xticks(), fontproperties=prop)
    # ax.set_yticklabels(ax.get_yticks(), fontproperties=prop)
    # ax.grid(b=True, which='major', color='white',
    #         linestyle='-', linewidth=0.5)
    # plot_url = py.plot_mpl(rot_fig)
    # ax.set_axisbelow(False)
    # lgd = plt.legend([pup_label, pdn_label, pmaj_label],
    #                  loc='center left', bbox_to_anchor=(1, 0.5))
    # plt.savefig(figure_dir+'rotframeFlip.eps',
    #             bbox_extra_artists=(lgd,), bbox_inches='tight')


    # plt.show()

if __name__ == "__main__":
    main()
