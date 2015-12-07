""" Short File Description

    Longer file description
"""
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import font_manager as font_manager
import argparse
# import plotly.plotly as py
import seaborn as sns


def main(infile='majorana_data',
         outfile='majorana',
         outdir='../../gfx'):
    # Sign in to plotly
    # py.sign_in('inJeans', 'hl1vdk1zv2')

    # Set some common figure properties
    # path = '/Users/miMac/Library/Fonts/AvenirNext-Regular.ttf'
    path = '/Users/miMacbookPro/Library/Fonts/avenir-next-regular.ttf'
    prop = font_manager.FontProperties(fname=path, size=10)
    figure_dim = (3.8, 2.25)

    print("Loading data file.")
    # Load wavefunction data
    k, nz, t, lab_prob, rot_prob = load_npz(infile+'.npz')
    tms = t * 1e3  # Time in miliseconds

    print("Creating figure in the lab frame.")
    # Create plots in the lab frame
    maj_prob = np.array([[np.exp(-0.5*k*np.pi)],
                         [1.-np.exp(-0.5*k*np.pi)]])

    lab_fig = plt.figure(num=1, figsize=figure_dim)
    sns.set_palette('Set2')
    set2_colours = sns.color_palette(palette=None)

    # Plot spin up
    plt.fill_between(tms, lab_prob[0, :]+lab_prob[1, :],
                     linewidth=0.0, facecolor=set2_colours[0])
    # Plot spin down
    plt.fill_between(tms, lab_prob[0, :],
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
    print("Saving lab frame figure as {0}"
          .format(outdir+outfile+'_lab_frame.eps'))
    lab_fig.savefig(outdir+outfile+'_lab_frame.eps', bbox_inches='tight')

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
    print("Saving lab frame legend as {0}"
          .format(outdir+'lab_frame_legend.eps'))
    lab_leg.savefig(outdir+'lab_frame_legend.eps', bbox_inches='tight')

    print("Creating plots in the rotating frame.")
    # Create plots in the co-rotating frame
    maj_prob_rot = 0.5 + nz*(-0.5+maj_prob[0])

    rot_fig = plt.figure(num=3, figsize=figure_dim)
    sns.set_palette('Set3')
    set3_colours = sns.color_palette(palette=None)
    # Plot spin up
    plt.fill_between(tms, rot_prob[0, :]+rot_prob[1, :],
                     linewidth=0.0, facecolor=set3_colours[4])
    # Plot spin down
    plt.fill_between(tms, rot_prob[0, :],
                     linewidth=0.0, facecolor=set3_colours[3])
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
    print("Saving rotating frame figure as {0}"
          .format(outdir+outfile+'_rot_frame.eps'))
    plt.savefig(outdir+outfile+'_rot_frame.eps', bbox_inches='tight')

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
    print("Saving rotating frame legend as {0}"
          .format(outdir+'rot_frame_legend.eps'))
    lab_leg.savefig(outdir+'rot_frame_legend.eps', bbox_inches='tight')

    plt.show()


def load_npz(filename='majorana_data.npz'):
    """ Read in the wavefunction data from the Schrodinger simulation """

    f = open(filename, 'r')
    npz_file = np.load(f)
    k = npz_file['k']
    nz = npz_file['nz']
    t = npz_file['t']
    lab_prob = npz_file['lab_prob']
    rot_prob = npz_file['rot_prob']
    f.close()
    return k, nz, t, lab_prob, rot_prob


def data_crunch_cli():
    parser = argparse.ArgumentParser(description="""A differential equation
                                     solver for the archetypeal Majorana spin
                                     flip problem.""")
    parser.add_argument('--input', '-i', dest='infile',
                        default='majorana_data',
                        help="""Name of the input npz file.""")
    parser.add_argument('--output', '-o', dest='outfile',
                        default='majorana',
                        help="""Name of the output file.""")
    parser.add_argument('--outdir', '-d', dest='outdir',
                        default='../../gfx',
                        help="""Name of the output file directory.""")
    args = parser.parse_args()

    print("*********************************************")
    print("* RUNNING MAJORANA PROBLEM DATA CRUNCHER    *")
    print("*********************************************")
    print(" ")
    print("File names are:")
    print("input filename - {0}.npz".format(args.infile))
    print("output filename prefix - {0}".format(args.outfile))
    print("output file directory - {0}".format(args.outdir))
    print(" ")
    print("*********************************************")
    print(" ")

    main(infile=args.infile,
         outfile=args.outfile,
         outdir=args.outdir)

if __name__ == "__main__":
    data_crunch_cli()
