""" Short File Description

    Longer file description
"""
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import font_manager as font_manager
import argparse
# import plotly.plotly as py
import seaborn as sns

from vector import Vector

# Physical constants
G_S = 0.5
MU_B = 9.27401e-24
H_BAR = 1.05457e-34
K_B = 1.3806503e-23

M_RB87 = 1.443160648e-25

def main(infile='ehrenfest_data',
         outfile='ehrenfest',
         outdir='../../gfx/'):
    # Sign in to plotly
    # py.sign_in('inJeans', 'hl1vdk1zv2')

    # Set some common figure properties
    path = '/Users/miMac/Library/Fonts/AvenirNext-Regular.ttf'
    # path = '/Users/miMacbookPro/Library/Fonts/avenir-next-regular.ttf'
    prop = font_manager.FontProperties(fname=path, size=10)
    figure_dim = (3.8, 2.25)

    print("Loading data file.")
    # Load wavefunction data
    k, t, B, pos, vel, lab_prob, rot_prob, potential = load_npz(infile+'.npz')
    tms = t * 1e3  # Time in miliseconds

    print("Creating plots in the rotating frame.")
    # Create plots in the co-rotating frame
    n = B[-1].unit()
    maj_prob = np.array([[np.exp(-0.5*k*np.pi)],
                         [1.-np.exp(-0.5*k*np.pi)]])
    maj_prob_rot = 0.5 + n.z*(-0.5+maj_prob[0])

    rot_fig = plt.figure(num=1, figsize=figure_dim)
    sns.set_palette('Set3')
    spin_colours = sns.color_palette(palette=None)
    sns.set_palette('Set2')
    maj_colours = sns.color_palette(palette=None)
    # Plot spin up
    plt.fill_between(tms, rot_prob[0, :]+rot_prob[1, :],
                     linewidth=0.0, facecolor=spin_colours[4])
    # Plot spin down
    plt.fill_between(tms, rot_prob[0, :],
                     linewidth=0.0, facecolor=spin_colours[3])
    # Plot majorana prediction
    plt.plot([tms[0], tms[-1]], [maj_prob_rot, maj_prob_rot],
             color=maj_colours[4], linewidth=0.75)

    plt.ylabel(r'$\langle\phi_{\uparrow,\downarrow}(t)\vert\phi_{\uparrow,\downarrow}(t)\rangle$',
               fontproperties=prop)
    plt.xlabel('time (ms)', fontproperties=prop)

    plt.axis([0, 3, 0, 1])
    ax = plt.gca()
    ax.set_xticklabels(ax.get_xticks(), fontproperties=prop)
    ax.set_yticklabels(ax.get_yticks(), fontproperties=prop)
    ax.grid(b=True, which='major', color='white',
            linestyle='-', linewidth=0.5)
    ax.set_axisbelow(False)
    print("Saving rotating frame figure as {0}"
          .format(outdir+outfile+'_spin.eps'))
    plt.savefig(outdir+outfile+'_spin.eps', bbox_inches='tight')

    # Plot rot frame legend
    lab_leg = plt.figure(num=2, figsize={figure_dim[0], 0.15})
    ax = plt.subplot()
    ax.set_axis_off()
    lines = {}
    # Dummy plot for spin-up
    plt.plot([], [], color=spin_colours[4], linewidth=10)
    pup_label = r'$\langle\phi_{\uparrow}\vert\phi_{\uparrow}\rangle$'
    # Dummy plot for legend
    plt.plot([], [], color=spin_colours[3], linewidth=10)
    pdn_label = r'$\langle\phi_{\downarrow}\vert\phi_{\downarrow}\rangle$'
    # Dummy plot for majorana prediction
    plt.plot([], [], color=maj_colours[4], linewidth=0.75)
    pmaj_label = r'$W\left(-\frac{1}{2},\frac{1}{2}\right)$'

    lgd = ax.legend([pup_label, pdn_label, pmaj_label], ncol=3)
    print("Saving rotating frame legend as {0}"
          .format(outdir+'spin_legend.eps'))
    lab_leg.savefig(outdir+'spin_legend.eps', bbox_inches='tight')

    # Create energy plots
    print("Creating energy plots.")
    kinetic = 0.5 * M_RB87 * np.array([v.norm()**2 for v in vel]) / K_B * 1e6
    potential = potential / K_B * 1e6
    max_potential = max(potential) # shift potential zero
    potential -= max_potential

    energy_fig = plt.figure(num=3, figsize=figure_dim)
    sns.set_palette('Set1')
    energy_colours = sns.color_palette(palette=None)
    # Plot potential energy
    plt.fill_between(tms, kinetic-potential,
                     linewidth=0.0, facecolor=energy_colours[4])
    # Plot kinetic energy
    plt.fill_between(tms, kinetic,
                     linewidth=0.0, facecolor=energy_colours[3])

    plt.ylabel(r'Energy ($\mu$K)',
               fontproperties=prop)
    plt.xlabel('time (ms)', fontproperties=prop)

    plt.axis([0, 3, 0, max(kinetic-potential)])
    ax = plt.gca()
    ax.set_xticklabels(ax.get_xticks(), fontproperties=prop)
    ax.set_yticklabels(ax.get_yticks(), fontproperties=prop)
    ax.grid(b=True, which='major', color='white',
            linestyle='-', linewidth=0.5)
    ax.set_axisbelow(False)
    print("Saving energy figure as {0}"
          .format(outdir+outfile+'_energy.eps'))
    plt.savefig(outdir+outfile+'_energy.eps', bbox_inches='tight')

    # Plot energy legend
    lab_leg = plt.figure(num=4, figsize={figure_dim[0], 0.15})
    ax = plt.subplot()
    ax.set_axis_off()
    lines = {}
    # Dummy plot for spin-up
    plt.plot([], [], color=energy_colours[4], linewidth=10)
    potential_label = r'potential ($\mu$K)'
    # Dummy plot for legend
    plt.plot([], [], color=energy_colours[3], linewidth=10)
    kinetic_label = r'kinetic ($\mu$K)'

    lgd = ax.legend([potential_label, kinetic_label], ncol=2)
    print("Saving energy legend as {0}"
          .format(outdir+'energy_legend.eps'))
    lab_leg.savefig(outdir+'energy_legend.eps', bbox_inches='tight')

    # Plot dynamics
    print("Creating dynamics plots.")
    dynamics_fig = plt.figure(num=5, figsize=figure_dim)
    z = [p.z*1e3 for p in pos]
    vz = [v.z for v in vel]
    plt.plot(tms, z, tms, vz)
    plt.ylabel(r'Position (m)',
               fontproperties=prop)
    plt.xlabel('time (ms)', fontproperties=prop)

    plt.axis([0, 3, 0, max(z)])
    ax = plt.gca()
    ax.set_xticklabels(ax.get_xticks(), fontproperties=prop)
    ax.set_yticklabels(ax.get_yticks(), fontproperties=prop)
    ax.grid(b=True, which='major', color='white',
            linestyle='-', linewidth=0.5)
    ax.set_axisbelow(False)
    print("Saving dynamics figure as {0}"
          .format(outdir+outfile+'_dynamics.eps'))
    plt.savefig(outdir+outfile+'_dynamics.eps', bbox_inches='tight')
    plt.show()


def load_npz(filename='ehrenfest_data.npz'):
    """ Read in the wavefunction data from the Schrodinger simulation """

    f = open(filename, 'r')
    npz_file = np.load(f)
    k = npz_file['k']
    t = npz_file['t']
    B = npz_file['B']
    pos = npz_file['pos']
    vel = npz_file['vel']
    lab_prob = npz_file['lab_prob']
    rot_prob = npz_file['rot_prob']
    potential = npz_file['potential']
    f.close()
    return k, t, B, pos, vel, lab_prob, rot_prob, potential


def data_crunch_cli():
    parser = argparse.ArgumentParser(description="""A differential equation
                                     solver for the archetypeal Majorana spin
                                     flip problem.""")
    parser.add_argument('--input', '-i', dest='infile',
                        default='ehrenfest_data',
                        help="""Name of the input npz file.""")
    parser.add_argument('--output', '-o', dest='outfile',
                        default='ehrenfest',
                        help="""Name of the output file.""")
    parser.add_argument('--outdir', '-d', dest='outdir',
                        default='../../gfx/',
                        help="""Name of the output file directory.""")
    args = parser.parse_args()

    print("*********************************************")
    print("* RUNNING EHRENFEST PROBLEM DATA CRUNCHER   *")
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
