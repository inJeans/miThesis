from __future__ import print_function, division

import math
import numpy as np
import scipy.stats as stats
import os
import time
import argparse

from vector import Vector

# Physical constants
G_S = 0.5
MU_B = 9.27401e-24
H_BAR = 1.05457e-34


def main(initial_time,
         Bx=2.0e-7,
         By=0.,
         C=0.0025,
         dt=1.e-6,
         filename='majorana_data'):
    start = time.clock()

    # Choose a start time that ensures the spin up direction will be 99.999%
    # aligned with the lab-frame z direction
    if not initial_time:
        initial_time = -math.sqrt((Bx**2+By**2) / (1/0.99999**2-1)) / C
    print('Start time = {0:.4} s'.format(initial_time))
    print(" ")

    # Set initial values
    nsteps = np.ceil(-initial_time / dt * 2).astype(int)
    t = np.zeros((nsteps))
    t[0] = initial_time
    Bz = -C*t[0]
    B = Vector((Bx, By, Bz))
    n = B.unit()

    print("Initialising simulation arrays")
    # Initialise wavefunction to be aligned with local magnetic field
    psi = np.zeros((2, nsteps), dtype=complex)
    psi[:, 0] = 0.5*np.array([[1 + n.x - 1j*n.y + n.z],
                              [1 + n.x + 1j*n.y - n.z]])[:, 0] / np.sqrt(1+n.x)

    P_up = spin_up_projection_operator(n)
    P_dn = spin_down_projection_operator(n)

    # Calculate inital spin projections (should be [1,0])
    rot_prob = np.zeros((2, nsteps))
    rot_prob[0, 0] = np.dot(psi[:, 0].conj().transpose(),
                            np.dot(P_up, psi[:, 0])).real
    rot_prob[1, 0] = np.dot(psi[:, 0].conj().transpose(),
                            np.dot(P_dn, psi[:, 0])).real

    print("Evolving wavefunction for {n} steps".format(n=nsteps))
    print(" ")
    # Evolve system for nsteps
    for i in range(1, nsteps):

        t[i] = t[i-1] + dt

        B.z = -C*t[i]
        n = B.unit()

        U = time_evolution_operator(B, dt)

        psi[:, i] = np.dot(U, psi[:, i-1])

        # Calculate spin projection at step n
        P_up = spin_up_projection_operator(n)
        P_dn = spin_down_projection_operator(n)

        rot_prob[0, i] = np.dot(psi[:, i].conj().transpose(),
                                np.dot(P_up, psi[:, i])).real
        rot_prob[1, i] = np.dot(psi[:, i].conj().transpose(),
                                np.dot(P_dn, psi[:, i])).real

    lab_prob = (psi.conj() * psi).real

    k = G_S*MU_B*(B.x**2+B.y**2) / (H_BAR*C)

    print("Probability of spin flip is {0}".format(np.exp(-0.5*k*np.pi)))
    print("This simulation predicts    {0}\n".format(lab_prob[0, -1]))

    end = time.clock()
    print("This simulation took {0}s.".format(end - start))

    print("Saving output file - {0}\n".format(filename+'.npz'))
    f = open(filename+'.npz', 'w')
    np.savez(f, k=k, nz=n.z, t=t, lab_prob=lab_prob, rot_prob=rot_prob)
    f.close()


def time_evolution_operator(B, dt):
    """ This function generates the time evolution operator given the magnetic
    field
    """
    n = B.unit()

    theta = G_S*MU_B*B.norm()*dt / 2. / H_BAR

    cosx = np.cos(theta)
    sinx = np.sin(theta)

    U = np.array([[cosx - 1j * n.z * sinx, -(n.y + 1j * n.x) * sinx],
                  [(n.y - 1j * n.x) * sinx, cosx + 1j * n.z * sinx]])

    return U


def spin_up_projection_operator(n):
    return 0.5*np.array([[1+n.z, n.x], [n.x, 1-n.z]])


def spin_down_projection_operator(n):
    return 0.5*np.array([[1-n.z, -n.x], [-n.x, 1+n.z]])


def majorana_cli():
    parser = argparse.ArgumentParser(description="""A differential equation
                                     solver for the archetypeal Majorana spin
                                     flip problem.""")
    parser.add_argument('--Bx', '-bx', dest='Bx', default=2.0e-7, type=float,
                        help="""A floating point value describing the strength
                                of the magnetic field in the "x" direction
                                (Tesla).""")
    parser.add_argument('--By', '-by', dest='By', default=0., type=float,
                        help="""A floating point value describing the strength
                                of the magnetic field in the "y" direction
                                (Tesla).""")
    parser.add_argument('--C', '-c', dest='C', default=0.0025, type=float,
                        help="""A floating point value describing the rate of
                                change of the magnetic field in the "z"
                                direction (Tesla/second). The field in the z
                                direction is given by Bz = -C * t, where t is
                                time (s).""")
    parser.add_argument('--delta_t', '-dt', dest='dt', default=1.e-6,
                        type=float,
                        help="""The computational time step to be used (s). The
                                optimal value for the parameter will change
                                depending on the values used for the magnetic
                                field.""")
    parser.add_argument('--initial_t', '-t0', dest='t0', type=float,
                        help="""The initial start time for the simualtion (s).
                                Keep in mind that this will need to be negative
                                if you want to simulate a zero crossing.""")
    parser.add_argument('--output', '-o', dest='filename',
                        default='majorana_data',
                        help="""Name of the output file.""")
    args = parser.parse_args()

    print("*********************************************")
    print("* RUNNING MAJORANA PROBLEM SIMULATOR        *")
    print("*********************************************")
    print(" ")
    print("Simulation parameters are:")
    print("Bx - {0} T".format(args.Bx))
    print("By - {0} T".format(args.By))
    print("C - {0} T/s".format(args.C))
    print("dt - {0} s".format(args.dt))
    print("t0 - {0} s".format(args.t0))
    print("output filename - {0}.npz".format(args.filename))
    print(" ")
    print("*********************************************")
    print(" ")

    main(initial_time=args.t0,
         Bx=args.Bx,
         By=args.By,
         C=args.C,
         dt=args.dt,
         filename=args.filename)

if __name__ == "__main__":
    majorana_cli()
