from __future__ import print_function, division

import math
import numpy as np
import scipy.stats as stats
import os
import time
import argparse

from progressbar import Bar, Percentage, ProgressBar

from vector import Vector

# Physical constants
G_S = 0.5
MU_B = 9.27401e-24
H_BAR = 1.05457e-34

M_RB87 = 1.443160648e-25

C = 2.5e-3# Approximate rate of change at the zero crossing


def main(initial_pos,
         Bx=2.0e-7,
         By=0.,
         dBz=2.5,
         dt=1.e-7,
         filename='ehrenfest_data'):
    start = time.clock()

    # Choose a start time that ensures the spin up direction will be 99.999%
    # aligned with the lab-frame z direction
    if not initial_pos:
        # initial_pos = -math.sqrt((Bx**2+By**2) / (1/0.99999**2-1)) / C
        initial_z = -np.sqrt((C**2*M_RB87**2) / (dBz**4*G_S**2*MU_B**2) -
                              2.*np.sqrt(Bx**2+By**2) * C*M_RB87 /
                             (dBz**3*G_S*MU_B))
        initial_pos = [0., 0., initial_z]
    print('Start pos = ({0:.4}, {0:.4}, {0:.4}) m'.format(initial_pos[0],
                                                          initial_pos[1],
                                                          initial_pos[2]))
    print(" ")

    # Set initial values
    # nsteps = np.ceil(-initial_time / dt * 2).astype(int)
    nsteps = int(3e4)

    print("Initialising simulation arrays")
    pos = np.ndarray(shape=(nsteps,),
                     dtype=Vector)
    pos[0] = Vector(initial_pos)
    vel = np.ndarray(shape=(nsteps,),
                     dtype=Vector)
    vel[0] = Vector([0., 0., 0.])
    acc = np.ndarray(shape=(nsteps,),
                     dtype=Vector)
    acc[0] = Vector([0., 0., 0.])

    t = np.zeros((nsteps))
    Bz = -dBz*pos[0].z
    B = np.ndarray(shape=(nsteps,),
                   dtype=Vector)
    B[0] = Vector((Bx, By, Bz))
    n = B[0].unit()

    # Initialise wavefunction to be aligned with local magnetic field
    psi = np.zeros((2, nsteps), dtype=complex)
    psi[:, 0] = 0.5*np.array([[1 + n.x - 1j*n.y + n.z],
                              [1 + n.x + 1j*n.y - n.z]])[:, 0] / np.sqrt(1+n.x)

    # Calculate inital spin projections (should be [1,0])
    P_up = spin_up_projection_operator(n)
    P_dn = spin_down_projection_operator(n)

    rot_prob = np.zeros((2, nsteps))
    rot_prob[0, 0] = np.dot(psi[:, 0].conj().transpose(),
                            np.dot(P_up, psi[:, 0])).real
    rot_prob[1, 0] = np.dot(psi[:, 0].conj().transpose(),
                            np.dot(P_dn, psi[:, 0])).real

    # Calculate intial potential energy
    V = potential_energy_operator(B[0])
    potential = np.zeros((nsteps))
    potential[0] = np.dot(psi[:, 0].conj().transpose(),
                          np.dot(V, psi[:, 0])).real

    zero_crossing = 0
    print("Evolving wavefunction for {n} steps".format(n=nsteps))
    print(" ")
    # Evolve system for nsteps
    pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=nsteps).start()
    for i in range(1, nsteps):

        t[i] = t[i-1] + dt

        acc[i] = ehrefest_force(dBz, pos[i-1], psi[:,i-1]) / M_RB87

        U = time_evolution_operator(B[i-1], dt)

        psi[:, i] = np.dot(U, psi[:, i-1])

        vel[i] = vel[i-1] + acc[i-1]*dt
        pos[i] = pos[i-1] + vel[i]*dt

        B[i] = Vector([Bx, By, -dBz*pos[i-1].z])
        n = B[i].unit()

        # Calculate spin projection at step n
        P_up = spin_up_projection_operator(n)
        P_dn = spin_down_projection_operator(n)

        rot_prob[0, i] = np.dot(psi[:, i].conj().transpose(),
                                np.dot(P_up, psi[:, i])).real
        rot_prob[1, i] = np.dot(psi[:, i].conj().transpose(),
                                np.dot(P_dn, psi[:, i])).real

        # Calculate expectation value of the potential energy
        V = potential_energy_operator(B[i])

        potential[i] = np.dot(psi[:, i].conj().transpose(),
                                np.dot(V, psi[:, i])).real

        if np.sign(pos[i-1].z) != np.sign(pos[i].z):
            zero_crossing = i-1

        pbar.update(i+1)

    pbar.finish()
    print(' ')

    lab_prob = (psi.conj() * psi).real

    k = G_S*MU_B*(Bx**2+By**2) / (H_BAR*vel[zero_crossing].z*dBz)

    print("Probability of spin flip is {0}".format(np.exp(-0.5*k*np.pi)))
    print("This simulation predicts    {0}\n".format(lab_prob[0, -1]))

    end = time.clock()
    print("This simulation took {0}s.".format(end - start))

    print("Saving output file - {0}\n".format(filename+'.npz'))
    f = open(filename+'.npz', 'w')
    np.savez(f, k=k, t=t, B=B, pos=pos, vel=vel,
             lab_prob=lab_prob, rot_prob=rot_prob,
             potential=potential)
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

def ehrefest_force(dBz, pos, psi):
    force = Vector([0., 0., 0.])

    force.x = np.dot(psi.conj().transpose(), 
                     np.dot(force_evolution_operator_x(dBz, pos), psi)).real
    force.y = np.dot(psi.conj().transpose(),
                     np.dot(force_evolution_operator_y(dBz, pos), psi)).real
    force.z = np.dot(psi.conj().transpose(), 
                     np.dot(force_evolution_operator_z(dBz, pos), psi)).real

    return force

def force_evolution_operator_x(dBz, pos):
    """ This function generates the force operator in the x direction 
    for a quadrupole potential given the position and field gradient in 
    the z direction.
    """
    Fx = -G_S*MU_B * np.array([[0., 0.],
                               [0., 0.]])

    return Fx

def force_evolution_operator_y(dBz, pos):
    """ This function generates the force operator in the y direction 
    for a quadrupole potential given the position and field gradient in 
    the z direction.
    """
    Fy = -G_S*MU_B * np.array([[0., 0.],
                               [0., 0.]])

    return Fy

def force_evolution_operator_z(dBz, pos):
    """ This function generates the force operator in the z direction 
    for a quadrupole potential given the position and field gradient in 
    the z direction.
    """
    Fz = -G_S*MU_B * np.array([[-dBz, 0.],
                               [0.,dBz]])

    return Fz


def spin_up_projection_operator(n):
    return 0.5*np.array([[1+n.z, n.x], [n.x, 1-n.z]])


def spin_down_projection_operator(n):
    return 0.5*np.array([[1-n.z, -n.x], [-n.x, 1+n.z]])

def potential_energy_operator(B):
    return -G_S * MU_B * np.array([[B.z, B.x-1j*B.y],
                                   [B.x+1j*B.y, -B.z]])


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
    parser.add_argument('--dBz', '-dbz', dest='dBz', default=2.5, type=float,
                        help="""A floating point value describing the gradient of
                                the magnetic field in the "z" direction (Tesla/m). 
                                The field in the z direction is given by Bz = -dBz * pos.z, 
                                where pos.z is the z coordinate position (m).""")
    parser.add_argument('--delta_t', '-dt', dest='dt', default=1.e-6,
                        type=float,
                        help="""The computational time step to be used (s). The
                                optimal value for the parameter will change
                                depending on the values used for the magnetic
                                field.""")
    parser.add_argument('--initial_pos', '-pos_0', dest='pos_0', type=float,
                        help="""The initial start position for the simualtion (m).
                                Keep in mind that this will need to be negative
                                if you want to simulate a zero crossing.""")
    parser.add_argument('--output', '-o', dest='filename',
                        default='ehrenfest_data',
                        help="""Name of the output file.""")
    args = parser.parse_args()

    print("*********************************************")
    print("* RUNNING EHRENFEST PROBLEM SIMULATOR       *")
    print("*********************************************")
    print(" ")
    print("Simulation parameters are:")
    print("Bx - {0} T".format(args.Bx))
    print("By - {0} T".format(args.By))
    print("dBz - {0} T/m".format(args.dBz))
    print("dt - {0} s".format(args.dt))
    print("pos_0 - {0} s".format(args.pos_0))
    print("output filename - {0}.npz".format(args.filename))
    print(" ")
    print("*********************************************")
    print(" ")

    main(initial_pos=args.pos_0,
         Bx=args.Bx,
         By=args.By,
         dBz=args.dBz,
         dt=args.dt,
         filename=args.filename)

if __name__ == "__main__":
    majorana_cli()
