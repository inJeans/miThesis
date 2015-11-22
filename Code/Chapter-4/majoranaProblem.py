from __future__ import print_function, division

import numpy as np
import scipy.stats as stats
import os
import time

from vector import Vector

# Physical constants
G_S = 0.5
MU_B = 9.27401e-24
H_BAR = 1.05457e-34

def main():
    start = time.clock()

    # Computational parameters
    dt = 1.e-6
    # nsteps = np.ceil(0.015 / dt * 2).astype(int)
    nsteps = 100

    t = np.zeros((nsteps))
    t[0] = -0.015

    # System parameters
    Bx = 2.0e-7
    By = 0.
    C = 0.0025

    # Set initial values
    Bz = -C*t[0]
    B = Vector((Bx, By, Bz))
    magB = np.sqrt(B.x**2 + B.y**2 + B.z**2)
    n = B / magB

    psi = np.zeros((2, nsteps), dtype=complex)
    psi[:, 0] = 0.5 * np.array([[1 + n.x + n.z],
                                [1 + n.x - n.z]])[:, 0] / np.sqrt(1+n.x)

    P_up = 0.5*np.array([[1+n.z, n.x], [ n.x, 1-n.z]])
    P_dn = 0.5*np.array([[1-n.z,-n.x], [-n.x, 1+n.z]])

    prob_up = np.zeros((nsteps))
    prob_up[0] = np.dot(psi[:, 0].conj().transpose(),
                        np.dot(P_up, psi[:, 0])).real
    prob_dn = np.zeros((nsteps))
    prob_dn[0] = np.dot(psi[:, 0].conj().transpose(),
                        np.dot(P_dn, psi[:, 0])).real

    # Evolve system for nsteps
    for i in range(1, nsteps):

        t[i] = t[i-1] + dt

        B.z = -C*t[i]

        U = time_evolution_operator(B, dt)

        psi[:, i] = np.dot(U, psi[:, i-1])

        P_up = 0.5*np.array([[1+n.z, n.x], [n.x,  1-n.z]])
        P_dn = 0.5*np.array([[1-n.z,-n.x], [-n.x, 1+n.z]])

        prob_up[i] = np.dot(psi[:, i].conj().transpose(),
                            np.dot(P_up, psi[:, i])).real
        prob_dn[i] = np.dot(psi[:, i].conj().transpose(),
                            np.dot(P_dn, psi[:, i])).real

    pops = (psi.conj() * psi).real

    k = G_S*MU_B*B.x**2 / (H_BAR*C)

    print("Probability of spin flip is {0}".format(np.exp(-0.5*k*np.pi)))
    print("This simulation predicts    {0}\n".format(pops[0, -1]))

    end = time.clock()
    print("This simulation took {0}s.".format(end - start))

    f = open('no_flip_data.npz', 'w')
    np.savez(f, k=k, nz=n.z, t=t, pops=pops, prob_up=prob_up, prob_dn=prob_dn)
    f.close()


def time_evolution_operator(B, dt):
    """ This function generates the time evolution operator given the magnetic
    field
    """
    magB = np.sqrt(B.x**2 + B.y**2 + B.z**2)
    n = B / magB
    print(B)
    print("{0}, {1}, {2}".format(B.x, B.y, B.z))
    # print(magB)
    # print(n)
    # print("{0}, {1}, {2}".format(B.x/magB, B.y/magB, B.z/magB))
    # print(np.sqrt(n.x**2 + n.y**2 + n.z**2))

    theta = G_S*MU_B*magB*dt / 2. / H_BAR

    cosx = np.cos(theta)
    sinx = np.sin(theta)

    U = np.array([[cosx - 1j * n.z * sinx, -(n.y + 1j * n.x) * sinx],
                  [(n.y - 1j * n.x) * sinx, cosx + 1j * n.z * sinx]])

    return U

if __name__ == "__main__":
    main()
