# -*- coding: utf-8 -*-
"""
Created on Mon Dec 29 10:10:37 2014

@author: jc3e13
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys

lib_path = os.path.abspath('../modules')
if lib_path not in sys.path:
    sys.path.append(lib_path)

import gravity_waves as gw


def u2D(x, z, t, k, m, om, u_0, U, phase_0=0.):
    phase = k*x + m*z - (om + U*k)*t + phase_0
    return np.real(u_0*np.exp(1j*phase))


def main():

    xg, zg = np.meshgrid(np.arange(-2000., 2000., 10),
                         np.arange(-3000., 0., 10.))

    N = 2e-3
    k = -2.*np.pi/2000.
    m = -2.*np.pi/2000.

    om = gw.omega(N, k, m)

    U = 0.0 # -om/k
    u_0 = 0.1

    phase_0 = 0. # np.pi

#    for phase_0 in np.linspace(0., 2*np.pi, 10):
#
#        ug = u2D(xg, zg, 0., k, m, om, u_0, U, phase_0)
#
#        plt.figure()
#        plt.pcolormesh(xg, zg, ug, cmap=plt.get_cmap('bwr'))
#        plt.colorbar()
#        plt.title("{:1.2f} rad".format(phase_0))

    for t in np.arange(0., 10000., 500.):

        ug = u2D(xg, zg, t, k, m, om, u_0, U, phase_0)

        plt.figure()
        plt.pcolormesh(xg, zg, ug, cmap=plt.get_cmap('bwr'))
        plt.colorbar()
        plt.title("{:g} s".format(t))


if __name__ == "__main__":
    main()


#    fig, axs = plt.subplots(1, 4, sharey=True, figsize=(14,6))
#    fig.suptitle('t = {:1.0f} s'.format(ts))
#    axs[0].set_ylabel('$z$ (m)')
#    C.append(axs[0].pcolormesh(xg, zg, 1e2*(u_xg + U), cmap=plt.get_cmap('bwr')))
#    axs[0].set_title('$U$ (cm s$^{-1}$)')
#
#    divider0 = make_axes_locatable(axs[0])
#    cax0 = divider0.append_axes("right", size="10%", pad=0.05)
#    plt.colorbar(C[0], cax=cax0)
#
#    C.append(axs[1].pcolormesh(xg, zg, 1e2*u_yg, cmap=plt.get_cmap('bwr')))
#    axs[1].set_title('$V$ (cm s$^{-1}$)')
#
#    divider1 = make_axes_locatable(axs[1])
#    cax1 = divider1.append_axes("right", size="10%", pad=0.05)
#    plt.colorbar(C[1], cax=cax1)
#
#    C.append(axs[2].pcolormesh(xg, zg, 1e2*u_zg, cmap=plt.get_cmap('bwr')))
#    axs[2].set_title('$W$ (cm s$^{-1}$)')
#
#    divider2 = make_axes_locatable(axs[2])
#    cax2 = divider2.append_axes("right", size="10%", pad=0.05)
#    plt.colorbar(C[2], cax=cax2)
#
#    C.append(axs[3].pcolormesh(xg, zg, 1e4*bg, cmap=plt.get_cmap('bwr')))
#    axs[3].set_title('$b$ ($10^{-4}$ m s$^{-2}$)')
#
#    divider3 = make_axes_locatable(axs[3])
#    cax3 = divider3.append_axes("right", size="10%", pad=0.05)
#    plt.colorbar(C[3], cax=cax3)
#
#    for i in xrange(4):
#        axs[i].set_xlabel('$x$ (m)')
#        axs[i].set_ylim(X.z_0, 0)
#        axs[i].set_xlim(np.min(xg), np.max(xg))
#
#        axs[i].plot(r[:idx, 0], r[:idx, 2], 'k--', linewidth=3)
#        axs[i].plot(r[idx, 0], r[idx, 2], 'yo', linewidth=3)
#
#
#    C[0].set_clim(-1e2*(X.u_0+U), 1e2*(X.u_0+U))
#    C[1].set_clim(-1e2*X.v_0, 1e2*X.v_0)
#    C[2].set_clim(-1e2*X.w_0, 1e2*X.w_0)
#    C[3].set_clim(-1e4*X.b_0, 1e4*X.b_0)