# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15 14:43:24 2014

@author: jc3e13
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
# from matplotlib.colors import LogNorm
# from matplotlib.ticker import LogFormatterMathtext

import os
import sys
import glob
import itertools

lib_path = os.path.abspath('../python')
if lib_path not in sys.path:
    sys.path.append(lib_path)

import emapex
import plotting_functions as pf
import float_advection_routines as far

try:
    print("Floats {} and {} exist!".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4976)
    E76.apply_w_model('../../data/EM-APEX/4976_fix_p0k0M_fit_info.p')
    E76.apply_strain('../../data/EM-APEX/4976_N2_ref_300dbar.p')
    E76.apply_isopycnal_displacement('../../data/EM-APEX/srho_4976_100mbin.p')
    E77 = emapex.EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4977)
    E77.apply_w_model('../../data/EM-APEX/4977_fix_p0k0M_fit_info.p')
    E77.apply_strain('../../data/EM-APEX/4977_N2_ref_300dbar.p')
    E77.apply_isopycnal_displacement('../../data/EM-APEX/srho_4977_100mbin.p')


# %% Script params.

# Bathymetry file path.
bf = os.path.abspath(glob.glob('../../data/sandwell_bathymetry/topo_*.img')[0])
# Figure save path.
sdir = '../figures/parameter_space_explorer'
if not os.path.exists(sdir):
    os.makedirs(sdir)
# Universal figure font size.
matplotlib.rc('font', **{'size': 9})

# %% ##########################################################################

# Parameters to check.
Xs = np.arange(-10000., 11000., 1000.)
Ys = np.arange(-10000., 11000., 1000.)
Zs = np.arange(-6000., 6500., 500.)

params = far.default_params
pfl = E77.get_profiles(26)
use = ~np.isnan(pfl.zef) & (pfl.zef < -600.)
bscale = 250.

zf = pfl.zef[use]
uf = pfl.interp(zf, 'zef', 'U_abs')
vf = pfl.interp(zf, 'zef', 'V_abs')
wf = pfl.interp(zf, 'zef', 'Ww')
bf = bscale*pfl.interp(zf, 'z', 'b')

ps = []
cost = []

for LX, LY, LZ in itertools.product(Xs, Ys, Zs):

    ps.append((LX, LY, LZ))

    if LX == 0. or LY == 0. or LZ == 0.:
        cost.append(1e10)
        continue

    X = far.model_verbose(LX, LY, LZ, 0., params)

    um = np.interp(zf, X.r[:, 2], X.u[:, 0])
    vm = np.interp(zf, X.r[:, 2], X.u[:, 1])
    wm = np.interp(zf, X.r[:, 2], X.u[:, 2])
    bm = bscale*np.interp(zf, X.r[:, 2], X.b)

    c = np.std(um - uf) + np.std(vm - vf) + np.std(wm - wf) + np.std(bm - bf)
    cost.append(c)

    print(LX, LY, LZ, c)

#    # Create subdirectory for saving.
#    subname = "X_{:g}_Y_{:g}_Z_{:g}".format(LX, LY, LZ)
#    s2dir = os.path.join(sdir, subname)
#    if not os.path.exists(s2dir):
#        os.makedirs(s2dir)

#    fig, axs = plt.subplots(1, 5, sharey=True, figsize=(14,6))
#    axs[0].plot(1e2*pfl.U_abs, pfl.zef, 'red')
#    plt.setp(axs[0].xaxis.get_majorticklabels(), rotation=60)
#    axs[1].plot(1e2*pfl.V_abs, pfl.zef, 'red')
#    plt.setp(axs[1].xaxis.get_majorticklabels(), rotation=60)
#    axs[2].plot(1e2*pfl.Ww, pfl.z, color='red')
#    plt.setp(axs[2].xaxis.get_majorticklabels(), rotation=60)
#    axs[3].plot(1e4*pfl.b, pfl.z, 'red')
#    plt.setp(axs[3].xaxis.get_majorticklabels(), rotation=60)
#    axs[4].plot((pfl.dist_ctd - np.nanmin(pfl.dist_ctd))*1000., pfl.z, 'red')
#    axs[4].set_xlabel('$x$ (m)')
#    plt.setp(axs[4].xaxis.get_majorticklabels(), rotation=60)
#
#    axs[0].set_ylabel('$z$')
#    axs[0].plot(1e2*X.u[:, 0], X.r[:, 2])
#    axs[0].set_xlabel('$u$ (cm s$^{-1}$)')
#    plt.setp(axs[0].xaxis.get_majorticklabels(), rotation=60)
#    axs[1].plot(1e2*X.u[:, 1] - 15, X.r[:, 2])
#    axs[1].set_xlabel('$v$ (cm s$^{-1}$)')
#    plt.setp(axs[1].xaxis.get_majorticklabels(), rotation=60)
#    axs[2].plot(1e2*X.u[:, 2], X.r[:, 2])
#    axs[2].set_xlabel('$w$ (cm s$^{-1}$)')
#    plt.setp(axs[2].xaxis.get_majorticklabels(), rotation=60)
#    axs[3].plot(1e4*X.b, X.r[:, 2])
#    axs[3].set_xlabel('$b$ ($10^{-4}$ m s$^{-2}$)')
#    plt.setp(axs[3].xaxis.get_majorticklabels(), rotation=60)
#    axs[4].plot(X.r[:,0], X.r[:, 2])
#    axs[4].set_xlabel('$x$ (m)')
#    plt.setp(axs[4].xaxis.get_majorticklabels(), rotation=60)
#    plt.ylim(X.z_0, 0)
#
#    pf.my_savefig(fig, 'model', 'data', s2dir)
#    plt.close()

#    k = X.k
#    l = X.l
#    m = X.m
#    om = X.om
#    N = X.N
#    f = X.f
#    U = X.U
#    w_0 = X.w_0
#    phi_0 = X.phi_0
#    r = X.r
#    t = X.t

###############################################################################

#    xg, zg = np.meshgrid(np.arange(-1500, np.max(r[:,0]), 50), np.arange(-1600, 25, 25))
#
#    om2 = om**2
#    f2 = f**2
#    K2 = k**2 + l**2 + m**2
#    N2 = N**2
#
#    for j, ts in enumerate(np.arange(0, X.t.max(), 500.)):
#
#        idx = t.searchsorted(ts)
#        C = []
#
#        u_xg = np.real(((k*om + 1j*l*f)/(om2 - f2))*phi_0*np.exp(1j*(k*xg + m*zg - (om + k*U)*ts)))
#        u_yg = np.real(((l*om - 1j*k*f)/(om2 - f2))*phi_0*np.exp(1j*(k*xg + m*zg - (om + k*U)*ts)))
#        u_zg = np.real(((-om*K2)/((N**2 - f2)*m))*phi_0*np.exp(1j*(k*xg + m*zg - (om + k*U)*ts)))
#        bg = np.real((1j*m*N2/(N2 - om2))*phi_0*np.exp(1j*(k*xg + m*zg - (om + k*U)*ts)))
#
#        fig, axs = plt.subplots(1, 4, sharey=True, figsize=(14,6))
#        fig.suptitle('t = {:1.0f} s'.format(ts))
#        axs[0].set_ylabel('$z$ (m)')
#        C.append(axs[0].pcolormesh(xg, zg, 1e2*(u_xg + U), cmap=plt.get_cmap('bwr')))
#        axs[0].set_title('$U$ (cm s$^{-1}$)')
#
#        divider0 = make_axes_locatable(axs[0])
#        cax0 = divider0.append_axes("right", size="10%", pad=0.05)
#        plt.colorbar(C[0], cax=cax0)
#
#        C.append(axs[1].pcolormesh(xg, zg, 1e2*u_yg, cmap=plt.get_cmap('bwr')))
#        axs[1].set_title('$V$ (cm s$^{-1}$)')
#
#        divider1 = make_axes_locatable(axs[1])
#        cax1 = divider1.append_axes("right", size="10%", pad=0.05)
#        plt.colorbar(C[1], cax=cax1)
#
#        C.append(axs[2].pcolormesh(xg, zg, 1e2*u_zg, cmap=plt.get_cmap('bwr')))
#        axs[2].set_title('$W$ (cm s$^{-1}$)')
#
#        divider2 = make_axes_locatable(axs[2])
#        cax2 = divider2.append_axes("right", size="10%", pad=0.05)
#        plt.colorbar(C[2], cax=cax2)
#
#        C.append(axs[3].pcolormesh(xg, zg, 1e4*bg, cmap=plt.get_cmap('bwr')))
#        axs[3].set_title('$b$ ($10^{-4}$ m s$^{-2}$)')
#
#        divider3 = make_axes_locatable(axs[3])
#        cax3 = divider3.append_axes("right", size="10%", pad=0.05)
#        plt.colorbar(C[3], cax=cax3)
#
#        for i in xrange(4):
#            axs[i].set_xlabel('$x$ (m)')
#            axs[i].set_ylim(X.z_0, 0)
#            axs[i].set_xlim(np.min(xg), np.max(xg))
#            axs[i].plot(r[:idx, 0], r[:idx, 2], 'k--', linewidth=3)
#            axs[i].plot(r[idx, 0], r[idx, 2], 'yo', linewidth=3)
#
#        C[0].set_clim(-1e2*(X.u_0+U), 1e2*(X.u_0+U))
#        C[1].set_clim(-1e2*X.v_0, 1e2*X.v_0)
#        C[2].set_clim(-1e2*X.w_0, 1e2*X.w_0)
#        C[3].set_clim(-1e4*X.b_0, 1e4*X.b_0)
#
#        pf.my_savefig(fig, 'model', 't{:1.0f}'.format(j), s2dir)
#        plt.close()
