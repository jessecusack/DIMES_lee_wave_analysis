# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 15:15:52 2016

@author: jc3e13
"""

import os
import glob
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import gsw

import emapex
from my_savefig import my_savefig
import gravity_waves as gw

try:
    print("Floats {} and {} exist!.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.load(4976)
    E77 = emapex.load(4977)
    E76.calculate_pressure_perturbation()
    E76.update_profiles()
    E77.calculate_pressure_perturbation()
    E77.update_profiles()

# %% Script params.

# Bathymetry file path.
bf = os.path.abspath(glob.glob('/noc/users/jc3e13/storage/smith_sandwell/topo_*.img')[0])
# Figure save path.
sdir = '../figures/animate'
if not os.path.exists(sdir):
    os.makedirs(sdir)
# Universal figure font size.
matplotlib.rc('font', **{'size': 9})

# %% Show...


def w_model(params, data):

    phi_0, X, Y, Z, phase_0 = params

    time, x, y, z, U, V, N, f = data

    k = 2*np.pi/X
    l = 2*np.pi/Y
    m = 2*np.pi/Z

    om = gw.omega(N, k, m, l, f)

    w = gw.w(x, y, z, time, phi_0, k, l, m, om, N, U=U, V=V, phase_0=phase_0)

    return w


pfl = E77.get_profiles(26)
x = pfl.x_ctd
nans = np.isnan(x)
x = x[~nans]
y = pfl.y_ctd[~nans]
z = pfl.z[~nans]
t = pfl.dUTC[~nans]
w = pfl.Ww[~nans]

# Wave parameters
X = -3974.57854378
Y = -1263.47814568
Z = -1809.46862232
phi_0 = 0.0332719148841
phase_0 = -3.2564615625377127

zmin = np.min(z)
zmax = -650.

U = np.mean(pfl.U_abs[pfl.zef < zmax])
V = np.mean(pfl.V_abs[pfl.zef < zmax])
N = np.nanmean(np.sqrt(pfl.N2))
f = gsw.f(pfl.lat_start)

xmin = np.min(x)
xmax = np.max(x)
xg, zg = np.meshgrid(np.arange(xmin, xmax, 200),
                     np.arange(zmin, zmax, 25))

for j, ts in enumerate(np.arange(10, 6000., 100.)):

    w_ = w[t < ts]
    x_ = x[t < ts]
    y_ = y[t < ts]
    z_ = z[t < ts]
    t_ = t[t < ts]

    yts = np.interp(ts, t, y)

    datag = [ts, xg, yts, zg, U, V, N, f]
    params = [phi_0, X, Y, Z, phase_0]
    wg = w_model(params, datag)
    datal = [t_, x_, y_, z_, U, V, N, f]
    wl = w_model(params, datal)

    fig = plt.figure(figsize=(3.125, 3))
    gs = gridspec.GridSpec(1, 3, width_ratios=[6, 6, 1])
    axs = np.array([plt.subplot(g) for g in gs])

    axs[0].set_ylim(zmin, zmax)
    axs[0].set_ylabel('$z$ (m)')

    fig.suptitle('t = {:1.0f} s'.format(ts))

    axs[0].plot(w_, z_, 'r', linewidth=2)
    axs[0].plot(wl, z_, 'k--', linewidth=2)
    axs[0].set_xlim(-0.3, 0.3)
    axs[0].set_xticks([-0.2, 0., 0.2])
    axs[0].set_xlabel('$w$ (m s$^{-1}$)')
    axs[0].set_ylim(zmin, zmax)
    axs[0].set_ylabel('$z$ (m)')

    C = axs[1].pcolormesh(xg, zg, wg, cmap=plt.get_cmap('bwr'))
    axs[1].plot(x_, z_, 'k--', linewidth=2)
    axs[1].plot(x_[-1], z_[-1], 'yo', markersize=6)
    axs[1].set_xlim(xmin, xmax)
    axs[1].set_xticks([0., 2000., 4000.] + x_[0])
    axs[1].set_xticklabels(['0', '2000', '4000'])
    axs[1].set_xlabel('Zonal distance (m)')
    axs[1].set_ylim(zmin, zmax)
    axs[1].set_yticklabels([])

    plt.colorbar(C, cax=axs[2], label='$w$ (m s$^{-1}$)')

    name = "t{:1.0f}".format(ts)
    my_savefig(fig, '4977_26', name, sdir, dtype='png')
    plt.close(fig)
