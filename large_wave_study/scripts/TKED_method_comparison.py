# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 15:20:57 2016

@author: jc3e13
"""

import os
import numpy as np
import scipy.signal as sig
import matplotlib
import matplotlib.pyplot as plt

import emapex
import TKED_parameterisations as TKED
import plotting_functions as pf
import misc_data_processing as mdp
import window as wdw


try:
    print("Floats {} and {} exist!.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.load(4976)
#    E76.generate_regular_grids(zmin=zmin, dz=dz)
    E77 = emapex.load(4977)
#    E77.generate_regular_grids(zmin=zmin, dz=dz)

# %% Script params.

# Figure save path.
sdir = '../figures/TKED_estimation'
if not os.path.exists(sdir):
    os.makedirs(sdir)
# Universal figure font size.
matplotlib.rc('font', **{'size': 8})

# %%

hpid = 26
Float = E77
#c = 0.146  # 4976
c = 0.123  # 4977
# eheight
#c = 0.192  # 4976
#c = 0.159  # 4977

pfl = Float.get_profiles(hpid)
nnan = ~np.isnan(pfl.P)
rho_1 = pfl.rho_1[nnan]
z = pfl.z[nnan]
zw = pfl.zw[nnan]
N2 = pfl.N2[nnan]
N2_ref = pfl.N2_ref[nnan]
w = pfl.Ww[nnan]
wz = pfl.Wz[nnan]
ws = pfl.Ws[nnan]

rho_1s = np.sort(rho_1)[::-1]


###############################################################################
# Thorpe
thorpe_scales, thorpe_disp, x_sorted, idxs = TKED.thorpe_scales(zw, rho_1, acc=3e-3)
eps_thorpe = 0.8*thorpe_scales**2 * N2_ref**(3./2.)

width = 200.
binned = wdw.window(zw, eps_thorpe, width=width, overlap=width/2)
eps_av = np.zeros(len(binned))
z_av = np.zeros(len(binned))
for i, (z_, ep_) in enumerate(binned):
    eps_av[i] = np.trapz(ep_, z_)/width
    z_av[i] = np.mean(z_)

###############################################################################
# LEM
zmin = np.ceil(np.min(zw))
zmax = np.floor(np.max(zw))
dz = 1.

x = np.arange(zmin, zmax, dz)
####
overlap = -1

#xvar = 'eheight'
#dx = 1.
#width = 10.
#lc = 30.
#btype = 'highpass'
#we = 0.001
###
xvar = 'timeeheight'
dx = 1.
width = 20.
lc = (100., 40.)
c = 1.
btype = 'highpass'
we = 0.001


if xvar == 'time':
    __, __, wp = Float.get_interp_grid(hpid, x, 'dUTC', 'Ww')
    __, __, N2p = Float.get_interp_grid(hpid, x, 'dUTC', 'N2_ref')
elif xvar == 'height':
    __, __, wp = Float.get_interp_grid(hpid, x, 'z', 'Ww')
    __, __, N2p = Float.get_interp_grid(hpid, x, 'z', 'N2_ref')
elif xvar == 'eheight':
    __, __, wp = Float.get_interp_grid(hpid, x, 'zw', 'Ww')
    __, __, N2p = Float.get_interp_grid(hpid, x, 'zw', 'N2_ref')
elif (xvar == 'timeheight') or (xvar == 'timeeheight'):
    # First low-pass in time.
    dt = 1.
    t = np.arange(0., 15000., dt)
    __, __, wt = Float.get_interp_grid(hpid, t, 'dUTC', 'Ww')
    xc = 1./lc[0]  # cut off wavenumber
    normal_cutoff = xc*dt*2.  # Nyquist frequency is half 1/dx.
    b, a = sig.butter(4, normal_cutoff, btype='lowpass')
    wf = sig.filtfilt(b, a, wt, axis=0)

    if xvar == 'timeheight':
        __, __, N2p = Float.get_interp_grid(hpid, x, 'z', 'N2_ref')
        __, __, it = Float.get_interp_grid(hpid, x, 'z', 'dUTC')
    elif xvar == 'timeeheight':
        __, __, N2p = Float.get_interp_grid(hpid, x, 'zw', 'N2_ref')
        __, __, it = Float.get_interp_grid(hpid, x, 'zw', 'dUTC')

    wp = np.interp(it, t, wf)

    btype = 'highpass'
    lc = lc[1]
else:
    raise ValueError("xvar must either be 'time', 'height', 'eheight' or "
                     "'timeheight'.")

xc = 1./lc  # cut off wavenumber
normal_cutoff = xc*dx*2.  # Nyquist frequency is half 1/dx.
b, a = sig.butter(4, normal_cutoff, btype=btype)
w_filt = sig.filtfilt(b, a, wp)

eps_lem, __, eps_lem_noise, noise_flag = \
    TKED.w_scales(wp, x, N2p, dx, width, overlap, lc, c, 0.2, btype, we, True)

###############################################################################
# VKE
width = 320.
overlap = width/2.
z_mid, eps_VKE = TKED.VKE_method(x, wp, width, overlap)

###############################################################################
fig, axs = plt.subplots(1, 6, sharey='row')
axs[0].plot(rho_1, zw, 'k')
axs[0].plot(rho_1s, zw, 'r')
axs[1].vlines(0., *axs[2].get_ylim())
axs[1].plot(N2, zw, 'k')
axs[1].plot(N2_ref, zw, 'r')
axs[2].vlines(0., *axs[2].get_ylim())
axs[2].plot(w_filt, x, 'k')
axs[3].plot(wz, zw, 'k')
axs[3].plot(ws, zw, 'r')
axs[4].plot(thorpe_scales, zw, 'k')
axs[5].semilogx(eps_lem_noise, x, 'grey')
axs[5].semilogx(eps_thorpe, zw, 'y', linestyle='none', marker='.')
axs[5].semilogx(eps_av, z_av, 'yo-')
axs[5].semilogx(eps_lem, x, 'k')
axs[5].semilogx(eps_VKE, z_mid, 'go-')