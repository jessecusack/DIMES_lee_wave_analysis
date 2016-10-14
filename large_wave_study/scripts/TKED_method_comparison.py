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
import gsw
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
# Estimate density measurement error.
Terr = 0.002
Tmean = 3.5
Serr = 0.002
Smean = 34.2
P = 1000.
P_ref = 1000.
lon = -60.
lat = -56.

Tr = np.random.normal(Tmean, Terr, 1000)
Sr = np.random.normal(Smean, Serr, 1000)

SAr = gsw.SA_from_SP(Sr, P, lon, lat)
CTr = gsw.CT_from_t(SAr, Tr, P)
rho1r = gsw.rho(SAr, CTr, P_ref)

#fig, axs = plt.subplots(1, 3)
#axs[0].hist(SAr - SAr.mean())
#axs[1].hist(CTr - CTr.mean())
#axs[2].hist(rho1r - rho1r.mean())

print("Standard dev rho1 {}".format(rho1r.std()))


# %%

hpid = 111
Float = E76
c = 0.146  # 4976
#c = 0.123  # 4977
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

###############################################################################
# %% Thorpe

rho_1_td, rho_1_bu, rho_1_av = TKED.intermediate_profile(rho_1, 1030., 2e-3)

R0 = 0.25
acc = 2e-3
thorpe_scales, thorpe_disp, Ls, R, rho_1s, __ = \
    TKED.thorpe_scales(-zw, rho_1, R0=R0, acc=acc, full_output=True)
L_o, L_neg, L_pos = Ls
eps_thorpe = 0.8*thorpe_scales**2 * N2_ref**(3./2.)

fig, axs = plt.subplots(1, 3, figsize=(6.5, 3), sharey='row')
axs[0].plot(rho_1, zw, label='original')
#axs[0].plot(rho_1_td, zw, label='int td')
#axs[0].plot(rho_1_bu, zw, label='int bu')
axs[0].plot(rho_1_av, zw, label='int av')
axs[0].plot(rho_1s, zw, label='sorted')
axs[0].legend(loc=0)
axs[1].plot(thorpe_disp, zw, 'yo-')
axs[1].plot(thorpe_scales, zw, 'k')
axs[1].plot(L_o, zw, 'b')
axs[1].plot(L_neg, zw, 'r')
axs[1].plot(L_pos, zw, 'g')
axs[2].plot(R, zw)
axs[2].vlines(R0, *axs[2].get_ylim())


width = 200.
binned = wdw.window(zw, eps_thorpe, width=width, overlap=0)
eps_av = np.zeros(len(binned))
z_av = np.zeros(len(binned))
for i, (z_, ep_) in enumerate(binned):
    eps_av[i] = np.trapz(ep_, z_)/width
    z_av[i] = np.mean(z_)

###############################################################################
# %% LEM
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
# %% VKE
width = 320.
overlap = width/2.
z_mid, eps_VKE = TKED.VKE_method(x, wp, width, overlap)

###############################################################################
# %% Plot
fig, axs = plt.subplots(1, 6, figsize=(6.5, 3), sharey='row')
axs[0].plot(rho_1, zw, 'k')
axs[0].plot(rho_1_av, zw, 'k:')
axs[0].plot(rho_1s, zw, 'r')
axs[1].vlines(0., *axs[2].get_ylim())
axs[1].plot(N2, zw, 'k')
axs[1].plot(N2_ref, zw, 'r')
axs[2].vlines(0., *axs[2].get_ylim())
axs[2].plot(w_filt, x, 'k')
axs[3].plot(wz, zw, 'k')
axs[3].plot(ws, zw, 'r')
axs[4].plot(thorpe_disp, zw, 'grey')
axs[4].plot(thorpe_scales, zw, 'k')
#axs[5].plot(R, zw)
#axs[5].vlines(R0, *axs[5].get_ylim())
axs[5].plot(np.log10(eps_lem_noise), x, 'grey')
axs[5].plot(np.log10(eps_thorpe), zw, 'y', linestyle='none', marker='.')
axs[5].plot(np.log10(eps_av), z_av, 'yo-')
axs[5].plot(np.log10(eps_lem), x, 'k')
axs[5].plot(np.log10(eps_VKE), z_mid, 'go-')
axs[5].grid()
