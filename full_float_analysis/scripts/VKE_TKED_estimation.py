# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 15:42:45 2016

@author: jc3e13
"""

import os
import numpy as np
import scipy.signal as sig
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
from glob import glob

import emapex
import misc_data_processing as mdp
from my_savefig import my_savefig


# Figure save path.
fsdir = '../figures/VKE_TKED'
if not os.path.exists(fsdir):
    os.makedirs(fsdir)
psdir = '../processed_data/'

# Universal figure font size.
matplotlib.rc('font', **{'size': 8})

wfi_files = glob(os.path.join(psdir, '*wfi*.p'))
N2_files = glob(os.path.join(psdir, '*N2_ref*200dbar.p'))
Floats = []
for floatID in emapex.FIDS_DIMES:
    load_file = emapex.find_file(floatID)
    Float = emapex.EMApexFloat(load_file, floatID, verbose=False)

    # Ugly and probably unnecessary but should work as a check...
    if any(str(floatID) in wfi_file for wfi_file in wfi_files):
        for wfi_file in wfi_files:
            if str(floatID) in wfi_file:
                break
    else:
        raise RuntimeError("Cannot find appropriate wfi file.")

    # Ugly and probably unnecessary but should work as a check...
    if any(str(floatID) in N2_file for N2_file in N2_files):
        for N2_file in N2_files:
            if str(floatID) in N2_file:
                break
    else:
        raise RuntimeError("Cannot find appropriate N2 file.")

    Float.apply_w_model(wfi_file)
    Float.apply_strain(N2_file)

    Floats.append(Float)


# %% Adiabatic level buoyancy frequency. (Run once!)

for Float in Floats:
    mdp.adiabatic_level_float(Float, 200., psdir)


# %% Check floats pressure for noise.

hpids = np.arange(10, 60)
tres = 44.
zres = 5.

dt = 1.

tmin = 1000.
tmax = 10000.

window = 'hanning'

for Float in Floats:
    if Float.floatID == 4089:
        continue

    t = np.arange(tmin, tmax, dt)
    ng, __, pt = Float.get_interp_grid(hpids, t, 'dUTC', 'P')

    PPts = []

    for i in ng[0, :]:
        f, PPt = sig.periodogram(pt[:, i], fs=1./dt, window=window)
        PPts.append(PPt)

    PPts = np.asarray(PPts).T

    # Chop interpolation noise.
    fuse = f < 1./tres
#    muse = m < 1./zres
    f, PPts = f[fuse], PPts[fuse, :]
#    m, Pzs = m[muse], Pzs[muse, :]
    PPts[0, :] = 0.

    fig, axs = plt.subplots(2, 1, sharex='col', figsize=(3.125, 5))
#    axs[0].loglog(1./f, np.median(Pts, axis=-1), linewidth=3., alpha=0.7)
    axs[1].loglog(1./f, np.median(PPts, axis=-1), linewidth=3., alpha=0.7)
    axs[1].set_xlim(1e1, 1e4)

    axs[0].set_ylabel('Vertical kinetic energy')
    axs[1].set_ylabel('Pressure variance')
    axs[1].set_xlabel('Time period (s)')
    name = "{:g}_w_p_spec".format(Float.floatID)
    my_savefig(fig, name, fsdir, ftype='png', fsize='single_col')
    plt.close(fig)


# %% VKE plots

zmin = -2000.
dz = 1.

c = 1.  # timeeheight
x = np.arange(zmin, 0., dz)
xvar = 'timeheight'
dx = 1.
hpids = np.arange(1, 1000)
width = 20.
lc = (100., 40.)
btype = 'highpass'
we = 0.001

for Float in Floats:
    if Float.floatID == 4089:
        continue

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)

    VKE, __, = mdp.w_scales_float(Float, hpids, xvar, x, width=width, lc=lc,
                                  c=1, btype=btype, we=we, ret_noise=False)

    ieps = 0.*np.zeros_like(idxs)

    if xvar == 'time':
        __, __, iZ = Float.get_interp_grid(hpids, x, 'dUTC', 'z')
        __, __, X = Float.get_interp_grid(hpids, x, 'dUTC', 'dist_ctd')
    if xvar == 'eheight' or xvar == 'timeeheight':
        __, __, it = Float.get_interp_grid(hpids, x, 'zw', 'dUTC')
        iZ = np.zeros_like(it)
        X = np.zeros_like(it)
        for i, pfl in enumerate(Float.get_profiles(hpids)):
            iZ[:, i] = pfl.interp(it[:, i], 'dUTC', 'z')
            X[:, i] = pfl.interp(it[:, i], 'dUTC', 'dist_ctd')
    elif xvar == 'height' or xvar == 'timeheight':
        __, __, iZ = Float.get_interp_grid(hpids, x, 'z', 'z')
        __, __, X = Float.get_interp_grid(hpids, x, 'z', 'dist_ctd')

    for i in xrange(len(idxs)):
        iuse = (iZ[:, i] < -100) & (iZ[:, i] > -1400)
        # The abs function accounts for problems with z being the wrong way.
        ieps[i] = np.abs(1025.*np.trapz(VKE[iuse, i], iZ[iuse, i]))

    Z = iZ.flatten(order='F')

    use = (Z < -10) & (Z > -1400)

    Z = Z[use]

    X = X.flatten(order='F')[use]
#    noise = noise.flatten(order='F')[use]
    LOG_EPS = (np.log10(VKE)).flatten(order='F')[use]


    # Plotting #
    fig = plt.figure(figsize=(3.125, 3))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 5])
    ax0 = plt.subplot(gs[1])
    ax1 = plt.subplot(gs[0])
    ax1.plot(Float.dist[idxs], 1000.*ieps, label=Float.floatID)

    step = 10
    sc = ax0.scatter(X[::step], Z[::step], s=5, c=LOG_EPS[::step],
                     edgecolor='none', cmap=plt.get_cmap('YlOrRd'), vmin=-10.,
                     vmax=-7, alpha=.5)

    ax1.set_ylabel('$P$')
#    ax1.yaxis.set_ticks(np.array([0., 5., 10.]))
    ax1.xaxis.set_ticks([])

    ax1.legend(loc='upper right', fontsize=7)

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.7])
    C = fig.colorbar(sc, cax=cbar_ax, extend='both')
    C.set_label(r'$\log_{10}(VKE)$')

    #    plt.clim(-11., -7.)
    ax0.set_xlim(np.min(X), np.max(X))

    ax0.set_xlabel('Distance (km)')
    ax0.set_ylabel('$z$ (km)')

    ax1.set_xlim(*ax0.get_xlim())

    fontdict = {'size': 10}
    plt.figtext(-0.05, 0.85, 'a)', fontdict=fontdict)
    plt.figtext(-0.05, 0.65, 'b)', fontdict=fontdict)

    name = "{:g}_VKE_section".format(Float.floatID)
    my_savefig(fig, name, fsdir, ftype='png', fsize='single_col')
    plt.close(fig)
