# -*- coding: utf-8 -*-
"""
Created on Fri Jan 23 16:25:26 2015

@author: jc3e13
"""

import sys
import os
import glob
import scipy as sp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
import scipy.signal as sig
from scipy.integrate import trapz

import gsw

lib_path = os.path.abspath('../modules')
if lib_path not in sys.path:
    sys.path.append(lib_path)

import emapex
import finescale as fs
import plotting_functions as pf
import sandwell
import window as wdw

zmin = -1500.
dz = 1.

try:
    print("Floats {} and {} exist!.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.load(4976)
    E76.generate_regular_grids(zmin=zmin, dz=dz)
    E77 = emapex.load(4977)
    E77.generate_regular_grids(zmin=zmin, dz=dz)

# %% Script params.

# Bathymetry file path.
bf = os.path.abspath(glob.glob('/noc/users/jc3e13/storage/smith_sandwell/topo_*.img')[0])
# Figure save path.
sdir = '../figures/TKED_estimation'
if not os.path.exists(sdir):
    os.makedirs(sdir)
# Universal figure font size.
matplotlib.rc('font', **{'size': 9})


# %% Experimental functions

def w_scales_experimental(w, z, N2, dz=5., c=0.5, eff=0.2, lc=30.):
    """Inputs should be regularly spaced."""

    # First we have to design the high pass filter the data. Beaird et. al.
    # 2012 use a forth order butterworth with a cutoff of 30m.
    mc = 1./lc  # cut off wavenumber (m-1)
    normal_cutoff = mc*dz*2.  # Nyquist frequency is half 1/dz.
    b, a = sig.butter(4, normal_cutoff, btype='highpass')

    # Filter the data.
    w_filt = sig.lfilter(b, a, w)

    w_wdws = wdw.window(z, w_filt, width=10., overlap=-1.)
    N2_wdws = wdw.window(z, N2, width=10., overlap=-1.)

    w_rms = np.zeros_like(z)
    N2_mean = np.zeros_like(z)

    for i, (w_wdw, N2_wdw) in enumerate(zip(w_wdws, N2_wdws)):
        w_rms[i] = np.std(w_wdw[1])
        N2_mean[i] = np.mean(N2_wdw[1])

    epsilon = c*np.sqrt(N2_mean)*w_rms**2
    kappa = eff*epsilon/N2_mean

    return epsilon, kappa


def w_scales_float_experimental(Float, hpids, c=0.5, eff=0.2, lc=30.):

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)

    w = Float.r_Ww[:, idxs]
    z = Float.r_z[:, 0]
    N2 = Float.r_N2_ref[:, idxs]

    dz = z[0] - z[1]

    epsilon = np.zeros_like(w)
    kappa = np.zeros_like(w)

    for i, (w_row, N2_row) in enumerate(zip(w.T, N2.T)):
        epsilon[:, i], kappa[:, i] = w_scales_experimental(w_row, z, N2_row,
                                                           dz, c, eff, lc)

    return epsilon, kappa


# %% Coefficient estimation using microstructure data

UK2_vmp = sp.io.loadmat('../../storage/DIMES/combined_jc054.mat',
                        variable_names=['vmp'])['vmp']

z_vmp = gsw.z_from_p(UK2_vmp['press'][0][0], UK2_vmp['startlat'][0][0][0])

fig1 = plt.figure()
plt.plot(UK2_vmp['startlon'][0][0][0], UK2_vmp['startlat'][0][0][0], 'ko')
plt.plot(UK2_vmp['startlon'][0][0][0][25:30], UK2_vmp['startlat'][0][0][0][25:30], 'ro')
plt.plot(E76.lon_start[:100], E76.lat_start[:100])
plt.plot(E77.lon_start[:100], E77.lat_start[:100])

z = np.arange(zmin, 0., dz)

hpids = np.arange(1, 100)

epsilon_76, kappa_76 = w_scales_float_experimental(E76, hpids, c=0.03)
epsilon_77, kappa_77 = w_scales_float_experimental(E77, hpids, c=0.04)

fig2 = plt.figure()
plt.semilogx(UK2_vmp['eps'][0][0][:, 25:30], z_vmp[:, 25:30], color='grey',
             alpha=0.2)
plt.semilogx(epsilon_76, z, color='red', alpha=0.2)
plt.semilogx(epsilon_77, z, color='red', alpha=0.2)

epsilon_vmp = UK2_vmp['eps'][0][0][:, 25:30].flatten()
z_vmp_flat = z_vmp[:, 25:30].flatten()
use = ~np.isnan(epsilon_vmp) & (z_vmp_flat > zmin)

bins = np.arange(-12., -5, 0.25)
fig, axs = plt.subplots(3, 1, sharex='col', sharey=True, figsize=(3, 6))
axs[0].hist(np.log10(epsilon_vmp[use]), bins=bins, color='blue', alpha=0.8,
            normed=True, label='VMP')
axs[1].hist(np.log10(epsilon_76.flatten()), bins=bins, color='red', alpha=0.8,
            normed=True, label='4976')
axs[2].hist(np.log10(epsilon_77.flatten()), bins=bins, color='red', alpha=0.8,
            normed=True, label='4977')
axs[2].set_xlabel('$\log_{10}(\epsilon)$ W kg$^{-1}$')

axs[1].text(-7, 0.6, '$c = 0.03$')
axs[2].text(-7, 0.6, '$c = 0.04$')

for ax in axs:
    ax.legend()

pf.my_savefig(fig, 'vmp', 'comparison', sdir, ftype='pdf', fsize='single_col')

# %% Start script
# Using large eddy method first.
# Different coefficient for each float.
cs = [0.03, 0.04]

for Float, c in zip([E76, E77], cs):

    hpids = np.arange(1, 100)
    __, idxs = Float.get_profiles(hpids, ret_idxs=True)

    bathy = sandwell.interp_track(Float.lon_start, Float.lat_start, bf)
    epsilon, kappa = w_scales_float_experimental(Float, hpids, c=c)

    ieps = 0.*np.zeros_like(idxs)

    iZ = Float.r_z[:, 0]
    iuse = (iZ < -100) & (iZ > -1400)

    for i in xrange(len(idxs)):
        ieps[i] = 1025.*trapz(epsilon[iuse, i], iZ[iuse])

    Z = (Float.r_z[:, idxs]).flatten()

    use = (Z < -100) & (Z > -1400)

    Z = Z[use]
    X = (Float.r_dist_ctd[:, idxs]).flatten()[use]
    LOG_EPS = (np.log10(epsilon)).flatten()[use]
    LOG_KAP = (np.log10(kappa)).flatten()[use]

    # Plotting #
    # Epsilon
    fig = plt.figure(figsize=(10, 4))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 5])
    ax0 = plt.subplot(gs[1])
    ax1 = plt.subplot(gs[0])

    ax1.plot(Float.dist[idxs], 1000.*ieps, color='black')
    ax1.set_ylabel('$P$ (mW m$^{-2}$)')
    ax1.yaxis.set_ticks(np.array([0., 10., 20.]))
    ax1.xaxis.set_ticks([])

    sc = ax0.scatter(X, Z, s=5, c=LOG_EPS, edgecolor='none',
                     cmap=plt.get_cmap('YlOrRd'), vmin=-11., vmax=-7)

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.7])
    C = fig.colorbar(sc, cax=cbar_ax, extend='both')
    C.set_label(r'$\log_{10}(\epsilon)$ (W kg$^{-1}$)')

#    plt.clim(-11., -7.)
    ax0.set_xlim(np.min(X), np.max(X))
    ax0.set_ylim(-5000., 0.)
    ax0.set_xlabel('Distance (km)')
    ax0.set_ylabel('$z$ (m)')
    ax0.fill_between(Float.dist, bathy, -5000., color='black', linewidth=2)

    ax1.set_xlim(*ax0.get_xlim())

    pf.my_savefig(fig, Float.floatID, 'epsilon_lem', sdir, fsize='double_col')

    ##
    # Kappa
    fig = plt.figure(figsize=(10, 4))
    plt.scatter(X, Z, s=5, c=LOG_KAP, edgecolor='none', cmap=plt.get_cmap('YlOrRd'))
    C = plt.colorbar(extend='both')
    C.set_label(r'$\log_{10}(\kappa_\rho)$ (m$^2$ s$^{-1}$)')
    plt.clim(-5., -3.)
    plt.xlim(np.min(X), np.max(X))
    plt.ylim(-5000., 0.)
    plt.xlabel('Distance (km)')
    plt.ylabel('$z$ (m)')
    plt.fill_between(Float.dist, bathy, -5000., color='black', linewidth=2)
    pf.my_savefig(fig, Float.floatID, 'kappa_lem', sdir, fsize='double_col')

    # U V W

    __, d_ef = Float.get_timeseries(hpids, 'dist_ef')
    __, d_ctd = Float.get_timeseries(hpids, 'dist_ctd')
    __, z = Float.get_timeseries(hpids, 'z')
    __, zef = Float.get_timeseries(hpids, 'zef')
    __, u = Float.get_timeseries(hpids, 'U_abs')
    __, v = Float.get_timeseries(hpids, 'V_abs')
    __, w = Float.get_timeseries(hpids, 'Ww')
    __, N2_ref = Float.get_timeseries(hpids, 'N2_ref')

    u[np.abs(u) > 1.5] = np.NaN
    v[np.abs(v) > 1.5] = np.NaN

    fig = plt.figure(figsize=(10, 4))
    plt.scatter(d_ef, zef, s=5., c=u, edgecolor='none', cmap=plt.get_cmap('bwr'))
    C = plt.colorbar(extend='both')
    C.set_label(r'$u$ (m s$^{-1}$)')
    plt.clim(-1., 1.)
    plt.xlim(np.min(X), np.max(X))
    plt.ylim(-5000., 0.)
    plt.xlabel('Distance (km)')
    plt.ylabel('$z$ (m)')
    plt.fill_between(Float.dist, bathy, -5000., color='black', linewidth=2)
    pf.my_savefig(fig, Float.floatID, 'u_rel', sdir, fsize='double_col')

    fig = plt.figure(figsize=(10, 4))
    plt.scatter(d_ef, zef, s=5., c=v, edgecolor='none', cmap=plt.get_cmap('bwr'))
    C = plt.colorbar(extend='both')
    C.set_label(r'$v$ (m s$^{-1}$)')
    plt.clim(-1., 1.)
    plt.xlim(np.min(X), np.max(X))
    plt.ylim(-5000., 0.)
    plt.xlabel('Distance (km)')
    plt.ylabel('$z$ (m)')
    plt.fill_between(Float.dist, bathy, -5000., color='black', linewidth=2)
    pf.my_savefig(fig, Float.floatID, 'v_rel', sdir, fsize='double_col')

    fig = plt.figure(figsize=(10, 4))
    plt.scatter(d_ctd, z, s=5., c=w, edgecolor='none', cmap=plt.get_cmap('bwr'))
    C = plt.colorbar(extend='both')
    C.set_label(r'$w$ (m s$^{-1}$)')
    plt.clim(-0.1, 0.1)
    plt.xlim(np.min(X), np.max(X))
    plt.ylim(-5000., 0.)
    plt.xlabel('Distance (km)')
    plt.ylabel('$z$ (m)')
    plt.fill_between(Float.dist, bathy, -5000., color='black', linewidth=2)
    pf.my_savefig(fig, Float.floatID, 'w', sdir, fsize='double_col')

    fig = plt.figure(figsize=(10, 4))
    plt.scatter(d_ctd, z, s=5., c=np.sqrt(N2_ref), edgecolor='none', cmap=plt.get_cmap('cool'))
    C = plt.colorbar(extend='both')
    C.set_label(r'$N$ (rad s$^{-1}$)')
    plt.clim(0.001, 0.003)
    plt.xlim(np.min(X), np.max(X))
    plt.ylim(-5000., 0.)
    plt.xlabel('Distance (km)')
    plt.ylabel('$z$ (m)')
    plt.fill_between(Float.dist, bathy, -5000., color='black', linewidth=2)
    pf.my_savefig(fig, Float.floatID, 'N2_ref', sdir, fsize='double_col')

# %% Using Thorpe scales
fs.thorpe_scales()

# %% Using finescale parameterisation

params = fs.default_params

params['plot_results'] = False
params['plot_profiles'] = False
params['plot_spectra'] = False
params['print_diagnostics'] = False
params['periodogram_params']['nfft'] = None
params['periodogram_params']['window'] = 'hanning'
params['m_0'] = 1./120.
params['m_c'] = 1./12.
params['bin_width'] = 200.
params['bin_overlap'] = 100.
params['apply_corrections'] = True
params['zmin'] = -1400
params['zmax'] = -100

E76_hpids = np.arange(1, 350)
E77_hpids = np.arange(1, 250)

for Float, hpids in zip([E76, E77], [E76_hpids, E77_hpids]):

    results = fs.analyse_float(Float, hpids, params)
    __, idxs = Float.get_profiles(hpids, ret_idxs=True)
    dists = Float.dist[idxs]

    z_meang = []
    EKg = []
    R_polg = []
    R_omg = []
    epsilong = []
    kappag = []

    for result in results:
        z_mean, EK, R_pol, R_om, epsilon, kappa = result

        z_meang.append(z_mean)
        EKg.append(EK)
        R_polg.append(R_pol)
        R_omg.append(R_om)
        epsilong.append(epsilon)
        kappag.append(kappa)

    z_meang = np.flipud(np.transpose(np.asarray(z_meang)))
    EKg = np.flipud(np.transpose(np.asarray(EKg)))
    R_polg = np.flipud(np.transpose(np.asarray(R_polg)))
    R_omg = np.flipud(np.transpose(np.asarray(R_omg)))
    epsilong = np.flipud(np.transpose(np.asarray(epsilong)))
    kappag = np.flipud(np.transpose(np.asarray(kappag)))
    # Plotting #

    ylims = (params['zmin'], params['zmax'])

    fig = plt.figure()
    plt.title('log10 R_pol')
    plt.pcolormesh(dists, z_meang[:,0], np.log10(R_polg), cmap=plt.get_cmap('bwr'))
    plt.clim(-1, 1)
    plt.colorbar()
    plt.ylim(ylims)
    plt.xlim(np.min(dists), np.max(dists))
    pf.my_savefig(fig, Float.floatID, 'R_pol_fs', sdir)

    fig = plt.figure()
    plt.title('log10 epsilon')
    plt.pcolormesh(dists, z_meang[:,0], np.log10(epsilong), cmap=plt.get_cmap('bwr'))
    plt.clim(-11, -7)
    plt.colorbar()
    plt.ylim(ylims)
    plt.xlim(np.min(dists), np.max(dists))
    pf.my_savefig(fig, Float.floatID, 'epsilon_fs', sdir)

#    plt.figure()
#    plt.title('log10 kappa')
#    for result, dist in zip(results, dists):
#        z_mean, EK, R_pol, R_om, epsilon, kappa = result
#        d = dist*np.ones_like(z_mean)
#        plt.scatter(d, z_mean, c=np.log10(kappa), edgecolor='none',
#                    cmap=plt.get_cmap('jet'))
#
#    plt.colorbar()
#    plt.ylim(ylims)
#    plt.xlim(np.min(dists), np.max(dists))
#    plt.savefig('../figures/finescale/kappa.png', bbox_inches='tight')
#
#    plt.figure()
#    plt.title('R_om')
#    for result, dist in zip(results, dists):
#        z_mean, EK, R_pol, R_om, epsilon, kappa = result
#        d = dist*np.ones_like(z_mean)
#        plt.scatter(d, z_mean, c=R_om, edgecolor='none',
#                    cmap=plt.get_cmap('jet'))
#
#    plt.colorbar()
#    plt.ylim(ylims)
#    plt.xlim(np.min(dists), np.max(dists))
#    plt.savefig('../figures/finescale/R_om.png', bbox_inches='tight')
#
#    plt.figure()
#    plt.title('log10 EK')
#    for result, dist in zip(results, dists):
#        z_mean, EK, R_pol, R_om, epsilon, kappa = result
#        d = dist*np.ones_like(z_mean)
#        plt.scatter(d, z_mean, c=np.log10(EK), edgecolor='none',
#                    cmap=plt.get_cmap('jet'))
#
#    plt.colorbar()
#    plt.ylim(ylims)
#    plt.xlim(np.min(dists), np.max(dists))
#    plt.savefig('../figures/finescale/EK.png', bbox_inches='tight')