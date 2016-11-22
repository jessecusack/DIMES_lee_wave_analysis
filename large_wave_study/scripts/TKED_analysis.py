# -*- coding: utf-8 -*-
"""
Created on Fri Jan 23 16:25:26 2015

@author: jc3e13
"""

import os
import glob
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.integrate import trapz
from scipy.stats import binned_statistic

import emapex
import misc_data_processing as mdp
import TKED_parameterisations as TKED
import plotting_functions as pf  # my_savefig
import sandwell

zmin = -1450.
dz = 1.

try:
    print("Floats {} and {} exist!.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.load(4976)
#    E76.generate_regular_grids(zmin=zmin, dz=dz)
    E77 = emapex.load(4977)
#    E77.generate_regular_grids(zmin=zmin, dz=dz)

# %% Script params.

# Bathymetry file path.
bf = os.path.abspath(glob.glob('../../../storage/smith_sandwell/topo_*.img')[0])
# Figure save path.
sdir = '../figures/TKED_estimation'
if not os.path.exists(sdir):
    os.makedirs(sdir)
# Universal figure font size.
matplotlib.rc('font', **{'size': 8})


# %% Start script
# Using large eddy method first.
# Different coefficient for each float.
#cs = [0.197, 0.158]  # time
#xvar = 'time'
#dx = 5.
#x = np.arange(0., 12000., dx)
#width = 120.
#lc = np.array([300., 120.])
#btype = 'bandpass'

#cs = [0.193, 0.160]  # height
#xvar = 'height'
#dx = 1.
#x = np.arange(-1500., 0., dx)
#width = 15.
#lc = np.array([40., 15.])
#btype = 'bandpass'

#cs = [0.176, 0.147] # timeheight
#xvar = 'timeheight'
#dx = 1.
#x = np.arange(-1450., -50, dx)
#width = 15.
#lc = np.array([100., 40.])
#btype = 'highpass'

#cs = [0.192, 0.159]  # eheight
#x = np.arange(zmin, 0., dz)
#xvar = 'eheight'
#dx = 1.
#width = 20.
#lc = np.array([40., 15.])
#btype = 'bandpass'

cs = [0.146, 0.123]  # timeeheight
x = np.arange(zmin, 0., dz)
xvar = 'timeeheight'
dx = 1.
hpids = np.arange(10, 52)
width = 20.
lc = (100., 40.)
btype = 'highpass'
we = 0.001

fig = plt.figure(figsize=(3.125, 3))
gs = gridspec.GridSpec(2, 1, height_ratios=[1, 5])
ax0 = plt.subplot(gs[1])
ax1 = plt.subplot(gs[0])

for Float, c in zip([E76, E77], cs):

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)

    epsilon, kappa, __, noise = mdp.w_scales_float(Float, hpids, xvar, x,
                                                   width=width, overlap=-1.,
                                                   lc=lc, c=c, btype=btype,
                                                   we=we, ret_noise=True)

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
        ieps[i] = np.abs(1025.*trapz(epsilon[iuse, i], iZ[iuse, i]))

    Z = iZ.flatten(order='F')

    use = (Z < -10) & (Z > -1400)

    Z = Z[use]

    X = X.flatten(order='F')[use]
    noise = noise.flatten(order='F')[use]
    LOG_EPS = (np.log10(epsilon)).flatten(order='F')[use]
    LOG_KAP = (np.log10(kappa)).flatten(order='F')[use]

    # Plotting #
    # Epsilon
    d = getattr(Float, 'dist_ctd')[:, idxs].flatten(order='F')

    tgps = getattr(Float, 'UTC_start')[idxs]
    lon = getattr(Float, 'lon_start')[idxs]
    lat = getattr(Float, 'lat_start')[idxs]
    tctd = getattr(Float, 'UTC')[:, idxs].flatten(order='F')
    nans = np.isnan(d) | np.isnan(tctd)
    tctd = tctd[~nans]
    dctd = d[~nans]
    lonctd = np.interp(tctd, tgps, lon)
    latctd = np.interp(tctd, tgps, lat)
    bathy = sandwell.interp_track(lonctd, latctd, bf)

    dbathymax = dctd[bathy.argmax()]

    dctd -= dbathymax
    X -= dbathymax

    ax1.plot(Float.dist[idxs] - dbathymax, 1000.*ieps, label=Float.floatID)

    LOG_EPS[noise] = np.NaN

    step = 10
    sc = ax0.scatter(X[::step], Z[::step], s=5, c=LOG_EPS[::step],
                     edgecolor='none', cmap=plt.get_cmap('YlOrRd'), vmin=-10.,
                     vmax=-7, alpha=.5)

ax1.set_ylabel('$P$ (mW m$^{-2}$)')
ax1.yaxis.set_ticks(np.array([0., 5., 10.]))
ax1.xaxis.set_ticks([])

ax1.legend(loc='upper right', fontsize=7)

ax0.fill_between(dctd[::100], bathy[::100],
                 np.nanmin(bathy), color='black', linewidth=2)
ax0.set_ylim(-4000., 0.)
ax0.yaxis.set_ticks(np.arange(-4000, 1000, 1000))
ax0.yaxis.set_ticklabels(['-4', '-3', '-2', '-1', '0'])

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.7])
C = fig.colorbar(sc, cax=cbar_ax, extend='both')
C.set_label(r'$\log_{10}(\epsilon)$ (W kg$^{-1}$)')

#    plt.clim(-11., -7.)
ax0.set_xlim(np.min(X), np.max(X))

ax0.set_xlabel('Distance from ridge top (km)')
ax0.set_ylabel('$z$ (km)')

ax1.set_xlim(*ax0.get_xlim())

fontdict={'size': 10}
plt.figtext(-0.05, 0.85, 'a)', fontdict=fontdict)
plt.figtext(-0.05, 0.65, 'b)', fontdict=fontdict)

pf.my_savefig(fig, 'both', 'epsilon_lem', sdir, ftype=('png', 'pdf'),
              fsize='single_col')


# %% Using Thorpe scales

def bin_weighted_average(x, y, bins):
    Nbins = len(bins)
    out = np.zeros(Nbins-1)
    bidxs = np.digitize(x, bins)
    for i in xrange(Nbins-1):
        inbin = bidxs == i
        out[i] = trapz(y[inbin], x[inbin])/(bins[i+1] - bins[i])
    return out

hpids = np.arange(10, 52, 2)
dbin = 200.
bins = np.arange(-1500., -100. + dbin, dbin)
eps_av = np.zeros((len(bins) - 1, len(hpids)))
z_av = np.zeros((len(bins) - 1, len(hpids)))
d_av = np.zeros((len(bins) - 1, len(hpids)))

fig = plt.figure(figsize=(3.125, 3))
gs = gridspec.GridSpec(2, 1, height_ratios=[1, 5])
ax0 = plt.subplot(gs[1])
ax1 = plt.subplot(gs[0])

for Float in [E76, E77]:

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)

    N2_ref = Float.N2_ref[:, idxs]
    zm = getattr(Float, 'z')[:, idxs]
    zw = getattr(Float, 'zw')[:, idxs]
    d = getattr(Float, 'dist_ctd')[:, idxs]

    ts, td, Nsq = mdp.thorpe_float(Float, hpids, zvar='zw', R0=0.3, acc=1.6e-3,
                                   use_int=True)

    eps_thorpe = (0.8*ts)**2 * Nsq**(3./2.)

#    LOG_EPS = np.log10(eps_thorpe).flatten(order='F')
#    Z = z.flatten(order='F')
#    X = Float.dist_ctd[:, idxs].flatten(order='F')

    ieps = 0.*np.zeros_like(idxs)

    for i in xrange(len(idxs)):
        eps_av[:, i] = bin_weighted_average(zm[:, i], eps_thorpe[:, i], bins)
        z_av[:, i], __, __ = binned_statistic(zm[:, i], zm[:, i], statistic=np.nanmean, bins=bins)
        d_av[:, i], __, __ = binned_statistic(zm[:, i], d[:, i], statistic=np.nanmean, bins=bins)

        use = zm[:, i] < -100.

        ieps[i] = np.abs(1025.*trapz(eps_thorpe[use, i], zw[use, i]))


    # Plotting #
    # Epsilon
    lon = getattr(Float, 'lon_start')[idxs]
    lat = getattr(Float, 'lat_start')[idxs]
    bathy = sandwell.interp_track(lon, lat, bf)
    dbathymax = Float.dist[idxs][bathy.argmax()]

    pdist = d_av - dbathymax

    ax1.plot(pdist[0, :], 1000.*ieps, label=Float.floatID)

    LOG_EPS_AV = np.ma.masked_invalid(np.log10(eps_av))

#    LOG_EPS[~np.isfinite(LOG_EPS)] = np.NaN

    sc = ax0.scatter(pdist.flatten(), z_av.flatten(), s=80., c=LOG_EPS_AV,
                     cmap=plt.get_cmap('YlOrRd'), vmin=-10., vmax=-7, alpha=.5)

#    step = 1
#    sc = ax0.scatter(X[::step], Z[::step], s=5, c=LOG_EPS[::step],
#                     edgecolor='none', cmap=plt.get_cmap('YlOrRd'), vmin=-10.,
#                     vmax=-7, alpha=.5)

ax1.set_ylabel('$P$ (mW m$^{-2}$)')
#ax1.yaxis.set_ticks(np.array([0., 5., 10.]))
ax1.xaxis.set_ticks([])

ax1.legend(loc='upper right', fontsize=7)

ax0.fill_between(pdist[0, :], bathy, np.nanmin(bathy), color='black', linewidth=2)
ax0.set_ylim(-4000., 0.)
ax0.yaxis.set_ticks(np.arange(-4000, 1000, 1000))
ax0.yaxis.set_ticklabels(['-4', '-3', '-2', '-1', '0'])

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.7])
C = fig.colorbar(sc, cax=cbar_ax, extend='both')
C.set_label(r'$\log_{10}(\epsilon)$ (W kg$^{-1}$)')

#plt.clim(-11., -7.)
ax0.set_xlim(np.min(pdist), np.max(pdist))

ax0.set_xlabel('Distance from ridge top (km)')
ax0.set_ylabel('$z$ (km)')

ax1.set_xlim(*ax0.get_xlim())

fontdict={'size': 10}
plt.figtext(-0.05, 0.85, 'a)', fontdict=fontdict)
plt.figtext(-0.05, 0.65, 'b)', fontdict=fontdict)

pf.my_savefig(fig, 'both', 'epsilon_thorpe', sdir, ftype=('png', 'pdf'),
              fsize='single_col')


# %% Using finescale parameterisation

params = TKED.default_params

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

# hpids is set further up the script.

for Float in [E76, E77]:

    results = mdp.analyse_float(Float, hpids, params)
    __, idxs = Float.get_profiles(hpids, ret_idxs=True)
    dists = Float.dist[idxs]
    bathy = sandwell.interp_track(Float.lon_start, Float.lat_start, bf)

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
    ieps = 0.*np.zeros_like(idxs)

    iZ = z_meang[:, 0]

    for i in xrange(len(idxs)):
        ieps[i] = 1025.*trapz(epsilong[::-1, i], iZ[::-1])


    Z = z_meang.flatten()
    X = np.asarray(len(z_meang[:, 0])*[Float.dist[idxs]])
    LOG_EPS = (np.log10(epsilong)).flatten()
    LOG_KAP = (np.log10(kappag)).flatten()

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

    pf.my_savefig(fig, Float.floatID, 'epsilon_fs', sdir, fsize='double_col')

#    ylims = (params['zmin'], params['zmax'])
#
#    fig = plt.figure()
#    plt.title('log10 R_pol')
#    plt.pcolormesh(dists, z_meang[:,0], np.log10(R_polg), cmap=plt.get_cmap('bwr'))
#    plt.clim(-1, 1)
#    plt.colorbar()
#    plt.ylim(ylims)
#    plt.xlim(np.min(dists), np.max(dists))
#    pf.my_savefig(fig, Float.floatID, 'R_pol_fs', sdir)
#
#    fig = plt.figure()
#    plt.title('log10 epsilon')
#    plt.pcolormesh(dists, z_meang[:,0], np.log10(epsilong), cmap=plt.get_cmap('bwr'))
#    plt.clim(-11, -7)
#    plt.colorbar()
#    plt.ylim(ylims)
#    plt.xlim(np.min(dists), np.max(dists))
#    pf.my_savefig(fig, Float.floatID, 'epsilon_fs', sdir)

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

# %% Combined Thorpe and LEM

def bin_weighted_average(x, y, bins):
    Nbins = len(bins)
    out = np.zeros(Nbins-1)
    bidxs = np.digitize(x, bins)
    for i in xrange(Nbins-1):
        inbin = bidxs == i
        out[i] = trapz(y[inbin], x[inbin])/(bins[i+1] - bins[i])
    return out

# Using large eddy method first.
# Different coefficient for each float.
#cs = [0.197, 0.158]  # time
#xvar = 'time'
#dx = 5.
#x = np.arange(0., 12000., dx)
#width = 120.
#lc = np.array([300., 120.])
#btype = 'bandpass'

#cs = [0.193, 0.160]  # height
#xvar = 'height'
#dx = 1.
#x = np.arange(-1500., 0., dx)
#width = 15.
#lc = np.array([40., 15.])
#btype = 'bandpass'

#cs = [0.176, 0.147] # timeheight
#xvar = 'timeheight'
#dx = 1.
#x = np.arange(-1450., -50, dx)
#width = 15.
#lc = np.array([100., 40.])
#btype = 'highpass'

#cs = [0.192, 0.159]  # eheight
#x = np.arange(zmin, 0., dz)
#xvar = 'eheight'
#dx = 1.
#width = 20.
#lc = np.array([40., 15.])
#btype = 'bandpass'

cs = [0.146, 0.123]  # timeeheight
x = np.arange(zmin, 0., dz)
xvar = 'timeeheight'
dx = 1.
hpids = np.arange(10, 38)  # Must start even.
width = 20.
lc = (100., 40.)
btype = 'highpass'
we = 0.001

hpids_t = np.arange(hpids[0], hpids[-1], 2)
dbin = 200.
bins = np.arange(-1500., -100. + dbin, dbin)
eps_av = np.zeros((len(bins) - 1, len(hpids_t)))
z_av = np.zeros((len(bins) - 1, len(hpids_t)))
d_av = np.zeros((len(bins) - 1, len(hpids_t)))

fig = plt.figure(figsize=(3.125, 4))
gs = gridspec.GridSpec(2, 1, height_ratios=[3, 5])
#gs.update(wspace=0.1)
ax0 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])
#ax2 = plt.subplot(gs[2])

colors = ['b', 'g']

for j, (Float, c) in enumerate(zip([E76, E77], cs)):

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)

    eps_lem, __, __, noise = mdp.w_scales_float(Float, hpids, xvar, x,
                                                width=width, overlap=-1.,
                                                lc=lc, c=c, btype=btype,
                                                we=we, ret_noise=True)

    __, idxs_t = Float.get_profiles(hpids_t, ret_idxs=True)

    N2_ref = Float.N2_ref[:, idxs_t]
    zm = getattr(Float, 'z')[:, idxs_t]
    zw = getattr(Float, 'zw')[:, idxs_t]
    d = getattr(Float, 'dist_ctd')[:, idxs_t]

    ts, td, Nsq = mdp.thorpe_float(Float, hpids_t, zvar='zw', R0=0.3, acc=1.6e-3,
                                   use_int=True)

    eps_thorpe = (0.8*ts)**2 * Nsq**(3./2.)

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

    ieps_lem = 0.*np.zeros_like(idxs)

    for i in xrange(len(idxs)):
        iuse = (iZ[:, i] < -100) & (iZ[:, i] > -1400)
        # The abs function accounts for problems with z being the wrong way.
        ieps_lem[i] = np.abs(1025.*trapz(eps_lem[iuse, i], iZ[iuse, i]))

    ieps_t = 0.*np.zeros_like(idxs_t)

    for i in xrange(len(idxs_t)):
        eps_av[:, i] = bin_weighted_average(zm[:, i], eps_thorpe[:, i], bins)
        z_av[:, i], __, __ = binned_statistic(zm[:, i], zm[:, i], statistic=np.nanmean, bins=bins)
        d_av[:, i], __, __ = binned_statistic(zm[:, i], d[:, i], statistic=np.nanmean, bins=bins)

        use = zm[:, i] < -100.

        ieps_t[i] = np.abs(1025.*trapz(eps_thorpe[use, i], zw[use, i]))


    Z = iZ.flatten(order='F')

    use = (Z < -10) & (Z > -1400)

    Z = Z[use]

    X = X.flatten(order='F')[use]
    noise = noise.flatten(order='F')[use]
    LOG_EPS = (np.log10(eps_lem)).flatten(order='F')[use]

    # Plotting #
    # Epsilon
    d_lem = getattr(Float, 'dist_ctd')[:, idxs].flatten(order='F')

    tgps = getattr(Float, 'UTC_start')[idxs]
    lon = getattr(Float, 'lon_start')[idxs]
    lat = getattr(Float, 'lat_start')[idxs]
    tctd = getattr(Float, 'UTC')[:, idxs].flatten(order='F')
    nans = np.isnan(d_lem) | np.isnan(tctd)
    tctd = tctd[~nans]
    dctd = d_lem[~nans]
    lonctd = np.interp(tctd, tgps, lon)
    latctd = np.interp(tctd, tgps, lat)
    bathy = sandwell.interp_track(lonctd, latctd, bf)

    dbathymax = dctd[bathy.argmax()]

    dctd -= dbathymax
    X -= dbathymax

    pdist = d_av - dbathymax

    ax0.plot(pdist[0, :], 1000.*ieps_t, color=colors[j], linestyle=':', label="{} Thorpe".format(Float.floatID))

    LOG_EPS_AV = np.ma.masked_invalid(np.log10(eps_av))

    sc = ax2.scatter(pdist.flatten(), z_av.flatten(), s=80., c=LOG_EPS_AV.flatten(),
                     cmap=plt.get_cmap('YlOrRd'), vmin=-10., vmax=-7, alpha=.5)

    ax0.plot(Float.dist[idxs] - dbathymax, 1000.*ieps_lem, color=colors[j], linestyle='-', label="{} LEM".format(Float.floatID))

    LOG_EPS[noise] = np.NaN

    step = 20
    sc = ax2.scatter(X[::step], Z[::step], s=5, c=LOG_EPS[::step],
                     edgecolor='none', cmap=plt.get_cmap('YlOrRd'), vmin=-10.,
                     vmax=-7, alpha=.5)


ax0.set_ylabel('$P$ (mW m$^{-2}$)')
ax0.yaxis.set_ticks(np.array([0., 12.5, 25.]))
ax0.xaxis.set_ticks([])
ax0.legend(loc='upper left', fontsize=7, ncol=1, bbox_to_anchor=(0.8, 1.02), borderaxespad=0.)

#ax1.xaxis.set_ticks([])
#ax1.set_ylim(-1500., 0.)
#ax1.yaxis.set_ticks(np.arange(-1500, 500, 500))
#ax1.yaxis.set_ticklabels(['-1.5', '-1', '-0.5', '0'])

ax2.fill_between(dctd[::100], bathy[::100],
                 np.nanmin(bathy), color='black', linewidth=2)
ax2.set_ylim(-4000., 0.)
ax2.yaxis.set_ticks(np.arange(-4000, 1000, 1000))
ax2.yaxis.set_ticklabels(['-4', '-3', '-2', '-1', '0'])

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.4])
C = fig.colorbar(sc, cax=cbar_ax, extend='both')
C.set_label(r'$\log_{10}(\epsilon)$ (W kg$^{-1}$)')

#    plt.clim(-11., -7.)
ax2.set_xlim(np.min(X), np.max(X))

ax2.set_xlabel('Distance from ridge top (km)')
#ax1.set_ylabel('$z$ (km)')
ax2.set_ylabel('$z$ (km)')

ax0.set_xlim(*ax2.get_xlim())
#ax1.set_xlim(*ax2.get_xlim())

fontdict={'size': 10}
plt.figtext(-0.05, 0.88, 'a)', fontdict=fontdict)
plt.figtext(-0.05, 0.53, 'b)', fontdict=fontdict)

pf.my_savefig(fig, 'both', 'epsilon_lem_thorpe_combined', sdir, ftype=('png', 'pdf'),
              fsize='single_col')



