"""
Plot a LOT of things.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.optimize as op
import gsw
import os
import sys
from scipy.interpolate import griddata
from scipy.integrate import cumtrapz
import scipy.signal as sig

lib_path = os.path.abspath('../modules')
if lib_path not in sys.path:
    sys.path.append(lib_path)

import sandwell
import emapex
import plotting_functions as pf


try:
    print("Floats {} and {} exist!.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.load(4976)
    E77 = emapex.load(4977)

plot_dir = '../figures/large_wave_analysis'
ftype = 'png'
bwr = plt.get_cmap('bwr')
bwr1 = mcolors.LinearSegmentedColormap.from_list(name='red_white_blue',
                                                 colors=[(0, 0, 1),
                                                         (1, 1, 1),
                                                         (1, 0, 0)],
                                                 N=40)


def my_savefig(fid, fname):
    fname = str(fid) + '_' + fname
    fname = os.path.join(plot_dir, fname) + '.' + ftype
    plt.savefig(fname, bbox_inches='tight')

# %% Sections of temperature, salinity and potential density.

hpids = np.arange(1, 600)

vars = ['T', 'S']  # , 'rho_1']
texvars = ['$T$', '$S$']  # , r'$\sigma_1$']
units = [r'$^\circ$C', 'PSU']  # , 'kg m$^{-3}$']
clims = [(2., 6.), (34., 34.6)]  # , (31.2, 32.3)]
cmaps = [plt.get_cmap(s) for s in ['rainbow', 'summer', 'spring']]

for Float in [E76, E77]:

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)
    z = getattr(Float, 'z')[:, idxs].flatten(order='F')
    d = getattr(Float, 'dist_ctd')[:, idxs].flatten(order='F')

    Xg, Zg, rho_1g = Float.get_griddata_grid(hpids, 'dist_ctd', 'z', 'rho_1')
    rho_1g -= 1000.

    for var, texvar, unit, clim, cmap in zip(vars, texvars, units, clims,
                                             cmaps):

        C = getattr(Float, var)[:, idxs].flatten(order='F')

        nans = np.isnan(z) | np.isnan(d) | np.isnan(C)

        plt.figure(figsize=(6, 4))
        plt.scatter(d[~nans], z[~nans], c=C[~nans], edgecolor='none',
                    cmap=cmap)
        cbar = plt.colorbar(orientation='horizontal', extend='both')
        cbar.set_label(texvar+' ('+unit+')')
        plt.clim(*clim)

        plt.contour(*Float.get_griddata_grid(hpids, 'dist_ctd', 'z', var),
                    N=6, colors='k')
        CS = plt.contour(Xg, Zg, rho_1g, N=6, colors='w')
#        plt.clabel(CS, inline=1, fontsize=8)

        plt.ylim(np.min(z[~nans]), np.max(z[~nans]))
        plt.xlim(np.min(d[~nans]), np.max(d[~nans]))
        plt.xlabel('Distance (km)')
        plt.ylabel('Depth (m)')
        title_str = ("Float {}").format(Float.floatID)
        plt.title(title_str)

        my_savefig(Float.floatID, var)

# %% Sections of horizontal velocity.

hpids = np.arange(1, 600)

for Float in [E76, E77]:

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)
    z = getattr(Float, 'zef')[:, idxs].flatten(order='F')
    d = getattr(Float, 'dist_ef')[:, idxs].flatten(order='F')

    for comp in ['U_abs', 'V_abs']:

        V = getattr(Float, comp)[:, idxs].flatten(order='F')

        plt.figure(figsize=(6, 4))
        plt.scatter(d, z, c=V, edgecolor='none', cmap=bwr)
        cbar = plt.colorbar(orientation='horizontal', extend='both')
        cbar.set_label('${}$ (m s$^{{-1}}$)'.format(comp[0]))
        plt.clim(-1., 1.)

        Xg, Zg, Vg = Float.get_griddata_grid(hpids, 'dist_ef', 'zef', comp)
        Vg = np.abs(Vg)
        plt.contour(Xg, Zg, Vg, levels=[0.5], colors='k',
                    linestyles='dashed')

        plt.ylim(np.nanmin(z), np.nanmax(z))
        plt.xlim(np.nanmin(Float.dist_ctd), np.nanmax(Float.dist_ctd))
        plt.xlabel('Distance (km)')
        plt.ylabel('Depth (m)')
        title_str = ("Float {}").format(Float.floatID)
        plt.title(title_str)

        my_savefig(Float.floatID, comp)

# %% A single profile showing Ws, Wf and Ww.

hpid = 88
Float = E77

pfl = Float.get_profiles(hpid)
plt.figure(figsize=(4, 6))
plt.plot(pfl.Ww*100., pfl.z, 'k')
plt.plot(pfl.Wz*100., pfl.z, color='grey')
plt.plot(pfl.Ws*100., pfl.z, 'r')

plt.xlabel('$W_w$, $W_f$, $W_f$ (cm s$^{-1}$)')
plt.xticks(rotation=35)
plt.ylabel('Depth (m)')
title_str = ("Float {}, half profile {}").format(Float.floatID, hpid)
plt.title(title_str)
plt.legend(['$W_w$', '$W_f$', '$W_f$'], loc=0)

my_savefig(Float.floatID, 'W_example')

# %% Sections of Ww.

hpids = np.arange(1, 600)

for Float in [E76, E77]:

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)
    z = getattr(Float, 'z')[:, idxs].flatten(order='F')
    d = getattr(Float, 'dist_ctd')[:, idxs].flatten(order='F')
    Ww = getattr(Float, 'Ww')[:, idxs].flatten(order='F')

    plt.figure(figsize=(12, 4))
    plt.scatter(d, z, c=Ww*100., edgecolor='none', cmap=bwr)

    cbar = plt.colorbar(orientation='horizontal', extend='both')
    cbar.set_label('$W_w$ (cm s$^{-1}$)')
    plt.clim(-5, 5)

    Xg, Zg, Wg = Float.get_griddata_grid(hpids, 'dist_ctd', 'z', 'Ww')
    Wg = np.abs(Wg)
    plt.contour(Xg, Zg, Wg, levels=[0.025], colors='k',
                linestyles='dashed')

    plt.ylim(np.nanmin(z), np.nanmax(z))
#    plt.xlim(np.nanmin(d), np.nanmax(d))
    plt.xlabel('Distance (km)')
    plt.ylabel('Depth (m)')
    plt.xlim(np.nanmin(d), np.nanmax(d))
    title_str = ("Float {}").format(Float.floatID)
    plt.title(title_str)

    my_savefig(Float.floatID, 'Ww_section')

# %% Zoomed in variables near the wave.

hpids = np.arange(10, 50)
bathy_file = '/noc/users/jc3e13/storage/smith_sandwell/topo_17.1.img'
vars = ['Ww', 'U_abs', 'V_abs']
zvars = ['z', 'zef', 'zef']
dvars = ['dist_ctd', 'dist_ef', 'dist_ef']
texvars = ['$W_w$', '$U$', '$V$']
clims = [(-10., 10.), (-100., 100.), (-100, 100.)]

for Float in [E76, E77]:

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)

    for var, zvar, dvar, texvar, clim in zip(vars, zvars, dvars, texvars,
                                             clims):

        V = getattr(Float, var)[:, idxs].flatten(order='F')
        z = getattr(Float, zvar)[:, idxs].flatten(order='F')
        d = getattr(Float, dvar)[:, idxs].flatten(order='F')

        plt.figure(figsize=(5, 7))
        plt.scatter(d, z, s=50, c=V*100., edgecolor='none', cmap=bwr)
        cbar = plt.colorbar(orientation='horizontal', extend='both')
        cbar.set_label(texvar+' (cm s$^{-1}$)')
        plt.clim(*clim)

        plt.xlim(np.nanmin(d), np.nanmax(d))
        plt.xlabel('Distance (km)')
        plt.ylabel('Depth (m)')
        title_str = ("Float {}").format(Float.floatID)
        plt.title(title_str)

        lons = Float.lon_start[idxs]
        lats = Float.lat_start[idxs]
        dist = Float.dist[idxs]
        bathy = sandwell.interp_track(lons, lats, bathy_file)
        plt.plot(dist, bathy, 'k', linewidth=3)

        plt.ylim(np.nanmin(bathy), np.nanmax(z))

        my_savefig(Float.floatID, var + '_mountain')

    pf.track_on_bathy(Float, hpids, bathy_file=bathy_file)
    my_savefig(Float.floatID, 'mountain_area_track')


# %% Fitting in only one direction. Vertical wavenumber.

def plane_wave(x, A, k, phase):
    return A*np.cos((k*x + phase))

E76_hpids = np.arange(31, 33)
E77_hpids = np.arange(23, 24)

for Float, hpids in zip([E76, E77], [E76_hpids, E77_hpids]):

    for pfl in Float.get_profiles(hpids):

        z = getattr(pfl, 'z')
        Ww = getattr(pfl, 'Ww')

        # Remove nans and the top 50 m.
        nans = np.isnan(z) | np.isnan(Ww) | (z > -50)
        z, Ww = z[~nans], Ww[~nans]

        # Setting up the spectral analysis.
        dz = np.abs(np.round(np.mean(np.diff(z))))
        dk = 1./dz
        iz = np.arange(np.min(z), np.max(z), dz)
        idxs = np.argsort(z)
        iWw = np.interp(iz, z[idxs], Ww[idxs])
        diWw = sig.detrend(iWw)

        # Doing the spectral analysis.
        m, Pw = sig.welch(iWw, fs=dk, detrend='linear')

        # Try fitting plane wave.
        popt, __ = op.curve_fit(plane_wave, iz, iWw,
                                p0=[0.15, 0.015, 0.])
        mfit = popt[1]/(2.*np.pi)

        plt.figure(figsize=(4,6))

#        plt.subplot(1, 2, 1)
        plt.plot(iWw*100., iz, plane_wave(iz, *popt)*100., iz)
        plt.xticks(rotation=45)
        plt.xlabel('$W_w$ (cm s$^{-1}$)')
        plt.ylabel('$z$ (m)')
        title = ("Float {}, profile {}\nm_fit {:1.1e} m$^{{-1}}$"
                 "\nlambda_fit {:4.0f} m"
                 "").format(Float.floatID, pfl.hpid[0], mfit, 1./mfit)
        plt.title(title)
        plt.legend(['$W_w$', 'fit'])
        my_savefig(Float.floatID, str(pfl.hpid) + '_mfit')
#        plt.subplot(1, 2, 2)
#        plt.loglog(m, Pw)
#        mmax = m[Pw.argmax()]
#        title = ("Float {}, profile {}\nm_max {:1.1e} m$^{{-1}}$"
#                 "\nlambda_max {:4.0f} m"
#                 "").format(Float.floatID, pfl.hpid, mmax, 1./mmax)
#        ylim = plt.ylim()
#        plt.plot(2*[mmax], ylim, 'r', 2*[mfit], ylim, 'g')
#        plt.title(title)
#        plt.xlabel('$m$ (m$^{-1}$)')

# %% Fitting in time only.

def plane_wave(x, A, k, phase):
    return A*np.cos((k*x + phase))

E76_hpids = np.arange(31, 33)
E77_hpids = np.arange(23, 24)

for Float, hpids in zip([E76, E77], [E76_hpids, E77_hpids]):

    for pfl in Float.get_profiles(hpids):

        z = getattr(pfl, 'z')
        t = getattr(pfl, 'UTC')
        t = 24.*60.*(t - np.nanmin(t))  # Convert to minutes.

        Ww = getattr(pfl, 'Ww')
        N2_ref = getattr(pfl, 'N2_ref')

        nans = np.isnan(z) | np.isnan(t) | np.isnan(Ww) | np.isnan(N2_ref) | (z > -50)
        z, t, Ww, N2_ref = z[~nans], t[~nans], Ww[~nans], N2_ref[~nans]

        # Setting up the spectral analysis.
        ts = 60.*t  # Convert to seconds.
        dt = np.round(np.mean(np.diff(ts)))
        fs = 1./dt
        it = np.arange(np.min(ts), np.max(ts), dt)
        iWw = np.interp(it, ts, Ww)
        diWw = sig.detrend(iWw)

        # Get an idea of the buoyancy frequency.
        N_mean = np.mean(np.sqrt(N2_ref[z < -200]))/(2.*np.pi)
        # The inertial frequency.
        fcor = np.abs(gsw.f(-57.5))/(2.*np.pi)

        # Perform the spectral analysis.
        freqs, Pw = sig.welch(iWw, fs=fs, detrend='linear')

        # Fit a curve.
        popt, __ = op.curve_fit(plane_wave, it, diWw,
                                p0=[0.1, 0.002, 0.])
        omfit = popt[1]/(2.*np.pi)
        period = 1/omfit/60.

        plt.figure(figsize=(4,6))

#        plt.subplot(1, 2, 1)
        plt.plot(diWw, it, plane_wave(it, *popt), it)
        plt.xticks(rotation=45)
        plt.xlabel('$W_w$ (m s$^{-1}$)')
        plt.ylabel('$t$ (s)')
        title = ("Float {}, profile {}\nperiod {} min").format(Float.floatID, pfl.hpid[0], period)
        plt.title(title)
        plt.legend(['$W_w$', 'fit'])
        my_savefig(Float.floatID, str(pfl.hpid) + '_ofit')

#        plt.subplot(1, 2, 2)
#        plt.loglog(freqs, Pw)
#        ylim = plt.ylim()
#        plt.plot(2*[N_mean], ylim, 'r', 2*[fcor], ylim, 'r',  2*[omfit], ylim, 'g')
#        plt.title(title)
#        plt.xlabel('$f$ (s$^{-1}$)')
#        title = ("Float {}, profile {}").format(Float.floatID, pfl.hpid[0])

# %% A few waves

hpids = np.arange(30, 33)
Float = E76
dttype = 'linear'
xticks = [-20, 0, 20]

for hpid in hpids:
    pfl = Float.get_profiles(hpid)

    uvnans = np.isnan(pfl.U) | np.isnan(pfl.V)
    wnans = np.isnan(pfl.Ww)

    fig, ax = plt.subplots(1, 3, sharey=True, figsize=(4, 6))
    title_str = ("Float {}, profile {}").format(Float.floatID, hpid)
    plt.suptitle(title_str)

    ax[0].plot(pfl.Ww[~wnans]*100., pfl.z[~wnans], 'b')
    ax[0].plot(2*[0.], [-1600, 0], 'k-')
    ax[0].set_ylabel('Depth (m)')
    ax[0].set_xlabel('$W_w$ (cm s$^{-1}$)')
    ax[0].grid()
    ax[0].set_xticks(ax[0].get_xlim())
    ax[0].set_xticklabels(ax[0].get_xticks(), rotation=45)

    ax[1].plot(pfl.U[~uvnans]*100., pfl.zef[~uvnans], 'b')
    ax[1].plot(2*[0.], [-1600, 0], 'k-')
    ax[1].set_xlabel('$U$ (cm s$^{-1}$)')
    ax[1].grid()
    ax[1].set_xticks(ax[1].get_xlim())
    ax[1].set_xticklabels(ax[1].get_xticks(), rotation=45)

    ax[2].plot(pfl.V[~uvnans]*100., pfl.zef[~uvnans], 'b')
    ax[2].plot(2*[0.], [-1600, 0], 'k-')
    ax[2].set_xlabel('$V$ (cm s$^{-1}$)')
    ax[2].grid()
    ax[2].set_xticks(ax[2].get_xlim())
    ax[2].set_xticklabels(ax[2].get_xticks(), rotation=45)

    my_savefig(Float.floatID, 'WU_profile_' + str(hpid))

# %% Fitting in all dimensions.

def plane_wave2(params, x):
    A, k, m, om, phi = params
    return A*np.cos(k*x[:, 0] + m*x[:, 1] + om*x[:, 2] + phi)


def wave_3(params, x):
    A, k, m, om, phi= params
    x[:, 0] = 0.
    z = cumtrapz(x[:, 3], x[:, 2], initial=0.) - \
        A/om*np.sin(k*x[:, 0] + m*x[:, 1] - om*x[:, 2] + phi)
    return A*np.cos(k*x[:, 0] + m*z - om*x[:, 2] + phi)


def cost(params, data, func, y):
    return (func(params, data) - y).flatten()

res = []

E76_hpids = np.arange(31, 32) # np.arange(31, 33)
E77_hpids = np.arange(26, 27) # np.arange(26, 28)

for Float, hpids in zip([E76, E77], [E76_hpids, E77_hpids]):

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)

    z = Float.z[:,idxs].flatten(order='F')
    x = Float.dist_ctd[:, idxs].flatten(order='F')*1000.
    t = Float.UTC[:, idxs].flatten(order='F')*86400.
    t -= np.nanmin(t)
    W = Float.Ww[:, idxs].flatten(order='F')
    Ws = Float.Ws[:, idxs].flatten(order='F')

    nans = np.isnan(z) | np.isnan(x) | np.isnan(t) | np.isnan(W) | (z > -600)

    data = np.array([x[~nans], z[~nans], t[~nans], Ws[~nans]]).T
    W2 = W[~nans]

    x0 = [0.15, 0.004, 0.006, 0.0004, 0.]

    fit = op.leastsq(cost, x0=x0, args=(data, wave_3, W2))[0]
    print(fit)
    freq = fit[3]/(2.*np.pi)
    period = 1./freq/60.
    wavenum = fit[2]/(2.*np.pi)
    wavelen = 1./wavenum

    print("Frequency = {} s-1, Period = {} min.".format(freq, period))
    print("Wavenumber = {} m-1, Wavelength = {} m.".format(wavenum, wavelen))
    res.append(fit)

    Wm = wave_3(fit, data)
    Wm2 = wave_3(fit, np.array([x, z, t, Ws]).T)
    Wm0 = wave_3(x0, data)

    plt.figure()
    plt.subplot(2, 1, 1)
    plt.plot(Wm*100., data[:,1], W*100., z)#, Wm0, data[:,1])
    plt.xlabel('$W_w$ (cm s$^{-1}$)')
    plt.ylabel('Depth (m)')
    plt.subplot(2, 1, 2)
    plt.plot(data[:,2], Wm*100., t, W*100.)#, data[:,2], Wm0)
    plt.xlabel('Time (s)')
    plt.ylabel('$W_w$ (cm s$^{-1}$)')

    my_savefig(Float.floatID, 'wave_3_fit')
# hpid 31 and 26
#[ -1.30611467e-01   4.00000000e-03   5.37679530e-03  -5.51902074e-04
#   4.99192551e+00  -1.07282887e+00   3.44287143e+00]
#Frequency = -8.78379431994e-05 s-1, Period = -189.743362146 min.
#Wavenumber = 0.000855743550121 m-1, Wavelength = 1168.57439341 m.
#[  1.56557333e-01   4.00000000e-03   6.24034100e-03   1.49703041e-03
#  -2.50568019e+00   6.85488386e-01  -8.23274676e+00]
#Frequency = 0.000238259789373 s-1, Period = 69.951655336 min.
#Wavenumber = 0.000993181117242 m-1, Wavelength = 1006.86569916 m.


# %%

bwr = plt.get_cmap('bwr')
E76_hpids = np.arange(20, 40) # np.arange(31, 33)
E77_hpids = np.arange(15, 35) # np.arange(26, 28)
bathy_file = '../../data/sandwell_bathymetry/topo_17.1.img'
vars = ['Ww']#, 'U_abs', 'V_abs']
zvars = ['z']#, 'zef', 'zef']
dvars = ['dist_ctd']#, 'dist_ef', 'dist_ef']
texvars = ['$W_w$']#, '$U$', '$V$']
clims = [(-10., 10.)]#, (-100., 100.), (-100, 100.)]

var_1_vals = np.linspace(-40., 40., 80)
var_2_vals = np.linspace(-1500, 0, 500)
Xg, Zg = np.meshgrid(var_1_vals, var_2_vals)

Wgs = []
ds = []
zs = []

for Float, hpids in zip([E76, E77], [E76_hpids, E77_hpids]):

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)

    for var, zvar, dvar, texvar, clim in zip(vars, zvars, dvars, texvars,
                                             clims):

        V = getattr(Float, var)[:, idxs].flatten(order='F')
        z = getattr(Float, zvar)[:, idxs].flatten(order='F')
        d = getattr(Float, dvar)[:, idxs].flatten(order='F')
        zs.append(z.copy())

        tgps = getattr(Float, 'UTC_start')[idxs]
        lon = getattr(Float, 'lon_start')[idxs]
        lat = getattr(Float, 'lat_start')[idxs]
        tctd = getattr(Float, 'UTC')[:, idxs].flatten(order='F')
        nans = np.isnan(d) | np.isnan(tctd)
        tctd = tctd[~nans]
        dctd = d[~nans]
        lonctd = np.interp(tctd, tgps, lon)
        latctd = np.interp(tctd, tgps, lat)
        bathy = sandwell.interp_track(lonctd, latctd, bathy_file)

        d -= dctd[bathy.argmax()]
        ds.append(d.copy())

        nans = np.isnan(d) | np.isnan(z) | np.isnan(V)

        Wg = griddata((d[~nans], z[~nans]), V[~nans], (Xg, Zg), method='linear')
        Wgs.append(Wg.copy())

        dctd -= dctd[bathy.argmax()]
        # Spectral analysis of bathymetry.
        dx = np.mean(np.diff(dctd))
        dk = 1./dx
        ix = np.arange(-40., 40., dx)
        ibathy = np.interp(ix, dctd, bathy)
        k, Pbathy = sig.welch(ibathy, fs=dk/1000., detrend='linear')

        plt.figure()
        mWg = np.ma.masked_where(np.isnan(Wg), Wg)
        plt.pcolormesh(Xg, Zg, mWg*100., cmap=bwr)
        plt.plot(dctd, bathy, 'k', linewidth=2)
        cbar = plt.colorbar(orientation='horizontal', extend='both')
        cbar.set_label(texvar+' (cm s$^{-1}$)')
        plt.clim(*clim)
        plt.scatter(d, z, s=1, edgecolor='none', color='grey')

        plt.ylabel('Depth (m)')
        title_str = ("Float {}").format(Float.floatID)
        plt.title(title_str)
        plt.xlim(-10., np.nanmax(d))
        plt.ylim(np.nanmin(bathy), np.nanmax(z))
        plt.grid()

        plt.figure()
        plt.scatter(d, z, s=50, c=V*100., edgecolor='none', cmap=bwr)
        cbar = plt.colorbar(orientation='horizontal', extend='both')
        cbar.set_label(texvar+' (cm s$^{-1}$)')
        plt.clim(*clim)

        plt.xlim(-10., np.nanmax(d))
        plt.xlabel('Distance (km)')
        plt.ylabel('Depth (m)')
        title_str = ("Float {}").format(Float.floatID)
        plt.title(title_str)

        plt.plot(dctd, bathy, 'k', linewidth=2)

        plt.ylim(np.nanmin(bathy), np.nanmax(z))

        plt.grid()

        my_savefig(Float.floatID, 'mountain_closeup')

        plt.figure(figsize=(5, 7))
        plt.loglog(k, Pbathy)



Wg_diff = Wgs[0] - Wgs[1]
mWg_diff = np.ma.masked_where(np.isnan(Wg_diff), Wg_diff)

plt.figure()
plt.pcolormesh(Xg, Zg, Wg_diff*100., cmap=bwr)
plt.plot(dctd, bathy, 'k', linewidth=2)
cbar = plt.colorbar(orientation='horizontal', extend='both')
cbar.set_label('$W_w$ difference (cm s$^{-1}$)')
plt.clim(*clim)
plt.scatter(ds[0], zs[0], s=1, edgecolor='none', color='grey')
plt.scatter(ds[0], zs[0], s=1, edgecolor='none', color='grey')

plt.ylabel('Depth (m)')
plt.xlim(np.nanmin(d), np.nanmax(d))
plt.ylim(np.nanmin(bathy), np.nanmax(z))