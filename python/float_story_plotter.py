"""
Plot a LOT of things.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import plotting_functions as pf
import emapex
import scipy.optimize as op
import pylab as pyl
import os
import sandwell

reload(emapex)

E76 = emapex.EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4976)
E76.apply_w_model('../../data/EM-APEX/4976_fix_p0k0M_fit_info.p')
E76.apply_strain('../../data/EM-APEX/4976_N2_ref_300dbar.p')

E77 = emapex.EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4977)
E77.apply_w_model('../../data/EM-APEX/4977_fix_p0k0M_fit_info.p')
E77.apply_strain('../../data/EM-APEX/4977_N2_ref_300dbar.p')

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

vars = ['T', 'S', 'rho_1']
texvars = ['$T$', '$S$', r'$\sigma_1$']
units = [r'$^\circ$C', 'PSU', 'kg m$^{-3}$']
clims = [(2., 6.), (34., 34.6), (31.2, 32.3)]
cmaps = [plt.get_cmap(s) for s in ['rainbow', 'summer', 'spring']]

for Float in [E76, E77]:

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)
    z = getattr(Float, 'z')[:, idxs].flatten(order='F')
    d = getattr(Float, 'dist_ctd')[:, idxs].flatten(order='F')

    for var, texvar, unit, clim, cmap in zip(vars, texvars, units, clims,
                                             cmaps):

        C = getattr(Float, var)[:, idxs].flatten(order='F')

        if var == 'rho_1':
            C -= 1000.

        nans = np.isnan(z) | np.isnan(d) | np.isnan(C)

        plt.figure(figsize=(6, 4))
        plt.scatter(d[~nans], z[~nans], c=C[~nans], edgecolor='none',
                    cmap=cmap)
        cbar = plt.colorbar(orientation='horizontal', extend='both')
        cbar.set_label(texvar+' ('+unit+')')
        plt.clim(*clim)

        plt.contour(*Float.get_griddata_grid(hpids, 'dist_ctd', 'z', var),
                    N=6, colors='k')

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

# %% A single profile showing Ws, Wf and Ww. Uses W hack.

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

for Float in E76, E77]:

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
bathy_file = '../../data/sandwell_bathymetry/topo_17.1.img'
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

# %%

E76_hpids = np.arange(27, 34)
E77_hpids = np.arange(23, 30)

vars = ['rWw', 'U_abs', 'V_abs']
texvars = ['$W_w$', '$U$', '$V$']

for Float, hpids in zip([E76, E77], [E76_hpids, E77_hpids]):
    for var, texvar in zip(vars, texvars):
        pf.depth_profile(Float, hpids, var, hold='on', dlim=[-1500., -100.])
        title_str = "Float {}".format(Float.floatID)
        plt.title(title_str)
        plt.xlabel(texvar+' (m s$^{-1}$)')
        plt.ylabel('Depth (m)')
        plt.legend(hpids, loc=0)

import scipy.optimize as op
import pylab as pyl

def plane_wave(x, A, k, phase, C):
    return A*np.cos(2*np.pi*(k*x + phase)) + C

dz = 5.
dk = 1./dz

E76_hpids = np.arange(27, 34)
E77_hpids = np.arange(23, 30)

var_names = ['rWw']  #['rU_abs', 'rV_abs', 'rWw']

for Float, hpids in zip([E76, E77], [E76_hpids, E77_hpids]):

    for pfl in Float.get_profiles(hpids):

        z = getattr(pfl, 'rz')

        for name in var_names:

            var = getattr(pfl, name)
            N = var.size

            # Try fitting plane wave.
            popt, __ = op.curve_fit(plane_wave, z, var,
                                    p0=[0.1, 1.6e-3, 0., 0.])
            mfit = popt[1]

            plt.figure()
            plt.subplot(1, 2, 1)

            plt.plot(var, z, plane_wave(z, *popt), z)
            plt.xticks(rotation=45)
            plt.xlabel('$W_w$ (m s$^{-1}$)')
            plt.ylabel('$z$ (m)')
            title = ("{}: Float {}, profile {}\nm_fit {:1.1e} m$^{{-1}}$"
                     "\nlambda_fit {:4.0f} m"
                     "").format(name, Float.floatID, pfl.hpid[0], mfit, 1./mfit)
            plt.title(title)

            plt.subplot(1, 2, 2)
            # Try calculating power spectral density.
            Pzz, ms = plt.psd(var, NFFT=N//2, Fs=dk,
                              noverlap=int(0.2*N//2),
                              detrend=pyl.detrend_linear)
            plt.gca().set_xscale('log')
            mmax = ms[Pzz.argmax()]
            title = ("{}: Float {}, profile {}\nm_max {:1.1e} m$^{{-1}}$"
                     "\nlambda_max {:4.0f} m"
                     "").format(name, Float.floatID, pfl.hpid, mmax, 1./mmax)
            ylim = plt.ylim()
            plt.plot(2*[mmax], ylim, 'r', 2*[mfit], ylim, 'g')
            plt.title(title)
            plt.xlabel('$m$ (m$^{-1}$)')


hpids = np.arange(10, 60)

vars = ['T', 'S', 'rho_1']
units = ['$^\\circ$C', 'PSU', 'kg m$^{-3}$']
texvars = ['$T$', '$S$', r'$\sigma_1$']
Tvals = np.arange(2.5, 5.5, 0.5)
Svals = np.arange(34.1, 34.6, 0.1)
rhovals = np.arange(1031.5, 1032.3, 0.1)
vals = [Tvals, Svals, rhovals]

for Float in [E76, E77]:
    __, idxs = Float.get_profiles(hpids, ret_idxs=True)
    dist = Float.dist[idxs]
    for var, texvar, unit, val in zip(vars, texvars, units, vals):

        __, __, var_1g = Float.get_interp_grid(hpids, val, var, 'z')
        plt.figure()
        plt.plot(dist, var_1g.T)
        plt.xlabel('Distance (km)')
        plt.ylabel('Depth (m)')
        plt.legend(str(val).strip('[]').split(), loc=1)
        title_str = ("Float {}, {} ({})").format(Float.floatID, texvar, unit)
        plt.title(title_str)


def plane_wave2(params, x):
    A, k, m, om = params
    return A*np.cos(2*np.pi*(k*x[:,0] + m*x[:,1] + om*x[:,2]))


def cost(params, data, func, y):
    return (func(params, data) - y).flatten()


res = []

E76_hpids = np.arange(27, 34) # np.arange(31, 33)
E77_hpids = np.arange(23, 30) # np.arange(26, 28)

for Float, hpids in zip([E76, E77], [E76_hpids, E77_hpids]):

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)

    z = Float.rz[:,idxs].flatten(order='F')
    x = Float.rdist_ctd[:, idxs].flatten(order='F')*1000.
    t = Float.rUTC[:, idxs].flatten(order='F')*86400.
    W = Float.rWw[:, idxs].flatten(order='F')

    nans = np.isnan(z) | np.isnan(x) | np.isnan(t) | np.isnan(W) | (z > -200)
    data = np.array([x[~nans], z[~nans], t[~nans]]).T
    W = W[~nans]

    x0 = [0.15, 0.001, 0.003, 0.]
    fit = op.leastsq(cost, x0=x0, args=(data, plane_wave2, W))[0]
    print(fit)
    res.append(fit)

    Wm = plane_wave2(fit, data)
    Wm0 = plane_wave2(x0, data)

    plt.figure()
    plt.subplot(3, 1, 1)
    plt.plot(data[:,0], Wm, data[:,0], W)#, data[:,0], Wm0)
    plt.subplot(3, 1, 2)
    plt.plot(Wm, data[:,1], W, data[:,1])#, Wm0, data[:,1])
    plt.subplot(3, 1, 3)
    plt.plot(data[:,2], Wm, data[:,2], W)#, data[:,2], Wm0)
