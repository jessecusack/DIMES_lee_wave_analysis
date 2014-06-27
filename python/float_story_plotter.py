"""
Plot a LOT of things.
"""

import numpy as np
import matplotlib.pyplot as plt
import plotting_functions as pf
import emapex
import scipy.optimize as op
import pylab as pyl
import os

reload(emapex)

E76 = emapex.EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4976)
E76.apply_w_model('../../data/EM-APEX/4976_fix_p0k0M_fit_info.p')
E76.apply_strain('../../data/EM-APEX/4976_N2_ref_300dbar.p')

E77 = emapex.EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4977)
E77.apply_w_model('../../data/EM-APEX/4977_fix_p0k0M_fit_info.p')
E77.apply_strain('../../data/EM-APEX/4977_N2_ref_300dbar.p')

plot_dir = '../figures/large_wave_analysis'
ftype = 'png'


def my_savefig(fid, fname):
    fname = str(fid) + '_' + fname
    fname = os.path.join(plot_dir, fname) + '.' + ftype
    plt.savefig(fname, bbox_inches='tight')

# %%

hpids = np.arange(1, 600)
z_vals = np.arange(-1500., 0., 10.)

vars = ['T', 'S', 'rho_1']
texvars = ['$T$', '$S$', r'$\sigma_1$']
units = [r'$^\circ$C', 'PSU', 'kg m$^{-3}$']

for Float in [E76, E77]:
    for var, texvar, unit in zip(vars, texvars, units):
        __, z, C = Float.get_interp_grid(hpids, z_vals, 'z', var)
        z = z.flatten(order='F')
        C = C.flatten(order='F')
        __, __, d = Float.get_interp_grid(hpids, z_vals, 'z',
                                          'dist_ctd_data')
        d = d.flatten(order='F')
        plt.figure(figsize=(6, 4))
        plt.scatter(d, z, c=C, edgecolor='none')
        plt.ylim(np.min(z), np.max(z))
        plt.xlim(np.min(d), np.max(d))
        plt.xlabel('Distance (km)')
        plt.ylabel('Depth (m)')
        title_str = ("Float {}").format(Float.floatID)
        plt.title(title_str)
        cbar = plt.colorbar(orientation='horizontal')
        cbar.set_label(texvar+' ('+unit+')')
        my_savefig(Float.floatID, var)


import matplotlib.colors as mcolors
rwb = mcolors.LinearSegmentedColormap.from_list(name='red_white_blue',
                                                colors=[(0, 0, 1),
                                                        (1, 1, 1),
                                                        (1, 0, 0)],
                                                N=40,
                                                )

hpids = np.arange(1, 600)
z_vals = np.arange(-1500., 0., 10.)

for Float in [E76, E77]:
    __, z, Ww = Float.get_interp_grid(hpids, z_vals, 'z', 'rWw')
    z = z.flatten(order='F')
    Ww = Ww.flatten(order='F')
    __, __, d = Float.get_interp_grid(hpids, z_vals, 'z',
                                      'dist_ctd_data')
    d = d.flatten(order='F')
    plt.figure()
    plt.scatter(d, z, c=Ww, edgecolor='none', cmap=rwb)
    plt.ylim(np.min(z), np.max(z))
    plt.xlim(np.min(d), np.max(d))
    plt.xlabel('Distance (km)')
    plt.ylabel('Depth (m)')
    plt.xlim(np.min(Float.dist), np.max(Float.dist))
    title_str = ("Float {}").format(Float.floatID)
    plt.title(title_str)
    cbar = plt.colorbar(orientation='horizontal', extend='both')
    cbar.set_label('$W_w$ (m s$^{-1}$)')
    plt.clim(-0.15, 0.15)

hpids = np.arange(10, 60)
z_vals = np.arange(-1500., 0., 10.)

for Float in [E76, E77]:
    __, z, Ww = Float.get_interp_grid(hpids, z_vals, 'z', 'rWw')
    z = z.flatten(order='F')
    Ww = Ww.flatten(order='F')
    __, __, d = Float.get_interp_grid(hpids, z_vals, 'z',
                                      'dist_ctd_data')
    d = d.flatten(order='F')
    plt.figure()
    plt.scatter(d, z, c=Ww, edgecolor='none', cmap=rwb)
    plt.ylim(np.min(z), np.max(z))
    plt.xlim(np.min(d), np.max(d))
    plt.xlabel('Distance (km)')
    plt.ylabel('Depth (m)')
    title_str = ("Float {}").format(Float.floatID)
    plt.title(title_str)
    cbar = plt.colorbar(orientation='horizontal', extend='both')
    cbar.set_label('$W_w$ (m s$^{-1}$)')
    plt.clim(-0.15, 0.15)

    pf.bathy_along_track(Float, hpids)
    plt.xlim(np.min(d), np.max(d))

    pf.track_on_bathy(Float, hpids)


for Float in [E76, E77]:
    for comp in ['U_abs', 'V_abs']:
        __, z, V = Float.get_interp_grid(hpids, z_vals, 'z', comp)
        z = z.flatten(order='F')
        V = V.flatten(order='F')
        __, __, d = Float.get_interp_grid(hpids, z_vals, 'z',
                                          'dist_ctd_data')
        d = d.flatten(order='F')
        plt.figure()
        plt.scatter(d, z, c=V, edgecolor='none', cmap=rwb)
        plt.ylim(np.min(z), np.max(z))
        plt.xlim(np.min(d), np.max(d))
        plt.xlabel('Distance (km)')
        plt.ylabel('Depth (m)')
        title_str = ("Float {}").format(Float.floatID)
        plt.title(title_str)
        cbar = plt.colorbar(orientation='horizontal', extend='both')
        cbar.set_label('$'+comp[0]+'$ (m s$^{-1}$)')
        plt.clim(-1.5, 1.5)

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
    x = Float.rdist_ctd_data[:, idxs].flatten(order='F')*1000.
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
