# -*- coding: utf-8 -*-
"""
Created on Tue Feb  3 10:40:31 2015

@author: jc3e13
"""

import numpy as np
import scipy as sp
import sys
import os
import glob
import matplotlib
import matplotlib.pyplot as plt
import pickle

lib_path = os.path.abspath('../modules')
if lib_path not in sys.path:
    sys.path.append(lib_path)

lib_path = os.path.abspath('../../ocean-tools')
if lib_path not in sys.path:
    sys.path.append(lib_path)

import sandwell
import emapex
import utils
import plotting_functions as pf
import coloured_noise as cn
import GM


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
bp = os.path.join(os.path.expanduser('~'), 'storage', 'smith_sandwell',
                  'topo_*.img')
bf = os.path.abspath(glob.glob(bp)[0])
# Figure save path.
sdir = '../figures/flux analysis'
if not os.path.exists(sdir):
    os.makedirs(sdir)
# Universal figure font size.
matplotlib.rc('font', **{'size': 8})


# %% Horizontal kinetic energy spectrum
hpids = np.arange(50, 150)
zres = 15.

dz = 1.

zmin = -1450.
zmax = -100.

window = 'hanning'
deg = 0

fig, axs = plt.subplots(1, 2, sharex=True, sharey='row', figsize=(6.5, 3))

for Float in [E76, E77]:

    z = np.arange(zmin, zmax, dz)
    __, z_, U = Float.get_interp_grid(hpids, z, 'z', 'U_abs')
    __, __, V = Float.get_interp_grid(hpids, z, 'z', 'V_abs')
    __, __, W = Float.get_interp_grid(hpids, z, 'z', 'Ww')
    U = utils.nan_detrend(z_, U, deg)
    V = utils.nan_detrend(z_, V, deg)

    m, PUs = sp.signal.periodogram(U, fs=1./dz, window=window, axis=0)
    __, PVs = sp.signal.periodogram(V, fs=1./dz, window=window, axis=0)
    __, VKE = sp.signal.periodogram(W, fs=1./dz, window=window, axis=0)

    # Chop interpolation noise.
    muse = m < 1./zres
    m, PUs, PVs, VKE = m[muse], PUs[muse, :], PVs[muse, :], VKE[muse, :]

    HKE = (PUs + PVs)/2.

    axs[0].loglog(1./m, np.median(HKE, axis=-1), linewidth=2., alpha=0.7,
              label=Float.floatID)
    axs[1].loglog(1./m, np.median(VKE, axis=-1), linewidth=2., alpha=0.7,
              label=Float.floatID)

# Set up GM spectrum
IWF = GM.GM(2.2e-3, -1.23e-4)
GMHKE = IWF.Sm(m, 'horiz_vel')
GMVKE = IWF.Sm(m, 'vert_vel')
axs[0].loglog(1./m, GMHKE, linewidth=2., label='GM')
axs[1].loglog(1./m, GMVKE, linewidth=2., label='GM')

# Set up noise estimate
Unoise = cn.noise(len(z), dz, -2)
Unoise *= 0.03/np.std(Unoise)
__, PUnoise = sp.signal.periodogram(Unoise, fs=1./dz, window=window)
PUnoise = PUnoise[muse]
axs[0].loglog(1./m, PUnoise, linewidth=2., label='noise')

Wnoise = cn.noise(len(z), dz, -2)
Wnoise *= 0.005/np.std(Wnoise)
__, PWnoise = sp.signal.periodogram(Wnoise, fs=1./dz, window=window)
PWnoise = PWnoise[muse]
axs[1].loglog(1./m, PWnoise, linewidth=2., label='noise')


axs[0].legend(loc=0)

for ax in axs:
    ax.set_xlim(10, 2000)
    ax.set_xlabel('Vertical wavelength (m)')

axs[0].set_ylabel('Horizontal kinetic energy density (m$^{3}$ s$^{-2}$)')
axs[1].set_ylabel('Vertical kinetic energy density (m$^{3}$ s$^{-2}$)')

pf.my_savefig(fig, 'both', 'UV_spec_GM', sdir, ftype=('png', 'pdf'),
              fsize='double_col')

# %%

np.random.seed(1029301)
# The portions of the profiles that contain the wave. Found by eye.
zlims = {4976: {25: (-1600, -100),
                26: (-1600, -100),
                27: (-1600, -100),
                28: (-1600, -100),
                29: (-1600, -200),  # 1
                30: (-1000, -200),  # 2
                31: (-1600, -400),  # 3
                32: (-1600, -400),  # 4
                33: (-1600, -100),
                34: (-1600, -100),
                35: (-1600, -100),
                36: (-1600, -100)},
         4977: {20: (-1600, -100),
                21: (-1600, -100),
                22: (-1600, -100),
                23: (-1600, -100),
                24: (-1600, -200),  # 1
                25: (-1400, -600),  # 2
                26: (-1600, -600),  # 3
                27: (-1200, -200),  # 4
                28: (-1600, -100),
                29: (-1600, -100),
                30: (-1600, -100),
                31: (-1600, -100)}}

#hpids_76 = np.array([28, 29, 30, 31, 32, 33, 34, 35])
#hpids_77 = np.array([23, 24, 25, 26, 27, 28, 29, 30])
hpids_76 = np.array([31, 32])
hpids_77 = np.array([26, 27])

N76 = hpids_76.size
N77 = hpids_77.size
N = np.sum((hpids_76.size, hpids_77.size))

rho0 = 1025.
# Detrend degree
deg = 1

#fig, ax1 = plt.subplots(1, 1)
#ax1.set_xlabel('Cov$(w, u)$')
#ax1.set_ylabel('$z$ (m)')

#fig, ax2 = plt.subplots(1, 1)
#ax2.set_xlabel('$u$ (m s$^{-1}$)')
#ax2.set_ylabel('$v$ (m s$^{-1}$)')
#ax2.set_xlim(-0.3, 0.3)
#ax2.set_ylim(-0.3, 0.3)


uwbar = np.zeros((2, N/2))
vwbar = np.zeros((2, N/2))
pwbar = np.zeros((2, N/2))
Uuwbar = np.zeros((2, N/2))
Vvwbar = np.zeros((2, N/2))
tau = np.zeros((2, N/2))
E = np.zeros((2, N/2))

Nreps = 200

uwbar_errs = np.zeros((2, N/2, Nreps))
vwbar_errs = np.zeros((2, N/2, Nreps))
pwbar_errs = np.zeros((2, N/2, Nreps))
tau_errs = np.zeros((2, N/2, Nreps))
E_errs = np.zeros((2, N/2, Nreps))


i = 0
for Float, hpids in zip([E76, E77], [hpids_76, hpids_77]):
    j = 0
    for hpid in hpids:

        zmin, zmax = zlims[Float.floatID][hpid]
        pfl = Float.get_profiles(hpid)

        use = (pfl.z > zmin) & (pfl.z < zmax)
        useef = (pfl.zef > zmin) & (pfl.zef < zmax)

        t = pfl.UTC[use]
        tef = pfl.UTCef[useef]
        z = pfl.zef[useef]
        u = pfl.U[useef]
        v = pfl.V[useef]
        w = pfl.Ww[use]
        b = pfl.b[use]
        pp = pfl.Pprime[use]
        N2 = pfl.N2_ref[use]

        w = np.interp(tef, t, w)
        b = np.interp(tef, t, b)
        pp = np.interp(tef, t, pp)
        N2 = np.interp(tef, t, N2)
        N2mean = np.mean(N2)

        U = np.mean(u)
        V = np.mean(v)

        u = utils.nan_detrend(z, u, deg)
        v = utils.nan_detrend(z, v, deg)
        pp = utils.nan_detrend(z, pp, 0)

        Ndata = tef.size

        # Correct for sign error in pp for downward profiles with large aspect
        # ratio. w is proxy for aspect.
        if (np.max(np.abs(w)) > 0.05) and (hpid % 2 != 0.):
            pp *= -1.

        if False:

            fig, axs = plt.subplots(1, 5, sharey='row')
            axs[0].plot(u, z, label='u')
            axs[0].set_xlabel('u')
            axs[1].plot(v, z, label='v')
            axs[1].set_xlabel('v')
            axs[2].plot(w, z, label='w')
            axs[2].set_xlabel('w')
            axs[3].plot(b, z, label='b')
            axs[3].set_xlabel('b')
            axs[4].plot(pp, z, label='pp')
            axs[4].set_xlabel('pp')
            axs[0].set_ylim(-1500., 0.)
            plt.title("Float {}. hpid {}.".format(pfl.floatID, pfl.hpid[0]))

        DT = np.max(tef) - np.min(tef)

        print("\nFloat {}. hpid {}.".format(pfl.floatID, pfl.hpid[0]))

        uwbar[i, j] = rho0*sp.integrate.trapz(w*u, tef)/DT
        vwbar[i, j] = rho0*sp.integrate.trapz(w*v, tef)/DT
        pwbar[i, j] = rho0*sp.integrate.trapz(w*pp, tef)/DT

        Uuwbar[i, j] = -U*uwbar[i, j]
        Vvwbar[i, j] = -V*vwbar[i, j]

        Efluxz = Uuwbar[i, j] + Vvwbar[i, j]

        tau[i, j] = np.sqrt(uwbar[i, j]**2 + vwbar[i, j]**2)

        kinetic = 0.5*rho0*sp.integrate.trapz(u**2 + v**2 + w**2, tef)/DT
        potential = 0.5*rho0*sp.integrate.trapz(b**2/N2mean, tef)/DT
        E[i, j] = kinetic + potential

        uwbar_e = np.zeros(Nreps)
        vwbar_e = np.zeros(Nreps)
        pwbar_e = np.zeros(Nreps)
        tau_e = np.zeros(Nreps)
        E_e = np.zeros(Nreps)
        # ERROR ESTIMATES
        for ii in xrange(Nreps):

            noise1 = cn.noise(Ndata, DT*60.*60.*24./Ndata, -2)
            noise1 *= 0.03/np.std(noise1)
            noise2 = cn.noise(Ndata, DT*60.*60.*24./Ndata, -2)
            noise2 *= 0.03/np.std(noise2)
            noise3 = cn.noise(Ndata, DT*60.*60.*24./Ndata, -2)
            noise3 *= 0.005/np.std(noise3)
            u_e = u + noise1
            v_e = v + noise2
            w_e = w + noise3

            uwbar_e[ii] = rho0*sp.integrate.trapz(w_e*u_e, tef)/DT
            vwbar_e[ii] = rho0*sp.integrate.trapz(w_e*v_e, tef)/DT
            pwbar_e[ii] = rho0*sp.integrate.trapz(w_e*pp, tef)/DT

            tau_e[ii] = np.sqrt(uwbar_e[ii]**2 + vwbar_e[ii]**2)

            kinetic_e = 0.5*rho0*sp.integrate.trapz(u_e**2 + v_e**2 + w_e**2, tef)/DT
            E_e[ii] = kinetic_e + potential

        uwbar_errs[i, j, :] = uwbar_e
        vwbar_errs[i, j, :] = vwbar_e
        pwbar_errs[i, j, :] = pwbar_e
        tau_errs[i, j, :] = tau_e
        E_errs[i, j, :] = E_e

        uwbar_mean = np.mean(uwbar_errs, axis=-1)
        vwbar_mean = np.mean(vwbar_errs, axis=-1)
        pwbar_mean = np.mean(pwbar_errs, axis=-1)
        tau_mean = np.mean(tau_errs, axis=-1)
        E_mean = np.mean(E_errs, axis=-1)

        uwbar_std = np.std(uwbar_errs, axis=-1)
        vwbar_std = np.std(vwbar_errs, axis=-1)
        pwbar_std = np.std(pwbar_errs, axis=-1)
        tau_std = np.std(tau_errs, axis=-1)
        E_std = 2*np.std(E_errs, axis=-1) #  2x because error in b, and equipartition of energy

        print("Reynolds stress: {:1.2f} +/- {:1.2f} N m-2".format(tau[i, j], tau_std[i, j]))
        print("Components: ({:1.2f}, {:1.2f}) +/- ({:1.2f}, {:1.2f}) N m-2".format(uwbar[i, j], vwbar[i, j], uwbar_std[i, j], vwbar_std[i, j]))
        print("Energy density: {:1.2f} +/- {:1.2f} J m-3".format(E[i, j], E_std[i, j]))
        print("Standard deviation of pressure perturbation: "
              "{:1.2f} m2 s-2".format(np.std(pp)))
        print("Vertical energy flux: {:1.2f} +/- {:1.2f} W m-2".format(pwbar[i, j], pwbar_std[i, j]))
        print("Vertical energy from from (Uuw, Vvw): ({:1.2f}, {:1.2f})"
              ", total {:1.2f} W m-2 ".format(Uuwbar[i, j], Vvwbar[i, j],
                                              Efluxz))

#        print(1025.*sp.integrate.trapz(v*u, z)/(z[0]-z[-1]))

#        ax1.plot(u*w, z, v*w, z)

#        ax2.quiver(u[:-1], v[:-1], u[1:]-u[:-1], v[1:]-v[:-1],
#                   scale_units='xy', angles='xy', scale=1,
#                   color=mpl.cm.Set1(i/N))

        j += 1
    i += 1

with open('/noc/users/jc3e13/storage/processed/Edens.p', 'wb') as f:
    pickle.dump(E_errs, f)
with open('/noc/users/jc3e13/storage/processed/Eflux.p', 'wb') as f:
    pickle.dump(pwbar_errs, f)
with open('/noc/users/jc3e13/storage/processed/Mflux.p', 'wb') as f:
    pickle.dump(tau_errs, f)

# %%

fig, axs = plt.subplots(1, 3, figsize=(3.125, 3))
fig.tight_layout()
#for ax in axs[1:]:
#    ax.yaxis.tick_right()
#    ax.yaxis.set_ticks_position('both')
#    ax.yaxis.set_label_position('right')

colors = ['blue', 'green', 'red', 'purple']

Ns = [N76, N77]
k = 0
for i in xrange(2):
    for j in xrange(Ns[i]):

        colprops={'color': colors[k]}

        axs[0].boxplot(E_errs[i, j, :], boxprops=colprops,
                       whiskerprops=colprops, capprops=colprops,
                       medianprops=colprops, showfliers=False,
                       labels=[''])
        axs[1].boxplot(pwbar_errs[i, j, :], boxprops=colprops,
                       whiskerprops=colprops, capprops=colprops,
                       medianprops=colprops, showfliers=False,
                       labels=[''])
        axs[2].boxplot(tau_errs[i, j, :], boxprops=colprops,
                       whiskerprops=colprops, capprops=colprops,
                       medianprops=colprops, showfliers=False,
                       labels=[''])

        k += 1

# The legend
labels = ['4976 P 31', '4976 P 32', '4977 P 26', '4977 P 27']
yloc = [38., 36., 34., 32.]
for i in xrange(4):
    fontprops={'color': colors[i], 'fontsize': 6}
    axs[0].text(0.55, yloc[i], labels[i], fontdict=fontprops)

labelpad = -1.5
axs[0].set_ylabel('Energy density $E$ (J m$^{-3}$)', labelpad=labelpad)
axs[0].set_ylim(0., 40.)
axs[0].grid(axis='y')
axs[1].set_ylabel(r"Vertical energy flux $\overline{p'w'}$ (W m$^{-2}$)", labelpad=labelpad)
axs[1].grid(axis='y')
axs[1].set_ylim(0., 2.)
axs[2].set_ylabel('Vertical momentum flux $F_M^{(z)}$ (N m$^{-2}$)', labelpad=labelpad)
axs[2].grid(axis='y')
axs[2].set_ylim(0., 12.)

pf.my_savefig(fig, 'both', 'observed_flux_boxplots', sdir, ftype=('png', 'pdf'),
              fsize='single_col')

# %% More figures

fig, axs = plt.subplots(4, 1, sharex=True, figsize=(3.125, 5))

ds = []

for Float, hpids in zip([E76, E77], [hpids_76, hpids_77]):

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)

    lon = getattr(Float, 'lon_start')[idxs]
    lat = getattr(Float, 'lat_start')[idxs]
    d = getattr(Float, 'dist')[idxs]

    bathy = sandwell.interp_track(lon, lat, bf)

    d -= d[bathy.argmax()]
    ds.append(d.copy())


axs[0].plot(ds[0], E[0, :], label='4976')
colprops = {'color': 'b'}
axs[0].boxplot(E_errs[0, :, :].T, positions=ds[0], boxprops=colprops,
               whiskerprops=colprops, capprops=colprops,
               medianprops=colprops, showfliers=False)
axs[0].plot(ds[1], E[1, :], label='4977')
colprops = {'color': 'g'}
axs[0].boxplot(E_errs[1, :, :].T, positions=ds[1], boxprops=colprops,
               whiskerprops=colprops, capprops=colprops,
                   medianprops=colprops, showfliers=False)
axs[0].set_ylabel("$E$ (J m$^{-3}$)")
axs[0].legend(loc=0)

axs[1].plot(ds[0], pwbar[0, :])
colprops = {'color': 'b'}
axs[1].boxplot(pwbar_errs[0, :, :].T, positions=ds[0], boxprops=colprops,
               whiskerprops=colprops, capprops=colprops,
               medianprops=colprops, showfliers=False)
axs[1].plot(ds[1], pwbar[1, :])
colprops = {'color': 'g'}
axs[1].boxplot(pwbar_errs[1, :, :].T, positions=ds[1], boxprops=colprops,
               whiskerprops=colprops, capprops=colprops,
               medianprops=colprops, showfliers=False)
axs[1].set_ylabel(r"$\overline{p'w'}$ (W m$^{-2}$)")
axs[1].hlines(0., *axs[1].get_xlim(), color='k')

#axs[2].plot(ds[0], uwbar[0, :], 'b:')
#axs[2].plot(ds[0], vwbar[0, :], 'b--')
axs[2].plot(ds[0], tau[0, :], 'b-')
colprops = {'color': 'b'}
axs[2].boxplot(tau_errs[0, :, :].T, positions=ds[0], boxprops=colprops,
               whiskerprops=colprops, capprops=colprops,
               medianprops=colprops, showfliers=False)
#axs[2].plot(ds[1], uwbar[1, :], 'g:')
#axs[2].plot(ds[1], vwbar[1, :], 'g--')
axs[2].plot(ds[1], tau[1, :], 'g-')
colprops = {'color': 'g'}
axs[2].boxplot(tau_errs[1, :, :].T, positions=ds[1], boxprops=colprops,
               whiskerprops=colprops, capprops=colprops,
               medianprops=colprops, showfliers=False)

axs[2].set_ylabel(r"$\tau$ (N m$^{-2}$)")

axs[3].set_xlabel('Distance from ridge top (km)', labelpad=0.06)
axs[3].set_ylabel('$z$ (m)')

axs[3].fill_between(d, bathy, np.nanmin(bathy), color='black',
                    linewidth=2)
axs[3].set_ylim(np.nanmin(bathy), np.nanmax(bathy))
axs[3].set_xlim(np.nanmin(d), np.nanmax(d))
axs[3].set_xticks([-10., 0., 10., 20.])
axs[3].set_xticklabels([-10., 0., 10., 20.])

pf.my_savefig(fig, 'both', 'fluxes', sdir, ftype=('png', 'pdf'), fsize='single_col')

## Just the horizontal component of the momentum flux.
#
#fig, axs = plt.subplots(2, 1, sharex=True, figsize=(3.125, 4))
#
#axs[0].plot(ds[0], uwbar[0, :], 'b', label='4976')
#axs[0].plot(ds[1], uwbar[1, :], 'g', label='4977')
#axs[0].set_ylabel("$<u'w'>$ (N m$^{-2}$)")
#axs[0].legend(loc=0)
#
#axs[1].set_xlabel('Distance from ridge top (km)', labelpad=0.06)
#axs[1].set_ylabel('$z$ (m)')
#axs[1].fill_between(d, bathy, np.nanmin(bathy), color='black',
#                    linewidth=2)
#axs[1].set_ylim(np.nanmin(bathy), np.nanmax(bathy))
#axs[1].set_xlim(np.nanmin(d), np.nanmax(d))
#
#pf.my_savefig(fig, 'both', 'M_flux_only', sdir, ftype='png',
#              fsize='single_col')