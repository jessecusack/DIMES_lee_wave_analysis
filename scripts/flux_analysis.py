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
import matplotlib as mpl
import matplotlib.pyplot as plt

lib_path = os.path.abspath('../modules')
if lib_path not in sys.path:
    sys.path.append(lib_path)

import sandwell
import emapex
import utils
import plotting_functions as pf


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
sdir = '../figures/flux analysis'
if not os.path.exists(sdir):
    os.makedirs(sdir)
# Universal figure font size.
matplotlib.rc('font', **{'size': 9})

# %%

# The portions of the profiles that contain the wave. Found by eye.
zlims = {4976: {25: (-1600, -100),
                26: (-1600, -100),
                27: (-1600, -100),
                28: (-1600, -100),
                29: (-1600, -200),  # 1
                30: (-1000, -200),  # 2
                31: (-1600, -600),  # 3
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

hpids_76 = np.array([25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36])
hpids_77 = np.array([20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31])

#hpids_76 = np.array([31])
#hpids_77 = np.array([26])

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

i = 0
for Float, hpids in zip([E76, E77], [hpids_76, hpids_77]):
    j = 0
    for hpid in hpids:

        zmin, zmax = zlims[Float.floatID][hpid]
        pfl = Float.get_profiles(hpid)

#        fig, axs = plt.subplots(1, 2)
#        axs[0].plot(pfl.Pprime, pfl.z)
#        axs[1].plot(pfl.Ww, pfl.z)

        use = (pfl.z > zmin) & (pfl.z < zmax)
        useef = (pfl.zef > zmin) & (pfl.zef < zmax)

#        fig = plt.figure()
#        plt.plot(pfl.U_abs[useef], pfl.zef[useef], 'b-', pfl.U_abs[~useef],
#                 pfl.zef[~useef], 'b--')
#        plt.plot(pfl.V_abs[useef], pfl.zef[useef], 'g-', pfl.V_abs[~useef],
#                 pfl.zef[~useef], 'g--')
#        plt.plot(pfl.Ww[use], pfl.z[use], 'r-', pfl.Ww[~use], pfl.z[~use],
#                 'r--')

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
        pp = utils.nan_detrend(z, pp, deg)

#        plt.figure()
#        plt.plot(u, z, v, z, w, z)
#        plt.xlim(-0.4, 0.4)
#        plt.ylim(-1500., 0.)
#        plt.title("Float {}. hpid {}.".format(pfl.floatID, pfl.hpid[0]))

        DT = np.max(tef) - np.min(tef)

        print("\nFloat {}. hpid {}.".format(pfl.floatID, pfl.hpid[0]))

        uwbar[i, j] = rho0*sp.integrate.trapz(w*u, tef)/DT
        vwbar[i, j] = rho0*sp.integrate.trapz(w*v, tef)/DT
        pwbar[i, j] = rho0*sp.integrate.trapz(w*pp, tef)/DT

        Uuwbar[i, j] = -U*uwbar[i, j]
        Vvwbar[i, j] = -V*vwbar[i, j]

        tau[i, j] = np.sqrt(uwbar[i, j]**2 + vwbar[i, j]**2)

        kinetic = 0.5*rho0*sp.integrate.trapz(u**2 + v**2 + w**2, tef)/DT
        potential = 0.5*rho0*sp.integrate.trapz(b**2/N2mean, tef)/DT
        E[i, j] = kinetic + potential

        print("Reynolds stress: {:1.2f} N m-2".format(tau[i, j]))
        print("Components: ({:1.2f}, {:1.2f}) N m-2".format(uwbar[i, j], vwbar[i, j]))
        print("Energy density: {:1.2f} J m-3".format(E[i, j]))
        print("Standard deviation of pressure perturbation: "
              "{:1.2f} m2 s-2".format(np.std(pp)))
        print("Vertical energy flux: {:1.2f} W m-2".format(pwbar[i, j]))
        print("Vertical energy from from (Uuw, Vvw): ({:1.2f}, {:1.2f})"
              "W m-2 ".format(Uuwbar[i, j], Vvwbar[i, j]))

#        print(1025.*sp.integrate.trapz(v*u, z)/(z[0]-z[-1]))

#        ax1.plot(u*w, z, v*w, z)

#        ax2.quiver(u[:-1], v[:-1], u[1:]-u[:-1], v[1:]-v[:-1],
#                   scale_units='xy', angles='xy', scale=1,
#                   color=mpl.cm.Set1(i/N))

        j += 1
    i += 1

fig, axs = plt.subplots(3, 1, sharex=True, figsize=(3.125, 5))

ds = []

for Float, hpids in zip([E76, E77], [hpids_76, hpids_77]):

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)

    lon = getattr(Float, 'lon_start')[idxs]
    lat = getattr(Float, 'lat_start')[idxs]
    d = getattr(Float, 'dist')[idxs]

    bathy = sandwell.interp_track(lon, lat, bf)

    d -= d[bathy.argmax()]
    ds.append(d.copy())

axs[2].set_xlabel('Distance from ridge top (km)', labelpad=0.06)
axs[2].set_ylabel('$z$ (m)')

axs[2].fill_between(d, bathy, np.nanmin(bathy), color='black',
                    linewidth=2)
axs[2].set_ylim(np.nanmin(bathy), np.nanmax(bathy))
axs[2].set_xlim(np.nanmin(d), np.nanmax(d))

axs[1].plot(ds[0], uwbar[0, :], 'b:')
axs[1].plot(ds[0], vwbar[0, :], 'b--')
axs[1].plot(ds[0], tau[0, :], 'b-')
axs[1].plot(ds[1], uwbar[1, :], 'g:')
axs[1].plot(ds[1], vwbar[1, :], 'g--')
axs[1].plot(ds[1], tau[1, :], 'g-')

axs[1].set_ylabel("$<u'w'>$ (N m$^{-2}$)")

axs[0].plot(ds[0], Uuwbar[0, :] + Vvwbar[0, :], label='4976')
axs[0].plot(ds[1], Uuwbar[1, :] + Vvwbar[1, :], label='4977')
axs[0].set_ylabel("$U<u'w'>$ (W m$^{-2}$)")
axs[0].legend(loc=0)

pf.my_savefig(fig, 'both', 'fluxes', sdir, ftype='png', fsize='single_col')
