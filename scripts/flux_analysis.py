# -*- coding: utf-8 -*-
"""
Created on Tue Feb  3 10:40:31 2015

@author: jc3e13
"""

import numpy as np
import scipy as sp
import sys
import os
import matplotlib as mpl
import matplotlib.pyplot as plt

lib_path = os.path.abspath('../modules')
if lib_path not in sys.path:
    sys.path.append(lib_path)

import emapex
import utils


try:
    print("Floats {} and {} exist!.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.load(4976)
    E77 = emapex.load(4977)


# %%

# The portions of the profiles that contain the wave. Found by eye.
zlims = {4976: {29: (-1600, -200),
                30: (-1000, -200),
                31: (-1600, -600),
                32: (-1600, -400)},
         4977: {24: (-1600, -200),
                25: (-1400, -600),
                26: (-1600, -600),
                27: (-1200, -200)}}

#hpids_76 = np.array([29, 30, 31, 32])
#hpids_77 = np.array([24, 25, 26, 27])

hpids_76 = np.array([31])
hpids_77 = np.array([26])

N = np.sum((hpids_76.size, hpids_77.size))

rho0 = 1025.
# Detrend degree
deg = 2

fig, ax1 = plt.subplots(1, 1)
ax1.set_xlabel('Cov$(w, u)$')
ax1.set_ylabel('$z$ (m)')

fig, ax2 = plt.subplots(1, 1)
ax2.set_xlabel('$u$ (m s$^{-1}$)')
ax2.set_ylabel('$v$ (m s$^{-1}$)')
ax2.set_xlim(-0.3, 0.3)
ax2.set_ylim(-0.3, 0.3)

i = 0.
for Float, hpids in zip([E76, E77], [hpids_76, hpids_77]):
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

        w = np.interp(tef, t, w)
        b = np.interp(tef, t, b)

        u = utils.nan_detrend(z, u, deg)
        v = utils.nan_detrend(z, v, deg)

        DT = np.max(tef) - np.min(tef)

        print("Float {}. hpid {}.".format(pfl.floatID, pfl.hpid))

        uwbar = rho0*sp.integrate.trapz(w*u, tef)/DT
        vwbar = rho0*sp.integrate.trapz(w*v, tef)/DT

        tau = np.sqrt(uwbar**2 + vwbar**2)

        kinetic = 0.5*rho0*sp.integrate.trapz(u**2 + v**2 + w**2, tef)/DT
        potential = 0.5*rho0*sp.integrate.trapz(b**2, tef)/DT
        E = kinetic + potential

        print("Reynolds stress: {:1.2f} N m-2".format(tau))
        print("Components: ({:1.2f}, {:1.2f}) N m-2".format(uwbar, vwbar))
        print("Energy density: {:1.2f} J m-3".format(E))

#        print(1025.*sp.integrate.trapz(v*u, z)/(z[0]-z[-1]))

        ax1.plot(u*w, z, v*w, z)

        ax2.quiver(u[:-1], v[:-1], u[1:]-u[:-1], v[1:]-v[:-1],
                   scale_units='xy', angles='xy', scale=1,
                   color=mpl.cm.Set1(i/N))

        i += 1

ax2.grid()