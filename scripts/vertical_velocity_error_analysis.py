# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 16:27:36 2014

@author: jc3e13
"""

import numpy as np
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
import os
import sys

lib_path = os.path.abspath('../modules')
if lib_path not in sys.path:
    sys.path.append(lib_path)

import emapex
import plotting_functions as pf

try:
    print("Floats {} and {}.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.load(4976)
    E76.generate_regular_grids()
    E77 = emapex.load(4977)
    E77.generate_regular_grids()

# Figure save path.
sdir = os.path.join('..', 'figures', 'vertical_velocity_analysis')
if not os.path.exists(sdir):
    os.makedirs(sdir)
# Universal figure font size.
matplotlib.rc('font', **{'size': 9})


# %% ##########################################################################

pfl = E77.get_profiles(82)
fig = plt.figure(figsize=(3.2, 5))
plt.plot(100.*pfl.Ww, pfl.z, label='$W_{water}$', color='black')
plt.plot(100.*pfl.Wz, pfl.z, label='$W_{float}$', color='grey')
plt.plot(100.*pfl.Ws, pfl.z, label='$W_{steady}$', color='red')
plt.legend(loc=0)
plt.xlabel('$W$ (cm s$^{-1}$)')
plt.ylabel('$z$ (m)')
plt.grid()
plt.title('Float {:}, Profile {:}'.format(int(pfl.floatID), int(pfl.hpid)))
pf.my_savefig(fig, 'example', 'w_model', sdir, ftype='png', fsize='single_col')

# %% ##########################################################################
scaling='density'

# Generate float vertical velocity with random depth error.
# Fake profile
dze = 5.
dt = dze/0.12  # 0.12 is a tpyical descent/ascent speed.
zvals = np.arange(-1400., -100, dze)
N = 50
Pe = np.empty((zvals.size/2, N))

# Use an arbitrary profile for the sizes.
dz = 5.
pfl = E76.get_profiles(155)
surface = pfl.r_z > -100.
z = pfl.r_z[~surface]
Pw76 = np.empty((z.size/2+1, N))
Pw77 = np.empty((z.size/2+1, N))

for i in xrange(N):
    ze = zvals + 0.16*np.random.rand(zvals.size)
    We = np.diff(ze)/dt
    me, Pe[:, i] = sp.signal.periodogram(We, 1./dze, scaling=scaling)

    # Pick a random profile.
    pfl76 = E76.get_profiles(10+i)
    pfl77 = E77.get_profiles(10+i)

    # Chop out the top 100 m.
    Ww76 = pfl76.r_Ww[~surface]
    Ww77 = pfl77.r_Ww[~surface]

    m, Pw76[:, i] = sp.signal.periodogram(Ww76, 1./dz, scaling=scaling)
    __, Pw77[:, i] = sp.signal.periodogram(Ww77, 1./dz, scaling=scaling)

mean_Pe = np.mean(Pe, axis=1)
mean_Pw76 = np.mean(Pw76, axis=1)
mean_Pw77 = np.mean(Pw77, axis=1)

# Remove the 0 wavenumber (infinite wavelength)
m, mean_Pw76, mean_Pw77, me, mean_Pe = m[1:], mean_Pw76[1:], mean_Pw77[1:], \
    me[1:], mean_Pe[1:]

plt.figure()
plt.loglog(m, mean_Pw76)
plt.loglog(m, mean_Pw77)
plt.loglog(me, mean_Pe)
plt.ylabel('Mean power spectral density (m$^3$ s$^{-2}$)')
plt.xlabel('Wavenumber (m$^{-1}$)')
plt.legend(['4976', '4977', 'noise'], loc=0)
plt.grid()

mmin = np.mean([m[mean_Pw76.argmin()], m[mean_Pw77.argmin()]])
pmin = np.mean([mean_Pw76.min(), mean_Pw77.min()])
print('Wavenumber at which noise begins to dominate: {} m-1'.format(mmin))
print('Wavelength: {} m'.format(1./mmin))
print('Power at minima: {:1.1e} m3 s-2'.format(pmin))
print('Error estimate: {} m s-1'.format(np.sqrt(pmin*m.max())))

# Try fitting power law to wavenumbers < 5e-2 m-1
#idxs = m < 5e-2
#pfit = np.polyfit(np.log(m[idxs]), np.log(Pw[idxs]), 1)
#plt.loglog(m, m**pfit[0]/np.mean(m[idxs]**pfit[0]/Pw[idxs]))

# %% What happens to the power spectra over time.

m, PE76 = sp.signal.periodogram(E76.r_Ww.T, 1./dz, scaling=scaling)
__, PE77 = sp.signal.periodogram(E77.r_Ww.T, 1./dz, scaling=scaling)