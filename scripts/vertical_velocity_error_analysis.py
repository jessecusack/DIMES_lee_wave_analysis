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
    E76.generate_regular_grids(dz=2.)
    E77 = emapex.load(4977)
    E77.generate_regular_grids(dz=2.)

# Figure save path.
sdir = os.path.join('..', 'figures', 'vertical_velocity_analysis')
if not os.path.exists(sdir):
    os.makedirs(sdir)
# Universal figure font size.
matplotlib.rc('font', **{'size': 9})


# %% ##########################################################################

pfl = E77.get_profiles(151)
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
dze = 2.
dt = dze/0.13  # 0.12 is a tpyical descent/ascent speed.
zvals = np.arange(-1400., -100, dze)
N = 50
Pe = np.empty((zvals.size/2, N))

# Use an arbitrary profile for the sizes.
dz = 2.
pfl = E76.get_profiles(155)
surface = pfl.r_z > -400.
z = pfl.r_z[~surface]
Pw76 = np.empty((z.size/2+1, N))
Pw77 = np.empty((z.size/2+1, N))

for i in xrange(N):
    ze = zvals + 0.15*np.random.rand(zvals.size)
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

# %% Investigation using the lomb-scargle method without interpolation and
# then with interpolation.

# Choose float and range of hpids.
hpids = np.arange(10, 70)
Float = E76
pfls = Float.get_profiles(hpids)

# Chop ends where model is invalid.
zmax = -30.
zmin = -1400.

# Choose frequencies to estimate power.
tdmax = 9000.
tdmin = 20.
dzmax = 1200.
dzmin = 4.
N = 100
om = np.logspace(np.log10(np.pi*2/tdmax), np.log10(np.pi*2/tdmin), N)
m = np.logspace(np.log10(np.pi*2/dzmax), np.log10(np.pi*2/dzmin), N)

om_pgrams = np.empty((om.size, pfls.size))
m_pgrams = np.empty((m.size, pfls.size))

for i, pfl in enumerate(pfls):
    nans = np.isnan(pfl.Ww)
    w = pfl.Ww[~nans]
    t = pfl.UTC[~nans]*24*60*60
    z = pfl.z[~nans]

    invalid = (z > zmax) & (z < zmin)
    w = w[~invalid]
    t = t[~invalid]
    z = z[~invalid]

    om_pgrams[:, i] = sp.signal.lombscargle(t, w, om)/om
    m_pgrams[:, i] = sp.signal.lombscargle(t, w, m)/m

f = om/(np.pi*2)
m = m/(np.pi*2)

fig1, ax1 = plt.subplots(1, 1)
ax1.loglog(1./f, om_pgrams, color='blue', alpha=0.1)
ax1.loglog(1./f, np.mean(om_pgrams, axis=-1), color='black',
         linewidth=3.)
ax1.set_xlabel('$\omega$ (rad s$^{-1})$')

fig2, ax2 = plt.subplots(1, 1)
plt.loglog(1./m, m_pgrams, color='blue', alpha=0.1)
plt.loglog(1./m, np.mean(m_pgrams, axis=-1), color='black',
         linewidth=3.)
plt.xlabel('$m$ (rad m$^{-1}$)')


######## Using Interpolation #########

# Choose float and range of hpids.
#hpids = np.arange(50, 150)
#Float = E76
#pfls = Float.get_profiles(hpids)
scaling = 'density'

# Chop ends where model is invalid.
zmax = -30.
zmin = -1400.
tmax = 9000.
tmin = 200.
tdmin = 20.
dzmin = 4.

# Interpolation points.
N = 2**11
ti, dt = np.linspace(tmin, tmax, N, retstep=True)
zi, dz = np.linspace(zmin, zmax, N, retstep=True)

om_pgrams = np.empty((ti.size/2 + 1, pfls.size))
m_pgrams = np.empty((zi.size/2 + 1, pfls.size))

for i, pfl in enumerate(pfls):
    nans = np.isnan(pfl.Ww)
    w = pfl.Ww[~nans]
    t = pfl.UTC[~nans]*24*60*60
    z = pfl.z[~nans]

    invalid = (z > zmax) & (z < zmin)
    w = w[~invalid]
    t = t[~invalid]
    z = z[~invalid]

    wit = np.interp(ti, t - t[0], w)
    wiz = np.interp(zi, z, w)

    f, om_pgrams[:, i] = sp.signal.periodogram(wit, fs=1./dt, scaling=scaling)
    m, m_pgrams[:, i] = sp.signal.periodogram(wiz, fs=1./dz, scaling=scaling)

# Chop interpolation noise.
fuse = f < 1./tdmin
muse = m < 1./dzmin
f, om_pgrams = f[fuse], om_pgrams[fuse, :]
m, m_pgrams = m[muse], m_pgrams[muse, :]
om_pgrams[0, :] = 0.
m_pgrams[0, :] = 0.

# Convert to radian units.
ax1.loglog(1./f, om_pgrams, color='red', alpha=0.1)
ax1.loglog(1./f, np.mean(om_pgrams, axis=-1), color='black',
         linewidth=3.)
ax1.set_xlabel('$T$ (s)')
ax1.grid()

ax2.loglog(1./m, m_pgrams, color='red', alpha=0.1)
ax2.loglog(1./m, np.mean(m_pgrams, axis=-1), color='black',
           linewidth=3.)
ax2.set_xlabel('$\lambda$ (m)')
ax2.grid()


def tick_func(ticks):
    return ["{:.0f}".format(tick) for tick in ticks]

ax1.set_xticklabels(tick_func(ax1.get_xticks()))
ax2.set_xticklabels(tick_func(ax2.get_xticks()))

#ax1t = ax1.twiny()
#ax1t.set_xticks(np.log10(ax1.get_xticks()))
#ax1t.set_xticklabels(tick_func(1./ax1.get_xticks()))
##ax1t.set_xbound(ax1.get_xbound())
#ax1t.set_xlabel("T (s)")
#
#ax2t = ax2.twiny()
#L_ticks = np.logspace(4, 0, 5)
#ax2t.set_xticks(np.log10(1./L_ticks))
#ax2t.set_xticklabels(tick_func(L_ticks))
##ax2t.set_xbound(ax2.get_xbound())
#ax2t.set_xlabel("$\lambda$ (m)")

# %% Error due to parameter uncertainty

Float = E76
wfi = Float.__wfi
data = [getattr(Float, data_name) for data_name in wfi.data_names]
iz, ix = Float.Ww.shape
ip = wfi.ps.shape[0]
w_set = np.empty((iz, ix, ip))

for i, p_set in enumerate(wfi.ps):
    w_set[:, :, i] = wfi.model_func(p_set, data, wfi.fixed)

plt.plot(np.std(w_set, axis=-1), z, color='black', alpha=0.1)
