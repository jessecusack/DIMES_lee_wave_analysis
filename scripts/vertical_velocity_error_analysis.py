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

lib_path = os.path.abspath('../../ocean-tools')
if lib_path not in sys.path:
    sys.path.append(lib_path)

import emapex
import plotting_functions as pf
import coloured_noise as cn

try:
    print("Floats {} and {}.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.load(4976)
    E76.generate_regular_grids(dz=1.)
    E77 = emapex.load(4977)
    E77.generate_regular_grids(dz=1.)

# Figure save path.
sdir = os.path.join('..', 'figures', 'vertical_velocity_analysis')
if not os.path.exists(sdir):
    os.makedirs(sdir)
# Universal figure font size.
matplotlib.rc('font', **{'size': 8})


# %% ##########################################################################
# Vertical velocity histograms
t1, w1 = E76.get_timeseries(np.arange(1, 600), 'Ww')
t2, w2 = E77.get_timeseries(np.arange(1, 600), 'Ww')
__, z1 = E76.get_timeseries(np.arange(1, 600), 'z')
__, z2 = E77.get_timeseries(np.arange(1, 600), 'z')

zmin = -100.  # Remove top number of m.
w_combined = 100.*np.hstack((w1[z1 < zmin], w2[z2 < zmin]))
meanw = np.mean(w_combined)
stdw = np.std(w_combined)
dist = sp.stats.norm(loc=meanw, scale=stdw)

print("Mean vertical velocity: {} +/- {}".format(meanw, stdw))

bins = np.arange(-8., 8.1, 0.1)

fig, ax = plt.subplots(1, 1, figsize=(3.125, 3))
n, __, patches = ax.hist(w_combined, bins=bins, normed=True, histtype='step',
                         color='black')
#ax.plot(bins, dist.pdf(bins), linestyle=':', color='black')
ax.set_xlabel('$w$ (cm s$^{-1}$)')
ax.set_xlim(np.min(bins), np.max(bins))
ax.annotate(r"Mean = {:1.1f} $\pm$ {:1.0f} cm s$^{{-1}}$".format(meanw, stdw),
            (-7., 0.4))
ax.set_yticks([])

pf.my_savefig(fig, 'both', 'w_hist', sdir, ftype='pdf', fsize='single_col')

# %% ##########################################################################

pfl = E77.get_profiles(131)
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
dze = 1.
dt = dze/0.13  # 0.12 is a tpyical descent/ascent speed.
zvals = np.arange(-1400., -100, dze)
N = 50
Pe = np.empty((zvals.size/2, N))

# Use an arbitrary profile for the sizes.
zmax = -100.
dz = 1.
pfl = E76.get_profiles(155)
surface = pfl.r_z > zmax
z = pfl.r_z[~surface]
Pw76 = np.empty((z.size/2+1, N))
Pw77 = np.empty((z.size/2+1, N))

for i in xrange(N):
    ze = zvals + 0.15*np.random.rand(zvals.size)
    We = np.diff(ze)/dt
    me, Pe[:, i] = sp.signal.periodogram(We, 1./dze, scaling=scaling)

    # Pick a random profile.
    pfl76 = E76.get_profiles(160+i)
    pfl77 = E77.get_profiles(160+i)

    # Chop out the top.
    Ww76 = pfl76.r_Ww[pfl76.r_z < zmax]
    Ww77 = pfl77.r_Ww[pfl77.r_z < zmax]

    m, Pw76[:, i] = sp.signal.periodogram(Ww76, 1./dz, scaling=scaling)
    __, Pw77[:, i] = sp.signal.periodogram(Ww77, 1./dz, scaling=scaling)

mean_Pe = np.mean(Pe, axis=1)
mean_Pw76 = np.mean(Pw76, axis=1)
mean_Pw77 = np.mean(Pw77, axis=1)

# Sampling every 2.5 m gives 1./5. as highest resolved wavenumber
mmax = 1./5.
use = m < mmax
use[0] = False

# Remove the 0 wavenumber (infinite wavelength)
m, mean_Pw76, mean_Pw77, me, mean_Pe = m[use], mean_Pw76[use], mean_Pw77[use],\
    me[use], mean_Pe[use]

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


def spectrum_func(params, om):
    c, S = params
    return c*om**(-5./3.) + S*om**2


def cost(params, data, om):
    return np.log10(data/spectrum_func(params, om))


# Choose float and range of hpids.
hpids = np.arange(50, 150)

# Chop ends where model is invalid.
zmax = -30.
zmin = -1440.
tmax = 10000.
tmin = 20.
# Interpolation points.
Ni = 2**11
# Choose frequencies to estimate power.
tdmax = 10000.
tdmin = 20.
dzmax = 1200.
dzmin = 4.
N = 100
om = np.logspace(np.log10(np.pi*2/tdmax), np.log10(np.pi*2/tdmin), N)
m = np.logspace(np.log10(np.pi*2/dzmax), np.log10(np.pi*2/dzmin), N)


fig1, ax1 = plt.subplots(1, 1, figsize=(3.125, 3))
#fig2, ax2 = plt.subplots(1, 1, figsize=(3.125, 3))

for j, Float in enumerate([E76, E77]):

    pfls = Float.get_profiles(hpids)

    om_pgrams = np.empty((om.size, pfls.size))
    m_pgrams = np.empty((m.size, pfls.size))

    for i, pfl in enumerate(pfls):
        nans = np.isnan(pfl.Ww)
        w = pfl.Ww[~nans]
        t = pfl.UTC[~nans]*24*60*60
        z = pfl.z[~nans]

        invalid = (z > zmax) | (z < zmin)
        w = w[~invalid]
        t = t[~invalid]
        z = z[~invalid]

        om_pgrams[:, i] = sp.signal.lombscargle(t, w, om)/om
        m_pgrams[:, i] = sp.signal.lombscargle(t, w, m)/m

    f = om/(np.pi*2)
    m = m/(np.pi*2)

    ax1.loglog(1./f, om_pgrams, color='grey', alpha=0.05)
    ax1.loglog(1./f, np.median(om_pgrams, axis=-1), color='black',
               linewidth=3.)



#    ax2.loglog(1./m, m_pgrams, color='grey', alpha=0.05)
#    ax2.loglog(1./m, np.median(m_pgrams, axis=-1), color='black',
#               linewidth=3.)


    use = 1./f < 300.
    data = np.median(om_pgrams, axis=-1)
    params0 = np.array([4e-6, 0.07])
    fit, __ = sp.optimize.leastsq(cost, params0, args=(data[use], om[use]))

    print("c = {:1.2e}".format(fit[0]))
    print("S = {:1.2e}".format(fit[1]))
    print("w variance = {:1.0e} cm/s".format(100.*np.sqrt(fit[1]*om.min())))

    # Only label one line
    if j == 0:
        ax1.loglog(1./f[use], spectrum_func(fit, om[use]), color='red',
                   label=r'$c \omega^{-\frac{5}{3}} + S\omega^2$')
    else:
        ax1.loglog(1./f[use], spectrum_func(fit, om[use]), color='red')

    def tick_func(ticks):
        return ["{:.0f}".format(tick) for tick in ticks]

ax1.set_xticklabels(tick_func(ax1.get_xticks()))
#ax2.set_xticklabels(tick_func(ax2.get_xticks()))

ax1.set_ylabel(r'$w$ power spectra (m$^2$ s$^{-1}$ rad$^{-1}$)')
ax1.set_xlabel('$T$ (s)')
ax1.grid()
ax1.legend(loc=0)

#ax2.set_xlabel('$\lambda$ (m)')
#ax2.grid()
#
#pf.my_savefig(fig1, 'both', 'w_spectrum', sdir, ftype='pdf',
#              fsize='single_col')

######## Using Interpolation #########

# Choose float and range of hpids.
#hpids = np.arange(50, 150)
#Float = E76
#pfls = Float.get_profiles(hpids)
scaling = 'density'
method = 'welch'
nperseg = 512

ti, dt = np.linspace(tmin, tmax, Ni, retstep=True)
zi, dz = np.linspace(zmin, zmax, Ni, retstep=True)

if method == 'welch':
    spec_func = sp.signal.welch
    nom = nperseg/2 + 1
    nm = nperseg/2 + 1
elif method == 'periodogram':
    sepc_func = sp.signal.periodogram
    nom = ti.size/2 + 1
    nm = zi.size/2 + 1

om_pgrams = np.empty((nom, pfls.size))
m_pgrams = np.empty((nm, pfls.size))

for i, pfl in enumerate(pfls):
    nans = np.isnan(pfl.Ww)
    w = pfl.Ww[~nans]
    t = pfl.UTC[~nans]*24*60*60
    z = pfl.z[~nans]

    t1 = t

    # If profile is descenting then z will not be monotonically increasing and
    # the interpolation will fail.
    if np.mean(np.diff(z)) < 0.:
        w = np.flipud(w)
        z = np.flipud(z)

    invalid = (z > zmax) | (z < zmin)
    w = w[~invalid]
    t = t[~invalid]
    z = z[~invalid]

    wit = np.interp(ti, t - t[0], w)
    wiz = np.interp(zi, z, w)

    f, om_pgrams[:, i] = spec_func(wit, fs=1./dt, nperseg=nperseg,
                                   scaling=scaling)
    m, m_pgrams[:, i] = spec_func(wiz, fs=1./dz, nperseg=nperseg,
                                  scaling=scaling)

# Chop interpolation noise.
fuse = f < 1./tdmin
muse = m < 1./dzmin
f, om_pgrams = f[fuse], om_pgrams[fuse, :]
m, m_pgrams = m[muse], m_pgrams[muse, :]
om_pgrams[0, :] = 0.
m_pgrams[0, :] = 0.

# Convert to radian units.
ax1.loglog(1./f, om_pgrams, color='red', alpha=0.1)
ax1.loglog(1./f, np.median(om_pgrams, axis=-1), color='black',
         linewidth=3.)


#ax2.loglog(1./m, m_pgrams, color='red', alpha=0.1)
#ax2.loglog(1./m, np.median(m_pgrams, axis=-1), color='black',
#           linewidth=3.)


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


for Float in [E76, E77]:

    wfi = Float.__wfi
    z = Float.z
    data = [getattr(Float, data_name) for data_name in wfi.data_names]
    iz, ix = Float.Ww.shape
    ip = wfi.ps.shape[0]
    w_set = np.empty((iz, ix, ip))

    print("Float {}".format(Float.floatID))
    for i in xrange(len(wfi.p)):
        pfit = wfi.p[i]
        pstd = np.std(wfi.ps[:, i])
        print("P{}: {} +/- {}".format(i, pfit, pstd))
    print("\n")

    for i, p_set in enumerate(wfi.ps):
        w_set[:, :, i] = wfi.model_func(p_set, data, wfi.fixed)

    plt.figure()
    plt.plot(np.std(w_set, axis=-1), z, color='black', alpha=0.1)

# %% Recreating w profiles with red noise

dz = 2.5
N = 600
beta = -2.  # w variance spectrum in wavenumber has slope close to -2

z = np.arange(0., N*dz, dz)
w = 0.01*cn.noise(N, dz, beta, std=0.5)

m, Pw = sp.signal.periodogram(w, 1./dz)
Pw[0] = 0.

Pyfit = lambda x, a, b: a*x + b
popt, __ = sp.optimize.curve_fit(Pyfit, np.log10(m[1:]), np.log10(Pw[1:]),
                                 p0=[1., 1.])
a, b = popt

fig, axs = plt.subplots(1, 2)
axs[0].loglog(m, Pw, 'k', label=None)
axs[0].loglog(m[1:], m[1:]**a*10**b, 'r',
              label="Fit exponent: {:1.0f}".format(popt[0]))
axs[1].plot(w, z, 'k')

axs[0].set_xlabel('Wavenumber')
axs[0].set_ylabel('Variance')
axs[1].set_xlabel('Velocity')
axs[1].set_ylabel('Depth')

axs[0].legend(loc=0)