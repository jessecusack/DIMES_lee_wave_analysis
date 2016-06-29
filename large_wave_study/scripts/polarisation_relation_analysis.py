# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 16:25:10 2014

@author: jc3e13
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import glob
from scipy.linalg import lstsq
import scipy.optimize as op

import utils
import emapex
from detect_peaks import detect_peaks
import plotting_functions as pf


try:
    print("Floats {} and {} exist!.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.load(4976)
    E77 = emapex.load(4977)

# %% Script params.

# Bathymetry file path.
bf = os.path.abspath(glob.glob('/noc/users/jc3e13/storage/smith_sandwell/topo_*.img')[0])
# Figure save path.
sdir = '../figures/polarisation_relation_analysis'
if not os.path.exists(sdir):
    os.makedirs(sdir)
# Universal figure font size.
matplotlib.rc('font', **{'size': 8})

# %%

E76_hpids = np.array([31, 32])
E77_hpids = np.array([26, 27])

for Float, hpids in zip([E76, E77], [E76_hpids, E77_hpids]):
    __, idxs = Float.get_profiles(hpids, ret_idxs=True)

    z = Float.z[:, idxs]
    zef = Float.zef[:, idxs]
    u = Float.U_abs[:, idxs]
    v = Float.V_abs[:, idxs]
    w = Float.Ww[:, idxs]
    b = Float.b[:, idxs]
    N = np.sqrt(Float.N2_ref[:, idxs])

    ud = utils.nan_detrend(zef, u)
    vd = utils.nan_detrend(zef, v)

    fig, axs = plt.subplots(1, 5, sharey=True, figsize=(16, 6))
    axs[0].set_ylabel('$z$ (m)')
    axs[0].plot(ud, zef)
    axs[0].set_xlabel('$u$')
    axs[1].plot(vd, zef)
    axs[1].set_xlabel('$v$')
    axs[2].plot(w, z)
    axs[2].set_xlabel('$w$')
    axs[3].plot(b, z)
    axs[3].set_xlabel('$b$')
    axs[4].plot(N, z)
    axs[4].set_xlabel('$N$')

    for ax in axs:
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)


# %%

E76_hpids = np.array([31, 32])
E77_hpids = np.array([26, 27])

for Float, hpids in zip([E76, E77], [E76_hpids, E77_hpids]):

    t, x = Float.get_timeseries(hpids, 'dist_ctd')
#    tef, xef = Float.get_timeseries(hpids, 'dist_ef')
#    __, u = Float.get_timeseries(hpids, 'U_abs')
#    __, v = Float.get_timeseries(hpids, 'V_abs')
    __, w = Float.get_timeseries(hpids, 'Ww')
#    __, b = Float.get_timeseries(hpids, 'b')
    __, z = Float.get_timeseries(hpids, 'z')

    posidxs = detect_peaks(w, mph=0.05, mpd=100.)
    negidxs = detect_peaks(w, mph=0.05, mpd=100., valley=True)
    pidxs = np.sort(np.hstack((negidxs, posidxs)))
#    TF = np.arange(len(pidxs))
#    tpeaks = t[pidxs]
#    pidxsef = np.searchsorted(tef, tpeaks)

#    t = utils.datenum_to_datetime(t)
#    tef = utils.datenum_to_datetime(tef)
    t *= 86400.
#    tef *= 86400.
    x *= 1000.
#    xef *= 1000.

    fig, axs = plt.subplots(3, 1, figsize=(12, 10))

    axs[0].plot(x, w)
    axs[0].plot(x[pidxs], w[pidxs], 'ro')
    axs[0].set_xlabel('$x$')
    axs[1].plot(t, w)
    axs[1].plot(t[pidxs], w[pidxs], 'ro')
    axs[1].set_xlabel('$t$')
    axs[2].plot(z, w)
    axs[2].plot(z[pidxs], w[pidxs], 'ro')
    axs[2].set_xlabel('$z$')

    for ax in axs.flatten():
        ax.grid()
        ax.set_ylabel('$w$')

    print("Float {}, peak coordinates:".format(Float.floatID))
    print("x = {}".format(x[pidxs] - x[pidxs][0]))
    print("z = {}".format(z[pidxs]))
    print("t = {}".format(t[pidxs] - t[pidxs][0]))
    print("u = ")
    print("v = ")
    print("w = ")
    print("b = ")

#    a = np.transpose(np.vstack((x[pidxs], z[pidxs], t[pidxs])))
#    b = np.pi*TF
#    x, resid, rank, s = lstsq(a, b)
#    print(x)
#    print("Horizontal wavelength: {:1.0f} m\nVertical wavelength: {:1.0f} m\n"
#          "Period: {:1.0f} min".format(np.pi*2/x[0], np.pi*2/x[1],
#                                       np.pi*2/x[2]/60.))

# %% NOTE: THIS METHOD IS OBSOLETE AND PROVEN NOT TO WORK
# Estimate position (x, z, t) of wave peaks from combination of profiles 31 and
# 32 from float 4976.
xq = np.array([0.,  516.66576652, 882., 1512.23911972,  1994.01441055,  2825.66675334,
  3584.75880071,  4620.14505228,  5247.30421987,  5964.80873041])
zq = np.array([-602.50953952,  -788.92712687, -1000., -1246.54503641, -1376.43262381, -1573.86360897,
 -1352.10223918, -1026.39032746, -807.29489018, -567.31845563,])
tq = np.array([0., 1081., 2000., 3164.00000763, 4172., 5911.99999237,
   7500., 9666., 10978., 12479.])

xqd = np.diff(xq)
zqd = np.diff(zq)
tqd = -np.diff(tq)
#U = 0.3
#xqd -= U*tqd
pis = np.pi*np.ones_like(tqd)

a = np.transpose(np.vstack((xqd, zqd, tqd)))
X, resid, rank, s = lstsq(a, pis)

# See how susceptible to random noise in X this result is.
Neq = len(a[:, 0])
N = 10000
X_dist = []
for i in xrange(N):
    b = a.copy()
    b[:, 0] += 100.*np.random.randn(Neq)
    Xt, __, __, __ = lstsq(b, pis)
    X_dist.append(Xt)

X_dist = np.pi*2./np.asarray(X_dist)
# Convert to minutes.
X_dist[:, 2] /= 60.
X_mean = np.mean(X_dist, axis=0)
X_std = np.std(X_dist, axis=0)

fig, axs = plt.subplots(1, 3, figsize=(10, 6))

axs[0].hist(X_dist[:,0], bins=np.arange(-20000, 20000, 500),
            normed=True, log=False, alpha=0.8);
axs[1].hist(X_dist[:,1], bins=np.arange(-20000, 20000, 500),
            normed=True, log=False, alpha=0.8);
axs[2].hist(X_dist[:,2], bins=np.arange(-300, 300, 10),
            normed=True, log=False, alpha=0.8, label='4976');


print(X)
print("Horizontal wavelength: {:1.0f} m\nVertical wavelength: {:1.0f} m\n"
      "Period: {:1.0f} min".format(np.pi*2/X[0], np.pi*2/X[1],
                                   np.pi*2/X[2]/60.))
print(X_mean)
print(X_std)
print(np.percentile(X_dist, [5, 25, 50, 75, 95] , axis=0))

# Estimate position (x, z, t) of wave peaks from combination of profiles 26 and
# 27 from float 4977.
xq = [0., 637.03576161, 1297.34685781, 1879.7361852, 5746.97095521, 6452.82920205, 6899.69383473, 7464.2015048]
zq = [-1.29372258e+03, -1.10985553e+03, -8.92841297e+02, -7.02023797e+02, -3.41960178e+02, -5.52210564e+02, -7.03930538e+02, -9.38402544e+02]
tq = [0., 1259., 2564., 3715.00000763, 11803., 13381., 14380., 15642.]
#step = [0., 0., 0., 0., 1., 1., 1.,]

xqd = np.diff(xq)
zqd = np.diff(zq)
tqd = -np.diff(tq)
#stepd = np.diff(step)
#U = 0.3
#xqd -= U*tqd
pis = np.pi*np.ones_like(tqd)

## Big step at surface, I think we missed 6 phase lines.
#pis[3] *= 6.
use = np.arange(len(pis)) != 3

a = np.transpose(np.vstack((xqd[use], zqd[use], tqd[use])))
pis = pis[use]
X, resid, rank, s = lstsq(a, pis)

# See how susceptible to random noise in X this result is.
Neq = len(a[:, 0])
N = 10000
X_dist = []
for i in xrange(N):
    b = a.copy()
    b[:, 0] += 100.*np.random.randn(Neq)
    Xt, __, __, __ = lstsq(b, pis)
    X_dist.append(Xt)

X_dist = np.pi*2./np.asarray(X_dist)
# Convert to minutes.
X_dist[:, 2] /= 60.
X_mean = np.mean(X_dist, axis=0)
X_std = np.std(X_dist, axis=0)

axs[0].hist(X_dist[:,0], bins=np.arange(-20000, 20000, 500),
            normed=True, log=False, alpha=0.4);
axs[1].hist(X_dist[:,1], bins=np.arange(-20000, 20000, 500),
            normed=True, log=False, alpha=0.4);
axs[2].hist(X_dist[:,2], bins=np.arange(-300, 300, 10),
            normed=True, log=False, alpha=0.4, label='4977');

axs[2].legend()

for ax in axs:
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=60)
    ax.set_yticks([])
    ax.grid()

axs[0].set_ylabel('Occurances (-)')
axs[0].set_xlabel('$\lambda_x$ (m)')
axs[1].set_xlabel('$\lambda_z$ (m)')
axs[2].set_xlabel('$T$ (min)')

pf.my_savefig(fig, 'both_matrix_inversion', 'hist', sdir, fsize='double_col')

print(X)
print("Horizontal wavelength: {:1.0f} m\nVertical wavelength: {:1.0f} m\n"
      "Period: {:1.0f} min".format(np.pi*2/X[0], np.pi*2/X[1],
                                   np.pi*2/X[2]/60.))
print(X_mean)
print(X_std)
print(np.percentile(X_dist, [5, 25, 50, 75, 95], axis=0))
# %%

# These are the measured amplitudes after detrending a 2nd order polynomial.
# Always bottom up.

#                    u,    v,    w, b
amps_31 = np.array([[0.17, 0.11, 0.17, 4e-4],
                    [-0.12, -0.08, -0.025, -4e-4],
                    [0.2, 0.045, 0.07, -3.6e-4]])
amps_32 = np.array([[0.1, 0.18, 0.19, 2e-4],
                    [-0.06, -0.2, -0.2, -5e-4],
                    [-0.05, 0.22, 0.1, 5e-4],
                    [-0.17, 0., -0.2, -6e-4]])
amps_26 = np.array([[-0.095, -0.075, 0.16, 4.5e-4],
                    [0.2, 0.16, -0.1, -4e-4],
                    [-0.1, -0.04, 0.17, 3e-4],
                    [0.11, 0.17, -0.21, -6e-4]])
amps_27 = np.array([[-0.13, 0., -0.08, -3e-4],
                    [0.1, 0., 0.12, 2.7e-4],
                    [-0.03, 0., -0.07, 2.3e-4]])

amps = np.vstack((amps_31, amps_32, amps_26, amps_27))

#                 G1    B1   G2   B1     B2     G1
#w_amp = np.array([-0.2, 0.2, 0.2, -0.23, -0.15, 0.15])
#b_amp = np.array([4e-4, 3e-4, 3e-4, 4e-4, 3.5e-4, 2e-4])
#u_amp = np.array([0.2, 0.2, 0.2, 0.15, 0.15, 0.1])
#v_amp = np.array([0.2, 0.1, 0.1, 0.25, 0.2, 0.15])

u_amp = amps[:, 0]
v_amp = amps[:, 1]
w_amp = amps[:, 2]
b_amp = amps[:, 3]

N = 2.2e-3
f = 1.2e-4

om = np.abs(w_amp*N**2/b_amp)

alpha2 = w_amp**2/(u_amp**2 + v_amp**2)
alpha = np.sqrt(alpha2)
om_alt = np.sqrt(alpha2*N**2/(1 + alpha2))

# With rotation.
#r1 = np.abs((1j*u_amp/v_amp*f + om)/(u_amp/v_amp*om - 1j*f))
# Without rotation.
#r2 = np.abs(v_amp/u_amp)

print("omega = {:1.2e} +/- {:1.2e}".format(np.mean(om_alt), np.std(om_alt)))
print("omega (buoyancy) = {:1.2e} +/- {:1.2e}".format(np.mean(om), np.std(om)))
print("alpha = {:1.1f} +/- {:1.1f}".format(np.mean(alpha), np.std(alpha)))

fig = plt.figure(figsize=(3.125, 3))

gs = gridspec.GridSpec(1, 2, width_ratios=[2,1])

axs = [plt.subplot(gs[0]), plt.subplot(gs[1])]

axs[1].yaxis.tick_right()
axs[1].yaxis.set_ticks_position('both')
axs[1].yaxis.set_label_position('right')

labels=[r'$\frac{\alpha^2}{1 + \alpha^2}N^2$', r'$|\frac{w_0}{b_0}|N^2$']
b = axs[0].boxplot([om_alt/N, om/N], labels=labels, showfliers=False);
[plt.setp(b[key], color='k') for key in b.keys()]
axs[0].set_ylabel('$\omega/N$ (-)')
axs[0].hlines(f/N, *axs[0].get_xlim())
axs[0].annotate('$f$', xy=(0.6, 1.1*f/N))
axs[0].hlines(1., *axs[0].get_xlim())
axs[0].annotate('$N$', xy=(0.6, 0.76))

b = axs[1].boxplot(alpha, labels=[r'$\frac{w_0^2}{u_0^2 + v_0^2}$'], showfliers=False);
[plt.setp(b[key], color='k') for key in b.keys()]
axs[1].set_ylabel(r'$\alpha$ (-)')

for ax in axs:
    ax.set_yscale('log')
    plt.setp(ax.get_xticklabels(), fontsize=12)

axs[0].set_yticks([1e-2, 1e-1, 1e0, 1e1])
axs[0].set_yticklabels(['0.01', '0.1', '1.0', '10'])
axs[1].set_yticks([1e-1, 1e0, 1e1])
axs[1].set_yticklabels(['0.1', '1.0', '10'])

fontdict = {'size':10}
plt.figtext(0, 0.9, 'a)', fontdict=fontdict)
plt.figtext(0.61, 0.9, 'b)', fontdict=fontdict)

pf.my_savefig(fig, 'both', 'omega_alpha_box', sdir, ftype=('png', 'pdf'),
              fsize='single_col')


# %%

# New approach that constrains the frequency using the wavenumber information.

U = 0.4
N = 2e-3

xq = np.array([0.,  516.66576652, 882., 1512.23911972,  1994.01441055,  2825.66675334,
  3584.75880071,  4620.14505228,  5247.30421987,  5964.80873041])
zq = np.array([-602.50953952,  -788.92712687, -1000., -1246.54503641, -1376.43262381, -1573.86360897,
 -1352.10223918, -1026.39032746, -807.29489018, -567.31845563,])
tq = np.array([0., 1081., 2000., 3164.00000763, 4172., 5911.99999237,
   7500., 9666., 10978., 12479.])

xqd = np.diff(xq)
zqd = np.diff(zq)
tqd = np.diff(tq)

def lin_model(params, dx, dz, dt, N, U):
    k, m = params
    om0 = np.sqrt(k**2 * N**2 /(k**2 + m**2))
    y = k*dx + m*dz - (om0 + U*k)*dt
    return y - np.pi

args = (xqd, zqd, tqd, N, U)

popt, __, info, __, __ = op.leastsq(lin_model, x0=[-0.0012, -0.0012],
                                    args=args, full_output=True)