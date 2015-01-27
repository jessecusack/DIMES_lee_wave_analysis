# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 16:25:10 2014

@author: jc3e13
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import sys
import glob
from scipy.linalg import lstsq
import gsw

lib_path = os.path.abspath('../python')
if lib_path not in sys.path:
    sys.path.append(lib_path)

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
bf = os.path.abspath(glob.glob('../../data/sandwell_bathymetry/topo_*.img')[0])
# Figure save path.
sdir = '../figures/polarisation_relation_analysis'
if not os.path.exists(sdir):
    os.makedirs(sdir)
# Universal figure font size.
matplotlib.rc('font', **{'size': 9})

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

# %% CRAZY MEAN MINUS plot for pfl 26

#N = 10
#N0_76 = 15
#N0_77 = 10
#E76_hpids = np.arange(N0_76, N0_76+N)
#E77_hpids = np.arange(N0_77, N0_77+N)
#
#dz = 1.
#z = np.arange(-1500, 0, dz)
#rho = []
#
#pfls = np.hstack((E76.get_profiles(E76_hpids), E77.get_profiles(E77_hpids)))
#
#for pfl in pfls:
#    rho.append(pfl.interp(z, 'z', 'rho_1'))
#
#rho = np.transpose(np.asarray(rho))
#mrho = np.mean(rho, axis=-1)

#axs[0].plot(mrho, z, 'red')
#axs[1].plot(mrho, z, 'red')

pfl = E77.get_profiles(26)
#srhop = utils.nan_interp(pfl.z, z, mrho)
#b = -gsw.grav(pfl.lat_start, pfl.P)*(pfl.rho_1 - srhop)/1031.

zmax = -650
use = pfl.z < zmax
useef = pfl.zef < zmax

fig, axs = plt.subplots(1, 4, sharey=True)
axs[0].set_ylabel('$z$ (m)')
#axs[0].plot(pfl.b, pfl.z, 'grey')
#axs[0].plot(b, pfl.z, 'red')
axs[0].plot(utils.nan_detrend(pfl.zef[useef], pfl.U_abs[useef]), pfl.zef[useef], 'red')
axs[0].set_xlabel('$U$ (m s$^{-1}$)')
plt.setp(axs[0].xaxis.get_majorticklabels(), rotation=45)
axs[1].plot(utils.nan_detrend(pfl.zef[useef], pfl.V_abs[useef]), pfl.zef[useef], 'red')
axs[1].set_xlabel('$V$ (m s$^{-1}$)')
plt.setp(axs[1].xaxis.get_majorticklabels(), rotation=45)
axs[2].plot(pfl.Ww[use], pfl.z[use], color='red')
axs[2].set_xlabel('$W$ (m s$^{-1}$)')
plt.setp(axs[2].xaxis.get_majorticklabels(), rotation=45)
axs[3].plot(pfl.b[use], pfl.z[use], color='red')
axs[3].set_xlabel('$b$ (m s$^{-2}$)')
plt.setp(axs[3].xaxis.get_majorticklabels(), rotation=45)

#pf.my_savefig(fig, '4977', 'UVWB', sdir, fsize='double_col')

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

#    a = np.transpose(np.vstack((x[pidxs], z[pidxs], t[pidxs])))
#    b = np.pi*TF
#    x, resid, rank, s = lstsq(a, b)
#    print(x)
#    print("Horizontal wavelength: {:1.0f} m\nVertical wavelength: {:1.0f} m\n"
#          "Period: {:1.0f} min".format(np.pi*2/x[0], np.pi*2/x[1],
#                                       np.pi*2/x[2]/60.))

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
pis = np.pi*np.ones_like(tqd)

a = np.transpose(np.vstack((xqd, zqd, tqd)))
X, resid, rank, s = lstsq(a, pis)

print(X)
print("Horizontal wavelength: {:1.0f} m\nVertical wavelength: {:1.0f} m\n"
      "Period: {:1.0f} min".format(np.pi*2/X[0], np.pi*2/X[1],
                                   np.pi*2/X[2]/60.))

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
pis = np.pi*np.ones_like(tqd)

# Big step at surface, I think we missed 6 phase lines.
pis[3] *= 6.

a = np.transpose(np.vstack((xqd, zqd, tqd)))
X, resid, rank, s = lstsq(a, pis)

print(X)
print("Horizontal wavelength: {:1.0f} m\nVertical wavelength: {:1.0f} m\n"
      "Period: {:1.0f} min".format(np.pi*2/X[0], np.pi*2/X[1],
                                   np.pi*2/X[2]/60.))

# %%


#                 G1    B1   G2   B1     B2     G1
w_amp = np.array([-0.2, 0.2, 0.2, -0.23, -0.15, 0.15])
b_amp = np.array([4e-4, 3e-4, 3e-4, 4e-4, 3.5e-4, 2e-4])
u_amp = np.array([0.2, 0.2, 0.2, 0.15, 0.15, 0.1])
v_amp = np.array([0.2, 0.1, 0.1, 0.25, 0.2, 0.15])
N = 1e-3
f = 1.2e-4

om = np.abs(w_amp*N**2/b_amp)

# With rotation.
r1 = np.abs((1j*u_amp/v_amp*f + om)/(u_amp/v_amp*om - 1j*f))
# Without rotation.
r2 = np.abs(v_amp/u_amp)

print("{:1.2e} +/- {:1.2e}".format(np.mean(om), np.std(om)))
print(r1)
print(r2)
print(np.mean(r2))

w0 = 0.2
U0 = 0.2
U1 = 0.5
h0 = 200
h1 = 1200

print(2*np.pi/(w0/(U0*h0)), 2*np.pi/(w0/(U1*h1)))