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
    print("Floats {} and {} exist!".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4976)
    E76.apply_w_model('../../data/EM-APEX/4976_fix_p0k0M_fit_info.p')
    E76.apply_strain('../../data/EM-APEX/4976_N2_ref_300dbar.p')
    E76.apply_isopycnal_displacement('../../data/EM-APEX/srho_4976_100mbin.p')
    E77 = emapex.EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4977)
    E77.apply_w_model('../../data/EM-APEX/4977_fix_p0k0M_fit_info.p')
    E77.apply_strain('../../data/EM-APEX/4977_N2_ref_300dbar.p')
    E77.apply_isopycnal_displacement('../../data/EM-APEX/srho_4977_100mbin.p')


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

    posidxs = detect_peaks(w, mph=0.08, mpd=100.)
    negidxs = detect_peaks(w, mph=0.08, mpd=100., valley=True)
    pidxs = np.hstack((negidxs, posidxs))
    TF = np.arange(len(pidxs))
    tpeaks = t[pidxs]
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
    print("x = {}".format(x[pidxs]))
    print("z = {}".format(z[pidxs]))
    print("t = {}".format(t[pidxs]))

    a = np.transpose(np.vstack((x[pidxs], z[pidxs], t[pidxs])))
    b = np.pi*TF
    x, resid, rank, s = lstsq(a, b)
    print(x)
    print("Horizontal wavelength: {:1.0f} m\nVertical wavelength: {:1.0f} m\n"
          "Period: {:1.0f} min".format(np.pi*2/x[0], np.pi*2/x[1],
                                       np.pi*2/x[2]/60.))

xq = np.array([98650.04863147, 99100.,   99645.62198467,  100959.04961829,  102753.52791723,
  104098.19159536,   98133.38286495,  100127.3972755,   101718.14166566,
  103380.68708482])
zq = np.array([-788.92712687, -1246.54503641, -1573.86360897, -1026.39032746,  -567.31845563,
  -602.50953952, -1376.43262381, -1352.10223918,  -807.29489018])
tq = np.array([6.34613303e+10,   6.34613324e+10,   6.34613352e+10,   6.34613389e+10,
   6.34613417e+10,   6.34613292e+10,   6.34613334e+10,   6.34613367e+10,
   6.34613402e+10])

# %%

a = np.array([[1300., -700., 25*60.],
              [0., -1200., 0.],
              [800., -900., 20*60.],
              [-600., -1400., -20*60.]])

b = np.array([np.pi, np.pi, 0., 0.])
x, resid, rank, s = lstsq(a, b)
print(x)
print("Horizontal wavelength: {:1.0f} m\nVertical wavelength: {:1.0f} m\n"
      "Period: {:1.0f} min".format(np.pi*2/x[0], np.pi*2/x[1],
                                   np.pi*2/x[2]/60.))

#                 G1    B1   G2   B1     B2     G1
w_amp = np.array([-0.2, 0.2, 0.2, -0.23, -0.15, 0.15])
b_amp = np.array([4e-4, 3e-4, 3e-4, 4e-4, 3.5e-4, 2e-4])
u_amp = np.array([0.2, 0.2, 0.2, 0.15, 0.15, 0.1])
v_amp = np.array([0.2, 0.1, 0.1, 0.25, 0.2, 0.15])
N = 0.0017
f = 1.2e-4

om = np.abs(w_amp*N**2/b_amp)

r = np.abs((1j*u_amp/v_amp*f + om)/(u_amp/v_amp*om - 1j*f))

print(np.mean(om))
print(r)

w0 = 0.2
U0 = 0.2
U1 = 0.4
h0 = 500
h1 = 1000

print(2*np.pi/(w0/(U0*h0)), 2*np.pi/(w0/(U1*h1)))