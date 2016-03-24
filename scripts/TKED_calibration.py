# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 18:33:56 2016

@author: jc3e13
"""

import sys
import os
import glob
import scipy as sp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#from matplotlib import gridspec
#import scipy.signal as sig
#from scipy.integrate import trapz

import gsw

lib_path = os.path.abspath('../modules')
if lib_path not in sys.path:
    sys.path.append(lib_path)

lib_path = os.path.abspath('../../ocean-tools')
if lib_path not in sys.path:
    sys.path.append(lib_path)

import emapex
import TKED_parameterisations as fs
import plotting_functions as pf
#import sandwell
#import window as wdw

zmin = -1500.
zmax = -100.
dz = 1.

try:
    print("Floats {} and {} exist!.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.load(4976)
    E76.generate_regular_grids(zmin=zmin, dz=dz)
    E77 = emapex.load(4977)
    E77.generate_regular_grids(zmin=zmin, dz=dz)

# %% Script params.

# Bathymetry file path.
bf = os.path.abspath(glob.glob('/noc/users/jc3e13/storage/smith_sandwell/topo_*.img')[0])
# Figure save path.
sdir = '../figures/TKED_estimation'
if not os.path.exists(sdir):
    os.makedirs(sdir)
# Universal figure font size.
matplotlib.rc('font', **{'size': 8})

# Load VMP data.
UK2_vmp = sp.io.loadmat('../../storage/DIMES/combined_jc054.mat',
                        variable_names=['vmp'])['vmp']


# %% Data map

z_vmp = gsw.z_from_p(UK2_vmp['press'][0][0], UK2_vmp['startlat'][0][0][0])

fig1 = plt.figure()
plt.plot(UK2_vmp['startlon'][0][0][0], UK2_vmp['startlat'][0][0][0], 'ko')
plt.plot(UK2_vmp['startlon'][0][0][0][25:30], UK2_vmp['startlat'][0][0][0][25:30], 'ro')
plt.plot(E76.lon_start[:100], E76.lat_start[:100])
plt.plot(E77.lon_start[:100], E77.lat_start[:100])

# %% Coefficient estimation using microstructure data

def w_scales_float(Float, hpids, xvar, dx=1., width=10., lc=30., c=1., eff=0.2,
                   btype='highpass', we=1e-3, ret_noise=False):

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)

    w = Float.r_Ww[:, idxs]
    N2 = Float.r_N2_ref[:, idxs]

    if xvar == 'time':
        x = np.arange(0., 10000., dx)
        __, __, w = Float.get_interp_grid(hpids, x, 'dUTC', 'Ww')
        __, __, N2 = Float.get_interp_grid(hpids, x, 'dUTC', 'N2_ref')
    elif xvar == 'height':
        x = Float.r_z[:, 0]
        w = Float.r_Ww[:, idxs]
        N2 = Float.r_N2_ref[:, idxs]
    else:
        raise ValueError("xvar must either be 'time' or 'height'.")

    epsilon = np.zeros_like(w)
    kappa = np.zeros_like(w)

    if ret_noise:
        epsilon_noise = np.zeros_like(w)
        noise_flag = np.zeros_like(w, dtype=bool)
        for i, (w_row, N2_row) in enumerate(zip(w.T, N2.T)):
            epsilon[:, i], kappa[:, i], epsilon_noise[:, i], noise_flag[:, i] = \
                fs.w_scales(w_row, x, N2_row, dx, width, lc, c, eff, btype, we,
                            ret_noise)
        return epsilon, kappa, epsilon_noise, noise_flag
    else:
        for i, (w_row, N2_row) in enumerate(zip(w.T, N2.T)):
            epsilon[:, i], kappa[:, i] = \
                fs.w_scales(w_row, x, N2_row, dx, width, lc, c, eff, btype)
        return epsilon, kappa

# %% Height

z = np.arange(zmin, 0., dz)
xvar = 'height'
dx = 1.
hpids = np.arange(50, 150)
width = 15.
lc = np.array([40., 15.])
c = 1.
btype = 'bandpass'
we = 0.001

epsilon_76, __, ep_noise_76, flag_76 = \
    w_scales_float(E76, hpids, xvar, dx=dx, width=width, lc=lc, c=c, btype=btype,
                   we=we, ret_noise=True)
epsilon_77, __, ep_noise_77, flag_77 = \
    w_scales_float(E77, hpids, xvar, dx=dx, width=width, lc=lc, c=c, btype=btype,
                   we=we, ret_noise=True)

epsilon_vmp = UK2_vmp['eps'][0][0][:, 25:30].flatten()
z_vmp = gsw.z_from_p(UK2_vmp['press'][0][0], UK2_vmp['startlat'][0][0][0])
z_vmp_flat = z_vmp[:, 25:30].flatten()
use = ~np.isnan(epsilon_vmp) & (z_vmp_flat > zmin) & (z_vmp_flat < zmax)
use_76 = np.array([list(z < zmax),]*len(hpids)).transpose()
use_77 = use_76

c_76 = np.median(epsilon_vmp[use])/np.median(epsilon_76[use_76])
c_77 = np.median(epsilon_vmp[use])/np.median(epsilon_77[use_77])

print("c(4976) = {:1.3f} and c(4977) = {:1.3f}.".format(c_76, c_77))

epsilon_76 *= c_76
epsilon_77 *= c_77
ep_noise_76 *= c_76
ep_noise_77 *= c_77

fig2 = plt.figure()
plt.semilogx(epsilon_76, z, color='red', alpha=0.2)
plt.semilogx(epsilon_77, z, color='red', alpha=0.2)
plt.semilogx(UK2_vmp['eps'][0][0][:, 25:30], z_vmp[:, 25:30], color='grey',
             alpha=0.5)
plt.semilogx(ep_noise_76, z, color='red')
plt.semilogx(ep_noise_77, z, color='red')
plt.ylim(-1500., 0.)

bins = np.arange(-12., -5, 0.25)

fig3, axs = plt.subplots(3, 1, sharex='col', figsize=(3, 6))
axs[0].hist(np.log10(epsilon_vmp[use]), bins=bins, color='blue',
            alpha=0.8, label='VMP')
axs[1].hist(np.log10(epsilon_76[use_76]).flatten(), bins=bins, color='red', alpha=0.8,
            label='4976')
axs[1].hist(np.log10(epsilon_76[~(flag_76 | ~use_76)]), bins=bins,
            color='green', alpha=0.8, label='4976 above noise')
axs[2].hist(np.log10(epsilon_77[use_77]).flatten(), bins=bins, color='red',
            alpha=0.8, label='4977')
axs[2].hist(np.log10(epsilon_77[~(flag_77 | ~use_77)]), bins=bins,
            color='green', alpha=0.8, label='4977 above noise')
axs[2].set_xlabel('$\log_{10}(\epsilon)$ W kg$^{-1}$')

axs[1].text(-7, 0.6, r'c = {:1.3f}'.format(c_76))
axs[2].text(-7, 0.6, r'c = {:1.3f}'.format(c_77))

for ax in axs:
    ax.legend()

# %% Time

xvar = 'time'
dx = 5.
hpids = np.arange(50, 150)
width = 120.
lc = np.array([300, 120])
c = 1.
btype = 'bandpass'
we = 0.001

t = np.arange(0., 10000., dx)
__, __, z_76 = E76.get_interp_grid(hpids, t, 'dUTC', 'z')
__, __, z_77 = E77.get_interp_grid(hpids, t, 'dUTC', 'z')

epsilon_76, __, ep_noise_76, flag_76 = \
    w_scales_float(E76, hpids, xvar, dx=dx, width=width, lc=lc, c=c, btype=btype,
                   we=we, ret_noise=True)
epsilon_77, __, ep_noise_77, flag_77 = \
    w_scales_float(E77, hpids, xvar, dx=dx, width=width, lc=lc, c=c, btype=btype,
                   we=we, ret_noise=True)

epsilon_vmp = UK2_vmp['eps'][0][0][:, 25:30].flatten()
z_vmp_flat = z_vmp[:, 25:30].flatten()
use = ~np.isnan(epsilon_vmp) & (z_vmp_flat > zmin) & (z_vmp_flat < zmax)
use_76 = z_76 < zmax
use_77 = z_77 < zmax

c_76 = np.median(epsilon_vmp[use])/np.median(epsilon_76[use_76])
c_77 = np.median(epsilon_vmp[use])/np.median(epsilon_77[use_77])

print("c(4976) = {:1.3f} and c(4977) = {:1.3f}.".format(c_76, c_77))

epsilon_76 *= c_76
epsilon_77 *= c_77
ep_noise_76 *= c_76
ep_noise_77 *= c_77

fig2 = plt.figure()
plt.semilogx(epsilon_76, z_76, color='red', alpha=0.2)
plt.semilogx(epsilon_77, z_77, color='red', alpha=0.2)
plt.semilogx(UK2_vmp['eps'][0][0][:, 25:30], z_vmp[:, 25:30], color='grey',
             alpha=0.5)
plt.semilogx(ep_noise_76, z_76, color='red')
plt.semilogx(ep_noise_77, z_77, color='red')
plt.ylim(-1500., 0.)

bins = np.arange(-12., -5, 0.25)

fig3, axs = plt.subplots(3, 1, sharex='col', figsize=(3, 6))
axs[0].hist(np.log10(epsilon_vmp[use]), bins=bins, color='blue', alpha=0.8,
            label='VMP')
axs[1].hist(np.log10(epsilon_76[use_76]).flatten(), bins=bins, color='red',
            alpha=0.8, label='4976')
axs[1].hist(np.log10(epsilon_76[~(flag_76 | ~use_76)]), bins=bins,
            color='green', alpha=0.8, label='4976 above noise')
axs[2].hist(np.log10(epsilon_77[use_77]).flatten(), bins=bins, color='red',
            alpha=0.8, label='4977')
axs[2].hist(np.log10(epsilon_77[~(flag_77 | ~use_77)]), bins=bins,
            color='green', alpha=0.8, label='4977 above noise')
axs[2].set_xlabel('$\log_{10}(\epsilon)$ W kg$^{-1}$')

axs[1].text(-7, 0.6, r'c = {:1.3f}'.format(c_76))
axs[2].text(-7, 0.6, r'c = {:1.3f}'.format(c_77))

for ax in axs:
    ax.legend()