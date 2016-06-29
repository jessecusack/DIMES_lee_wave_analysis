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
from matplotlib import gridspec
from scipy.integrate import trapz

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
import sandwell
#import window as wdw

zmin = -1450.
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
bf = os.path.abspath(glob.glob('../../storage/smith_sandwell/topo_*.img')[0])
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
# Height

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
    fs.w_scales_float(E76, hpids, xvar, z, width=width, lc=lc, c=c,
                      btype=btype, we=we, ret_noise=True)
epsilon_77, __, ep_noise_77, flag_77 = \
    fs.w_scales_float(E77, hpids, xvar, z, width=width, lc=lc, c=c,
                      btype=btype, we=we, ret_noise=True)

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

# %% Effectice height

z = np.arange(zmin, 0., dz)
xvar = 'eheight'
dx = 1.
hpids = np.arange(50, 150)
width = 15.
lc = np.array([40., 15.])
c = 1.
btype = 'bandpass'
we = 0.001

epsilon_76, __, ep_noise_76, flag_76 = \
    fs.w_scales_float(E76, hpids, xvar, z, width=width, lc=lc, c=c,
                      btype=btype, we=we, ret_noise=True)
epsilon_77, __, ep_noise_77, flag_77 = \
    fs.w_scales_float(E77, hpids, xvar, z, width=width, lc=lc, c=c,
                      btype=btype, we=we, ret_noise=True)

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

pf.my_savefig(fig3, 'eheight', 'lem_hist', sdir, ftype='png')

# %% Time

xvar = 'time'
dx = 5.
hpids = np.arange(50, 150)
width = 120.
lc = np.array([300, 120])
c = 1.
btype = 'bandpass'
we = 0.001

t = np.arange(0., 11000., dx)
__, __, z_76 = E76.get_interp_grid(hpids, t, 'dUTC', 'z')
__, __, z_77 = E77.get_interp_grid(hpids, t, 'dUTC', 'z')

epsilon_76, __, ep_noise_76, flag_76 = \
    fs.w_scales_float(E76, hpids, xvar, t, width=width, lc=lc, c=c,
                      btype=btype, we=we, ret_noise=True)
epsilon_77, __, ep_noise_77, flag_77 = \
    fs.w_scales_float(E77, hpids, xvar, t, width=width, lc=lc, c=c,
                      btype=btype, we=we, ret_noise=True)

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

# %% timeheight

z = np.arange(zmin, 0., dz)
xvar = 'timeheight'
dx = 1.
hpids = np.arange(50, 150)
width = 15.
lc = (100., 40.)
c = 1.
btype = 'highpass'
we = 0.001

epsilon_76, __, ep_noise_76, flag_76 = \
    fs.w_scales_float(E76, hpids, xvar, z, width=width, lc=lc, c=c,
                      btype=btype, we=we, ret_noise=True)
epsilon_77, __, ep_noise_77, flag_77 = \
    fs.w_scales_float(E77, hpids, xvar, z, width=width, lc=lc, c=c,
                      btype=btype, we=we, ret_noise=True)

vmp_pfls = slice(24, 29)
#vmp_pfls = slice(1, 50)
epsilon_vmp = UK2_vmp['eps'][0][0][:, vmp_pfls].flatten()
z_vmp = gsw.z_from_p(UK2_vmp['press'][0][0], UK2_vmp['startlat'][0][0][0])
z_vmp_flat = z_vmp[:, vmp_pfls].flatten()
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
plt.semilogx(UK2_vmp['eps'][0][0][:, vmp_pfls], z_vmp[:, vmp_pfls],
             color='grey', alpha=0.5)
plt.semilogx(ep_noise_76, z, color='red')
plt.semilogx(ep_noise_77, z, color='red')
plt.ylim(-1500., 0.)

pf.my_savefig(fig2, 'both', 'lem_vmp_pfls', sdir, ftype='png')

bins = np.arange(-12., -5, 0.25)

fig3, axs = plt.subplots(3, 1, sharex='col', figsize=(3, 6))
axs[0].hist(np.log10(epsilon_vmp[use]), bins=bins, color='blue',
            alpha=0.8, label='VMP')
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

pf.my_savefig(fig3, 'timeheight', 'lem_hist', sdir, ftype='png')

# %% Time and effective height.

z = np.arange(zmin, 0., dz)
xvar = 'timeeheight'
dx = 1.
hpids = np.arange(50, 150)
width = 15.
lc = (100., 40.)
c = 1.
btype = 'highpass'
we = 0.001

epsilon_76, __, ep_noise_76, flag_76 = \
    fs.w_scales_float(E76, hpids, xvar, z, width=width, lc=lc, c=c,
                      btype=btype, we=we, ret_noise=True)
epsilon_77, __, ep_noise_77, flag_77 = \
    fs.w_scales_float(E77, hpids, xvar, z, width=width, lc=lc, c=c,
                      btype=btype, we=we, ret_noise=True)

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

pf.my_savefig(fig3, 'timeeheight', 'lem_hist', sdir, ftype='png')

# %% Full trajectory of dissipation.

#cs = [0.197, 0.158]  # time
#xvar = 'time'
#dx = 5.
#width = 120.
#lc = np.array([300., 120.])

#cs = [0.193, 0.160]  # height
#xvar = 'height'
#dx = 1.
#width = 15.
#lc = np.array([40., 15.])

#cs = [0.176, 0.147] # timeheight
#x = np.arange(zmin, 0., dz)
#xvar = 'timeheight'
#dx = 1.
#hpids = np.arange(10, 500)
#width = 15.
#lc = (100., 40.)
#btype = 'highpass'
#we = 0.001

#cs = [0.192, 0.159]  # eheight
#x = np.arange(zmin, 0., dz)
#xvar = 'eheight'
#dx = 1.
#width = 20.
#lc = np.array([40., 15.])

z = np.arange(zmin, 0., dz)  # timeeheight
x = np.arange(zmin, 0., dz)
xvar = 'timeeheight'
dx = 1.
hpids = np.arange(50, 150)
width = 20.
lc = (100., 40.)
c = 1.
btype = 'highpass'
we = 0.001

epsilon_76, __, ep_noise_76, flag_76 = \
    fs.w_scales_float(E76, hpids, xvar, z, width=width, lc=lc, c=c,
                      btype=btype, we=we, ret_noise=True)
epsilon_77, __, ep_noise_77, flag_77 = \
    fs.w_scales_float(E77, hpids, xvar, z, width=width, lc=lc, c=c,
                      btype=btype, we=we, ret_noise=True)

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

pf.my_savefig(fig3, xvar, 'lem_hist', sdir, ftype='png')

####

Float = E77
c = c_77
hpids = np.arange(10, 70)

fig = plt.figure(figsize=(7, 4))
gs = gridspec.GridSpec(2, 1, height_ratios=[1, 5])
ax0 = plt.subplot(gs[1])
ax1 = plt.subplot(gs[0])

__, idxs = Float.get_profiles(hpids, ret_idxs=True)

epsilon, kappa = fs.w_scales_float(Float, hpids, xvar, x, width=width,
                                   lc=lc, c=c, btype=btype, we=we,
                                   ret_noise=False)

ieps = 0.*np.zeros_like(idxs)

if xvar == 'time':
    t = np.arange(0., 10000., dx)
    __, __, iZ = Float.get_interp_grid(hpids, t, 'dUTC', 'z')
    __, __, X = Float.get_interp_grid(hpids, t, 'dUTC', 'dist_ctd')
if xvar == 'eheight' or xvar == 'timeeheight':
    __, __, it = Float.get_interp_grid(hpids, x, 'zw', 'dUTC')
    iZ = np.zeros_like(it)
    X = np.zeros_like(it)
    for i, pfl in enumerate(Float.get_profiles(hpids)):
        iZ[:, i] = pfl.interp(it[:, i], 'dUTC', 'z')
        X[:, i] = pfl.interp(it[:, i], 'dUTC', 'dist_ctd')
elif xvar == 'height' or xvar == 'timeheight':
    __, __, iZ = Float.get_interp_grid(hpids, x, 'z', 'z')
    __, __, X = Float.get_interp_grid(hpids, x, 'z', 'dist_ctd')

for i in xrange(len(idxs)):
    iuse = (iZ[:, i] < -100) & (iZ[:, i] > -1400)
    # The abs function accounts for problems with z being the wrong way.
    ieps[i] = np.abs(1025.*trapz(epsilon[iuse, i], iZ[iuse, i]))

print("Max integrated dissipation: {}".format(np.max(ieps)))

Z = iZ.flatten()

use = (Z < -100) & (Z > -1400)

Z = Z[use]

X = X.flatten()[use]
LOG_EPS = (np.log10(epsilon)).flatten()[use]
LOG_KAP = (np.log10(kappa)).flatten()[use]

# Plotting #
# Epsilon
d = getattr(Float, 'dist_ctd')[:, idxs].flatten(order='F')

tgps = getattr(Float, 'UTC_start')[idxs]
lon = getattr(Float, 'lon_start')[idxs]
lat = getattr(Float, 'lat_start')[idxs]
tctd = getattr(Float, 'UTC')[:, idxs].flatten(order='F')
nans = np.isnan(d) | np.isnan(tctd)
tctd = tctd[~nans]
dctd = d[~nans]
lonctd = np.interp(tctd, tgps, lon)
latctd = np.interp(tctd, tgps, lat)
bathy = sandwell.interp_track(lonctd, latctd, bf)


ax1.plot(Float.dist[idxs], 1000.*ieps)

sc = ax0.scatter(X, Z, s=5, c=LOG_EPS,
                 edgecolor='none', cmap=plt.get_cmap('bwr'), vmin=-11.,
                 vmax=-7)

ax1.set_ylabel('$P$ (mW m$^{-2}$)')
ax1.yaxis.set_ticks(np.array([0., 5., 10., 15]))
ax1.xaxis.set_ticks([])

ax0.fill_between(dctd[::100], bathy[::100],
                 np.nanmin(bathy), color='black', linewidth=2)
ax0.set_ylim(np.nanmin(bathy), 0.)
ax0.yaxis.set_ticks(np.arange(-4000, 1000, 1000))

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.7])
C = fig.colorbar(sc, cax=cbar_ax, extend='both')
C.set_label(r'$\log_{10}(\epsilon)$ (W kg$^{-1}$)')

#    plt.clim(-11., -7.)
ax0.set_xlim(np.min(X), np.max(X))
ax0.set_ylim(-5000., 0.)
ax0.set_xlabel('Distance from ridge top (km)')
ax0.set_ylabel('$z$ (m)')

ax1.set_xlim(*ax0.get_xlim())

#pf.my_savefig(fig, '4977', 'epsilon_lem_full', sdir, ftype='png', fsize='double_col')

# %% High and low energy spectra
isort = ieps.argsort()
imin = isort[1]
imax = isort[-6]

iZmax = iZ[:, imax]
iZmin = iZ[:, imin]

if xvar == 'time':
    x = np.arange(0., 10000., dx)
    __, __, w = Float.get_interp_grid(hpids, x, 'dUTC', 'Ww')
elif xvar == 'height':
    x = np.arange(-1500., 0., dx)
    __, __, w = Float.get_interp_grid(hpids, x, 'z', 'Ww')

wmax = w[:, imax]
wmin = w[:, imin]

nperseg = 600
noverlap = nperseg/2
m, Pmax = sp.signal.welch(wmax, dx, nperseg=nperseg, noverlap=noverlap)
m, Pmin = sp.signal.welch(wmin, dx, nperseg=nperseg, noverlap=noverlap)

use = m < 1./15
m, Pmax, Pmin = m[use], Pmax[use], Pmin[use]

fig, axs = plt.subplots(2, 1)
axs[0].loglog(1./m, Pmax)
axs[0].loglog(1./m, Pmin)
axs[1].semilogx(1./m, Pmax/Pmin)