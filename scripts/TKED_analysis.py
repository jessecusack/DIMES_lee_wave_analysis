# -*- coding: utf-8 -*-
"""
Created on Fri Jan 23 16:25:26 2015

@author: jc3e13
"""

import numpy as np
import sys
import os
import glob
import matplotlib
import matplotlib.pyplot as plt

lib_path = '/noc/users/jc3e13/emapex/python'
if lib_path not in sys.path:
    sys.path.append(lib_path)

import emapex
import finescale as fs
import plotting_functions as pf
import sandwell


try:
    print("Floats {} and {} exist!".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.EMApexFloat('/noc/users/jc3e13/data/EM-APEX/allprofs11.mat', 4976)
    E76.apply_w_model('/noc/users/jc3e13/data/EM-APEX/4976_fix_p0k0M_fit_info.p')
    E76.apply_strain('../../data/EM-APEX/4976_N2_ref_300dbar.p')
    E76.apply_isopycnal_displacement('../../data/EM-APEX/srho_4976_100mbin.p')
    E76.generate_regular_grids(zmin=-1500., dz=1.)
    E77 = emapex.EMApexFloat('/noc/users/jc3e13/data/EM-APEX/allprofs11.mat', 4977)
    E77.apply_w_model('/noc/users/jc3e13/data/EM-APEX/4977_fix_p0k0M_fit_info.p')
    E77.apply_strain('../../data/EM-APEX/4977_N2_ref_300dbar.p')
    E77.apply_isopycnal_displacement('../../data/EM-APEX/srho_4977_100mbin.p')
    E77.generate_regular_grids(zmin=-1500., dz=1.)


# %% Script params.

# Bathymetry file path.
bf = os.path.abspath(glob.glob('../../data/sandwell_bathymetry/topo_*.img')[0])
# Figure save path.
sdir = '../figures/TKED_estimation'
if not os.path.exists(sdir):
    os.makedirs(sdir)
# Universal figure font size.
matplotlib.rc('font', **{'size': 9})


# %% Start script
# Using large eddy method first.

for Float in [E76, E77]:

    bathy = sandwell.interp_track(Float.lon_start, Float.lat_start, bf)

    epsilon, kappa = fs.w_scales_float(Float)

    fig = plt.figure(figsize=(10, 6))

    Z = Float.r_z.flatten()

    use = (Z < -100) & (Z > -1400)

    Z = Z[use]
    X = Float.r_dist_ctd.flatten()[use]
    LOG_EPS = (np.log10(epsilon)).flatten()[use]
    LOG_KAP = (np.log10(kappa)).flatten()[use]

    plt.scatter(X, Z, s=5, c=LOG_EPS, edgecolor='none', cmap=plt.get_cmap('bwr'))
    C = plt.colorbar(extend='both')
    C.set_label(r'$\log_{10}(\epsilon)$ (W kg$^{-1}$)')
    plt.clim(-11, -7)
    plt.xlim(np.min(X), np.max(X))
    plt.ylim(np.nanmin(bathy), 0.)
    plt.xlabel('Distance (km)')
    plt.ylabel('$z$ (m)')
    plt.fill_between(Float.dist, bathy, np.nanmin(bathy), color='black',
                     linewidth=2)
    pf.my_savefig(fig, Float.floatID, 'epsilon', sdir)

    ##

    fig = plt.figure(figsize=(10, 6))
    plt.scatter(X, Z, s=5, c=LOG_KAP, edgecolor='none', cmap=plt.get_cmap('bwr'))
#    plt.pcolormesh(Float.r_dist_ctd, Float.r_z, np.log10(kappa))
    C = plt.colorbar(extend='both')
    C.set_label(r'$\log_{10}(\kappa_\rho)$ (m$^2$ s$^{-1}$)')
    plt.clim(-5, -3)
    plt.xlim(np.min(X), np.max(X))
    plt.ylim(np.nanmin(bathy), 0.)
    plt.xlabel('Distance (km)')
    plt.ylabel('$z$ (m)')
    plt.fill_between(Float.dist, bathy, np.nanmin(bathy), color='black',
                     linewidth=2)
    pf.my_savefig(fig, Float.floatID, 'kappa', sdir)

# %% Using finescale parameterisation

params = fs.default_params

params['plot_results'] = False
params['plot_profiles'] = False
params['plot_spectra'] = False
params['print_diagnostics'] = False
params['periodogram_params']['nfft'] = None
params['periodogram_params']['window'] = 'hanning'
params['m_0'] = 1./120.
params['m_c'] = 1./12.
params['bin_width'] = 200.
params['bin_overlap'] = 100.
params['apply_corrections'] = True
params['zmin'] = -1400
params['zmax'] = -100

hpids = np.arange(10, 300)
results = fs.analyse_float(E76, hpids, params)
__, idxs = E76.get_profiles(hpids, ret_idxs=True)
dists = E76.dist[idxs]


ylims = (params['zmin'], params['zmax'])


pf.scatter_section(E76, hpids, 'Ww', cmap=plt.get_cmap('bwr'))
plt.clim(-0.1, 0.1)
plt.savefig('../figures/finescale/Ww.png', bbox_inches='tight')

plt.figure()
plt.title('log10 R_pol')
for result, dist in zip(results, dists):
    z_mean, EK, R_pol, R_om, epsilon, kappa = result
    d = dist*np.ones_like(z_mean)
    # PCOLORMESH DOESN'T WORK UNTIL GRID MADE...
    plt.pcolormesh(d, z_mean, np.log10(R_pol), cmap=plt.get_cmap('bwr'))

plt.clim(-1, 1)
plt.colorbar()
plt.ylim(-1400., -100.)
plt.xlim(np.min(dists), np.max(dists))
plt.savefig('../figures/finescale/R_pol.png', bbox_inches='tight')

plt.figure()
plt.title('log10 epsilon')
for result, dist in zip(results, dists):
    z_mean, EK, R_pol, R_om, epsilon, kappa = result
    d = dist*np.ones_like(z_mean)
    plt.scatter(d, z_mean, c=np.log10(epsilon), edgecolor='none',
                cmap=plt.get_cmap('jet'))

plt.colorbar()
plt.ylim(-1400., -100.)
plt.xlim(np.min(dists), np.max(dists))
plt.savefig('../figures/finescale/epsilon.png', bbox_inches='tight')

plt.figure()
plt.title('log10 kappa')
for result, dist in zip(results, dists):
    z_mean, EK, R_pol, R_om, epsilon, kappa = result
    d = dist*np.ones_like(z_mean)
    plt.scatter(d, z_mean, c=np.log10(kappa), edgecolor='none',
                cmap=plt.get_cmap('jet'))

plt.colorbar()
plt.ylim(-1400., -100.)
plt.xlim(np.min(dists), np.max(dists))
plt.savefig('../figures/finescale/kappa.png', bbox_inches='tight')

plt.figure()
plt.title('R_om')
for result, dist in zip(results, dists):
    z_mean, EK, R_pol, R_om, epsilon, kappa = result
    d = dist*np.ones_like(z_mean)
    plt.scatter(d, z_mean, c=R_om, edgecolor='none',
                cmap=plt.get_cmap('jet'))

plt.colorbar()
plt.ylim(-1400., -100.)
plt.xlim(np.min(dists), np.max(dists))
plt.savefig('../figures/finescale/R_om.png', bbox_inches='tight')

plt.figure()
plt.title('log10 EK')
for result, dist in zip(results, dists):
    z_mean, EK, R_pol, R_om, epsilon, kappa = result
    d = dist*np.ones_like(z_mean)
    plt.scatter(d, z_mean, c=np.log10(EK), edgecolor='none',
                cmap=plt.get_cmap('jet'))

plt.colorbar()
plt.ylim(-1400., -100.)
plt.xlim(np.min(dists), np.max(dists))
plt.savefig('../figures/finescale/EK.png', bbox_inches='tight')