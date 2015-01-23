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

for Float in [E76, E77]:

    epsilon, kappa = fs.w_scales_float(Float)

    fig = plt.figure()
    plt.scatter(Float.r_dist_ctd.flatten(), Float.r_z.flatten(),
                c=(np.log10(epsilon)).flatten(), edgecolor='none')
    plt.colorbar()
    plt.clim(-11, -5)
    plt.xlim(np.min(Float.r_dist_ctd), np.max(Float.r_dist_ctd))
    plt.ylim(-1400, -100)
    plt.xlabel('Distance (km)')
    plt.ylabel('$z$ (m)')
    plt.title('log10 turbulent kinetic energy dissipation (W/kg)')
    pf.my_savefig(fig, Float.floatID, 'epsilon', sdir)

    fig = plt.figure()
    plt.scatter(Float.r_dist_ctd.flatten(), Float.r_z.flatten(),
                c=(np.log10(kappa)).flatten(), edgecolor='none')
    plt.colorbar()
    plt.clim(-6, -3)
    plt.xlim(np.min(Float.r_dist_ctd), np.max(Float.r_dist_ctd))
    plt.ylim(-1400, -100)
    plt.xlabel('Distance (km)')
    plt.ylabel('$z$ (m)')
    plt.title('log10 diapycnal diffusivity (m2/s)')
    pf.my_savefig(fig, Float.floatID, 'kappa', sdir)
