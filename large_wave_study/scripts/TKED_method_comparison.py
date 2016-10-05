# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 15:20:57 2016

@author: jc3e13
"""

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import emapex
import TKED_parameterisations as TKED
import plotting_functions as pf

zmin = -1450.
dz = 1.

try:
    print("Floats {} and {} exist!.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.load(4976)
#    E76.generate_regular_grids(zmin=zmin, dz=dz)
    E77 = emapex.load(4977)
#    E77.generate_regular_grids(zmin=zmin, dz=dz)

# %% Script params.

# Figure save path.
sdir = '../figures/TKED_estimation'
if not os.path.exists(sdir):
    os.makedirs(sdir)
# Universal figure font size.
matplotlib.rc('font', **{'size': 8})

# %%

pfl = E77.get_profiles(352)
nnan = ~np.isnan(pfl.P)
rho_1 = pfl.rho_1[nnan]
z = pfl.z[nnan]
N2 = pfl.N2[nnan]
w = pfl.Ww[nnan]
wz = pfl.Wz[nnan]
ws = pfl.Ws[nnan]

rho_1s = np.sort(rho_1)[::-1]

thorpe_scales, thorpe_disp, x_sorted, idxs = TKED.thorpe_scales(z, rho_1)

fig, axs = plt.subplots(1, 5, sharey='row')
axs[0].plot(rho_1, z)
axs[0].plot(rho_1s, z)
axs[1].vlines(0., *axs[2].get_ylim())
axs[1].plot(N2, z)
axs[2].vlines(0., *axs[2].get_ylim())
axs[2].plot(w, z)
axs[3].plot(wz, z)
axs[3].plot(ws, z)
axs[4].plot(thorpe_scales, z)