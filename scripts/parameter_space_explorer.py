# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15 14:43:24 2014

@author: jc3e13
"""

import os
import sys
import glob
import itertools
import pickle

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
# from matplotlib.colors import LogNorm
# from matplotlib.ticker import LogFormatterMathtext

lib_path = os.path.abspath('../modules')
if lib_path not in sys.path:
    sys.path.append(lib_path)

import emapex
import plotting_functions as pf
import float_advection_routines as far


try:
    print("Floats {} and {} exist!.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.load(4976)
    E77 = emapex.load(4977)

# %% Script params.

# Bathymetry file path.
bf = os.path.abspath(glob.glob('/noc/users/jc3e13/storage/smith_sandwell/topo_*.img')[0])
# Figure save path.
sdir = '../figures/parameter_space_explorer'
if not os.path.exists(sdir):
    os.makedirs(sdir)
# Universal figure font size.
matplotlib.rc('font', **{'size': 9})

# %% ##########################################################################

# Parameters to check.
Xs = np.arange(-10000., 11000., 1000.)
Ys = np.arange(-10000., 11000., 1000.)
Zs = np.arange(-6000., 6500., 500.)
Phases = np.linspace(0., 2.*np.pi, 10)

params = far.default_params
params['z_0'] = -1600.
pfl = E76.get_profiles(32)
use = ~np.isnan(pfl.zef) & (pfl.zef < -400.)
bscale = 250.

zf = pfl.zef[use]
uf = pfl.interp(zf, 'zef', 'U_abs')
vf = pfl.interp(zf, 'zef', 'V_abs')
wf = pfl.interp(zf, 'zef', 'Ww')
bf = bscale*pfl.interp(zf, 'z', 'b')

ps = []
cost = []

for LX, LY, LZ, Phase in itertools.product(Xs, Ys, Zs, Phases):

    ps.append((LX, LY, LZ, Phase))

    if LX == 0. or LY == 0. or LZ == 0.:
        cost.append(1e10)
        continue

    X = far.model_verbose(LX, LY, LZ, Phase, params)

    um = np.interp(zf, X.r[:, 2], X.u[:, 0])
    vm = np.interp(zf, X.r[:, 2], X.u[:, 1])
    wm = np.interp(zf, X.r[:, 2], X.u[:, 2])
    bm = bscale*np.interp(zf, X.r[:, 2], X.b)

    c = np.std(um - uf) + np.std(vm - vf) + np.std(wm - wf) + np.std(bm - bf)
    cost.append(c)

    print(LX, LY, LZ, Phase, c)

with open('pfl32_param_search_cost.p', 'wb') as f:
    pickle.dump(c, f)

with open('pfl_32param_search_params.p', 'wb') as f:
    pickle.dump(ps, f)
