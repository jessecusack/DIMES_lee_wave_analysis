# -*- coding: utf-8 -*-
"""
Created on Thu May 15 12:07:00 2014

@author: jc3e13
"""

import numpy as np
import matplotlib
import os

import emapex
import vertical_velocity_fitter as vvf


# Figure save path.
sdir = os.path.join('..', 'figures', 'vertical_velocity_analysis')
if not os.path.exists(sdir):
    os.makedirs(sdir)
# Universal figure font size.
matplotlib.rc('font', **{'size': 9})


###############################################################################

try:
    print("Floats {} and {} exist!".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.load(4976, apply_w=False)
    E77 = emapex.load(4977, apply_w=False)

# %% ##########################################################################

# cf_key = 'diffsq'
# params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
# fixed = [None, None, None, None, None, None, None]
# wfi = vvf.fitter(E76, params0, fixed, model=model, profiles='all',
#                 cf_key=cf_key)
# E76.apply_w_model(wfi)
# vvf.assess_w_fit(E76, str(E76.floatID))
# print(E76.__wfi.p)
#
# cf_key = 'diffsq'
# params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
# fixed = [None, None, None, None, None, None, None]
# wfi = vvf.fitter(E77, params0, fixed, model=model, profiles='all',
#                 cf_key=cf_key)
# E77.apply_w_model(wfi)
# vvf.assess_w_fit(E77, str(E77.floatID))
# print(E77.__wfi.p)

# %% ##########################################################################

cf_key = 'diffsq'
params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
fixed = [None, None, None, None, None, None, 27.179]
wfi = vvf.fitter(E76, params0, fixed, model=model, profiles='all',
                 cf_key=cf_key)
E76.apply_w_model(wfi)
vvf.assess_w_fit(E76, str(E76.floatID)+'_fix_M')
print(E76.__wfi['p'])

cf_key = 'diffsq'
params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
fixed = [None, None, None, None, None, None, 27.179]
wfi = vvf.fitter(E77, params0, fixed, model=model, profiles='all',
                 cf_key=cf_key)
E77.apply_w_model(wfi)
vvf.assess_w_fit(E77, str(E77.floatID)+'_fix_M')
print(E77.__wfi['p'])

# %% ##########################################################################

cf_key = 'diffsq'
params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
fixed = [None, None, None, 2000., None, 16, 27.179]
Plims = (60., 1500.)
hpids = np.arange(50, 151)

wfi = vvf.fitter(E76, hpids, params0, fixed, Plims=Plims, profiles='all',
                 cf_key=cf_key)
E76.apply_w_model(wfi)
vvf.assess_w_fit(E76, str(E76.floatID)+'_fix_p0k0M')
print(E76.__wfi['p'])

wfi = vvf.fitter(E77, hpids, params0, fixed, Plims=Plims, profiles='all',
                 cf_key=cf_key)
E77.apply_w_model(wfi)
vvf.assess_w_fit(E77, str(E77.floatID)+'_fix_p0k0M')
print(E77.__wfi['p'])

# %% ##########################################################################

cf_key = 'diffsq'
params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
fixed = [None, None, None, None, 1.156e-6, None, 27.179]
wfi = vvf.fitter(E76, params0, fixed, model=model, profiles='all',
                 cf_key=cf_key)
E76.apply_w_model(wfi)
vvf.assess_w_fit(E76, str(E76.floatID)+'_fix_alphakM')
print(E76.__wfi['p'])

cf_key = 'diffsq'
params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
fixed = [None, None, None, None, 1.156e-6, None, 27.179]
wfi = vvf.fitter(E77, params0, fixed, model=model, profiles='all',
                 cf_key=cf_key)
E77.apply_w_model(wfi)
vvf.assess_w_fit(E77, str(E77.floatID)+'_fix_alphakM')
print(E77.__wfi['p'])

# %% ##########################################################################
#
# cf_key = 'diffsq'
# params0 = np.array([3e-2, 5e-2, 3e-6, 2e+3, 1e-6, 16., 27.179])
# fixed = [None, None, None, 2e+3, None, None, 27.179]
# wfi = vvf.fitter(E76, params0, fixed, model=model, profiles='all',
#                 cf_key=cf_key)
# E76.apply_w_model(wfi)
# vvf.assess_w_fit(E76, str(E76.floatID)+'_fix_p0M')
# print(E76.__wfi.p)
#
# cf_key = 'diffsq'
# params0 = np.array([3e-2, 5e-2, 2e-6, 2e+3, 1e-6, 16., 27.179])
# fixed = [None, None, None, 2e+3, None, None, 27.179]
# wfi = vvf.fitter(E77, params0, fixed, model=model, profiles='all',
#                 cf_key=cf_key)
# E77.apply_w_model(wfi)
# vvf.assess_w_fit(E77, str(E77.floatID)+'_fix_p0M')
# print(E77.__wfi.p)
#
# %% ##########################################################################
#
# cf_key = 'diffsq'
# params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
# fixed = [None, None, None, None, None, None, None]
# wfi = vvf.fitter(E76, params0, fixed, model=model, profiles='updown',
#                 cf_key=cf_key)
# E76.apply_w_model(wfi)
# vvf.assess_w_fit(E76, str(E76.floatID))
# print(E76.__wfi.p)
#
# cf_key = 'diffsq'
# params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
# fixed = [None, None, None, None, None, None, None]
# wfi = vvf.fitter(E77, params0, fixed, model=model, profiles='updown',
#                 cf_key=cf_key)
# E77.apply_w_model(wfi)
# vvf.assess_w_fit(E77, str(E77.floatID))
# print(E77.__wfi.p)
#
# %% ##########################################################################
#
# cf_key = 'diffsq'
# params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
# fixed = [None, None, None, None, 1.156e-6, None, 27.179]
# wfi = vvf.fitter(E76, params0, fixed, model=model, profiles='updown',
#                 cf_key=cf_key)
# E76.apply_w_model(wfi)
# vvf.assess_w_fit(E76, str(E76.floatID)+'_fix_alphakM')
# print(E76.__wfi.p)
#
# cf_key = 'diffsq'
# params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
# fixed = [None, None, None, None, 1.156e-6, None, 27.179]
# wfi = vvf.fitter(E77, params0, fixed, model=model, profiles='updown',
#                 cf_key=cf_key)
# E77.apply_w_model(wfi)
# vvf.assess_w_fit(E77, str(E77.floatID)+'_fix_alphakM')
# print(E77.__wfi.p)

