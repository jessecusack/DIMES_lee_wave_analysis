# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 16:48:46 2014

@author: jc3e13
"""

import scipy.signal as sig
import numpy as np
import emapex
import pickle


try:
    
    print("Floats {} and {}.".format(E76.floatID, E77.floatID))
    
except NameError:
    
    E76 = emapex.EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4976)
    E76.apply_w_model('../../data/EM-APEX/4976_fix_p0k0M_fit_info.p')
    
    E77 = emapex.EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4977)
    E77.apply_w_model('../../data/EM-APEX/4977_fix_p0k0M_fit_info.p')

    for Float, fstr in zip([E76, E77], ['76', '77']):
        with open('../../data/EM-APEX/49'+fstr+'_N2_ref_100dbar.p', 'rb') as f:
            N2_ref = pickle.load(f)
            setattr(Float, 'N2_ref', N2_ref)
            setattr(Float, 'strain_z', (Float.N2 - N2_ref)/N2_ref)
            Float.update_profiles()
    
# Do this for one profile initially.
P78 = E77.get_profiles(70)

zef = P78.zef
U_abs = P78.U_abs
V_abs = P78.V_abs

nans = np.isnan(zef) | np.isnan(U_abs) | np.isnan(V_abs)

zef = zef[~nans]
U_abs = U_abs[~nans]
V_abs = V_abs[~nans]

zmin, zmax = np.min(zef), np.max(zef)
dkmin = 1./np.round(zmax - zmin)

dz = 5.
dk = 1./dz
k = np.linspace(dkmin, dk, 470)

Pv = sig.lombscargle(zef, V_abs**2, k)
