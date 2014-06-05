# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 16:48:46 2014

@author: jc3e13
"""

import scipy.signal as sig
import emapex
import pickle


try:
    print("Floats {} and {}.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4976)
    E77 = emapex.EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4977)
    
with open('../../data/EM-APEX/4976_N2_ref_100dbar.p', 'rb') as f:
    N2_ref = pickle.load(f)
    setattr(E76, 'N2_ref', N2_ref)
with open('../../data/EM-APEX/4977_N2_ref_100dbar.p', 'rb') as f:
    N2_ref = pickle.load(f)
    setattr(E77, 'N2_ref', N2_ref)
    
    