# -*- coding: utf-8 -*-
"""
This script loads floats and imports useful modules into the namespace.

Created on Thu May 29 14:52:32 2014

@author: jc3e13
"""

import numpy as np
import matplotlib.pyplot as plt
import plotting_functions as pf
import emapex
import gsw
import os
import sandwell
from scipy.interpolate import griddata
from scipy.integrate import cumtrapz
import scipy.signal as sig

reload(emapex)

E76 = emapex.EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4976)
E76.apply_w_model('../../data/EM-APEX/4976_fix_p0k0M_fit_info.p')
E76.apply_strain('../../data/EM-APEX/4976_N2_ref_300dbar.p')
E76.apply_isopycnal_displacement('../../data/EM-APEX/srho_4976_100mbin.p')

E77 = emapex.EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4977)
E77.apply_w_model('../../data/EM-APEX/4977_fix_p0k0M_fit_info.p')
E77.apply_strain('../../data/EM-APEX/4977_N2_ref_300dbar.p')
E77.apply_isopycnal_displacement('../../data/EM-APEX/srho_4977_100mbin.p')
