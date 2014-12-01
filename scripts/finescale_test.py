# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 10:30:46 2014

@author: jc3e13
"""

import os
import sys

lib_path = os.path.abspath('../python')
if lib_path not in sys.path:
    sys.path.append(lib_path)

import finescale as fs
import GM79
import emapex

reload(fs)
reload(GM79)
reload(emapex)

try:
    print("Float {} exist!".format(E76.floatID))
except NameError:
    E76 = emapex.EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4976)
    E76.apply_w_model('../../data/EM-APEX/4976_fix_p0k0M_fit_info.p')
    E76.apply_strain('../../data/EM-APEX/4976_N2_ref_300dbar.p')
    E76.apply_isopycnal_displacement('../../data/EM-APEX/srho_4976_100mbin.p')

# %%

params = fs.default_params

params['plot_results'] = True
params['plot_profiles'] = True
params['plot_spectra'] = False
params['print_diagnostics'] = True
params['periodogram_params']['nfft'] = None
params['periodogram_params']['window'] = 'hanning'
params['m_0'] = 1./120.
params['m_c'] = 1./12.
params['bin_width'] = 200.
params['bin_overlap'] = 100.
params['apply_corrections'] = True
params['zmin'] = -1400
params['zmax'] = -200

results = fs.analyse_profile(E76.get_profiles(99), params)
