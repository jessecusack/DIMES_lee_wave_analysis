# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 10:30:46 2014

@author: jc3e13
"""

import finescale as fs
import GM79

reload(fs)
reload(GM79)

params = fs.default_params
periodogram_params = fs.default_periodogram_params

params['plot_results'] = False
params['plot_profiles'] = False
params['plot_spectra'] = False
params['periodogram_params']['nfft'] = None
params['periodogram_params']['window'] = 'hanning'
params['m_0'] = 1./120.
params['m_c'] = 1./12.
params['bin_width'] = 300.
params['bin_overlap'] = 250.

results = fs.analyse_float(E76, np.arange(1,60), params)
