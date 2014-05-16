# -*- coding: utf-8 -*-
"""
Created on Thu May 15 12:07:00 2014

@author: jc3e13
"""

import vertical_velocity_model as vvm
import numpy as np
import plotting_functions as pf

reload(vvm)
reload(pf)

model = '1'
cf_key = 'diffsq'
params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
fixed = [None, None, None, None, 1.156e-6, None, 27.179]
vvm.fitter(E76, params0, fixed, model=model, profiles='all', cf_key=cf_key)
pf.assess_vvm_fit(E76)

model = '1'
cf_key = 'diffsq'
params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
fixed = [None, None, None, None, 1.156e-6, None, 27.179]
vvm.fitter(E77, params0, fixed, model=model, profiles='all', cf_key=cf_key)
pf.assess_vvm_fit(E77)





#params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16.])
#vvm.fitter(E77, params0, model=model, profiles='down', cf_key=cf_key)
#plt.title('Float 4977, model ' + model + ', down profiles only.')
#
#params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16.])
#vvm.fitter(E77, params0, model=model, profiles='up', cf_key=cf_key)
#plt.title('Float 4977, model ' + model + ', up profiles only.')


#params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16.])
#vvm.fitter(E76, params0, model=model, profiles='down', cf_key=cf_key)
#plt.title('Float 4976, model ' + model + ', down profiles only.')
#
#
#params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16.])
#vvm.fitter(E76, params0, model=model, profiles='up', cf_key=cf_key)
#plt.title('Float 4976, model ' + model + ', up profiles only.')