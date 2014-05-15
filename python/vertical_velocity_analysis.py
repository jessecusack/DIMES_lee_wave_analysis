# -*- coding: utf-8 -*-
"""
Created on Thu May 15 12:07:00 2014

@author: jc3e13
"""

import vertical_velocity_model as vvm
import numpy as np
import matplotlib.pyplot as plt

reload(vvm)


cf_key = 'diffsq'
model = '1'


params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16.])
vvm.fitter(E76, params0, model=model, profiles='down', cf_key=cf_key)
plt.title('Float 4976, model ' + model + ', down profiles only.')

params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16.])
vvm.fitter(E76, params0, model=model, profiles='up', cf_key=cf_key)
plt.title('Float 4976, model ' + model + ', up profiles only.')

params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16.])
vvm.fitter(E76, params0, model=model, profiles='all', cf_key=cf_key)
plt.title('Float 4976, model ' + model + ', all profiles.')


params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16.])
vvm.fitter(E77, params0, model=model, profiles='down', cf_key=cf_key)
plt.title('Float 4977, model ' + model + ', down profiles only.')

params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16.])
vvm.fitter(E77, params0, model=model, profiles='up', cf_key=cf_key)
plt.title('Float 4977, model ' + model + ', up profiles only.')

params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16.])
vvm.fitter(E77, params0, model=model, profiles='all', cf_key=cf_key)
plt.title('Float 4977, model ' + model + ', all profiles.')