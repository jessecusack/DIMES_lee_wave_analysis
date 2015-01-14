# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 12:58:56 2014

@author: jc3e13
"""

import numpy as np
import sys
import matplotlib.pyplot as plt

lib_path = '/noc/users/jc3e13/emapex/python'
pymc3_path = '/noc/users/jc3e13/envs/my_root/lib/python2.7/site-packages/pymc-3.0-py2.7.egg'
if lib_path not in sys.path:
    sys.path.append(lib_path)

# We want pymc 2.3
if pymc3_path in sys.path:
    sys.path.remove(pymc3_path)

import pymc
import emapex
import float_advection_routines as far
import utils

try:
    print("Floats {} and {} exist!".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.EMApexFloat('/noc/users/jc3e13/data/EM-APEX/allprofs11.mat', 4976)
    E76.apply_w_model('/noc/users/jc3e13/data/EM-APEX/4976_fix_p0k0M_fit_info.p')
    E76.apply_strain('/noc/users/jc3e13/data/EM-APEX/4976_N2_ref_300dbar.p')
    E76.apply_isopycnal_displacement('/noc/users/jc3e13/data/EM-APEX/srho_4976_100mbin.p')
    E77 = emapex.EMApexFloat('/noc/users/jc3e13/data/EM-APEX/allprofs11.mat', 4977)
    E77.apply_w_model('/noc/users/jc3e13/data/EM-APEX/4977_fix_p0k0M_fit_info.p')
    E77.apply_strain('/noc/users/jc3e13/data/EM-APEX/4977_N2_ref_300dbar.p')
    E77.apply_isopycnal_displacement('/noc/users/jc3e13/data/EM-APEX/srho_4977_100mbin.p')

pfl26 = E77.get_profiles(26)

zmax = -650.
use = ~np.isnan(pfl26.z) & (pfl26.z < zmax)
zf = pfl26.z[use]
wf = pfl26.Ww[use]
uf = pfl26.U_abs[use]
vf = pfl26.V_abs[use]
bf = pfl26.b[use]
zmin = np.min(zf)

uf = utils.nan_detrend(zf, uf, 0)
vf = utils.nan_detrend(zf, vf, 0)


data_stack = np.hstack((uf, vf, wf, 250*bf))


def model():

    # Priors.
    sig = pymc.Uniform('sig', 0.0, 5., value=0.01)
    X = pymc.Uniform('X', -2e4, 2e4, value=-2e3)
    Y = pymc.Uniform('Y', -1e5, 1e5, value=-4e3)
    Z = pymc.Uniform('Z', -2e4, 2e4, value=-2e3)

    @pymc.deterministic()
    def wave_model(zf=zf, X=X, Y=Y, Z=Z):
        return far.model_pymc(zf, X, Y, Z)

    # Likelihood
    y = pymc.Normal('y', mu=wave_model, tau=1./sig**2, value=data_stack, observed=True)

    return locals()

M = pymc.MCMC(model(), db='pickle', dbname='trace.p')
samples = 100000
burn = 40000
thin = 6
M.sample(samples, burn, thin)
pymc.Matplot.plot(M, common_scale=False)

plt.figure()
plt.plot(data_stack)
plt.plot(far.model_pymc(zf, np.median(M.trace('X')[:]), np.median(M.trace('Y')[:]),
                      np.median(M.trace('Z')[:])))
plt.savefig('output.png')

plt.figure()
plt.plot(data_stack)
plt.plot(far.model_pymc(zf, 5000., 15000, 8000, 0.))
