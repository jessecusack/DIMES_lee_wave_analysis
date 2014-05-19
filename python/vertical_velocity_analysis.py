# -*- coding: utf-8 -*-
"""
Created on Thu May 15 12:07:00 2014

@author: jc3e13
"""

import emapex
import vertical_velocity_model as vvm
import numpy as np
import plotting_functions as pf
import matplotlib.pyplot as plt

reload(vvm)
reload(pf)


def assess_vvm_fit(Float):
    """ """

    vfi = Float.__vfi
    hpids = vfi.hpids
    floatID = Float.floatID

    # Histogram of vertical water velocity.

    Ww = Float.rWw.flatten(order='F')

    Ww_mean = np.nanmean(Ww)
    Ww_std = np.nanstd(Ww)

    plt.figure()
    bins = np.arange(-0.15, 0.16, 0.005)
    Ww_hist, bins, patches = plt.hist(Ww, bins=bins, histtype='stepfilled')
    plt.setp(patches, 'facecolor', 'b', 'alpha', 0.75)
    plt.xlim(np.min(bins), np.max(bins))
    plt.xlabel('$W_w$ (m s$^{-1}$)')
    title_str = ("Float {}, mean = {:1.2e} m s$^{{-1}}$, std = {:1.2e} "
                 "m s$^{{-1}}$").format(floatID, Ww_mean, Ww_std)
    plt.title(title_str)

    # Time series of different velocity measures.

    hpid1 = hpids[0]

    plt.figure()
    N = 4
    time, Ww = Float.get_timeseries(np.arange(hpid1, hpid1+N), 'rWw')
    __, Wz = Float.get_timeseries(np.arange(hpid1, hpid1+N), 'rWz')
    __, Ws = Float.get_timeseries(np.arange(hpid1, hpid1+N), 'rWs')
    plt.plot(time, Ww)
    plt.plot(time, Wz)
    plt.plot(time, Ws)
    plt.ylabel('$W_w$, $W_f$, $W_s$ (m s$^{-1}$)')
    plt.xlabel('Time')
    title_str = ("Float {}, half profiles {}").format(floatID, hpids[0:N])
    plt.title(title_str)
    plt.legend(['$W_w$', '$W_f$', '$W_s$'])

    # Distance section of water velocity.

#    dist_section(Float, np.arange(1, 1000), 'rWw')
#    plt.xlabel('Distance (km)')
#    plt.ylabel('Depth (m)')
#    plt.xlim(np.min(Float.dist), np.max(Float.dist))
#    title_str = ("Float {}").format(floatID)
#    plt.title(title_str)
#    cbar = plt.colorbar(orientation='horizontal')
#    cbar.set_label('$W_w$ (m s$^{-1}$)')


try:
    print("Floats {} and {}.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4976)
    E77 = emapex.EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4977)


model = '1'
cf_key = 'diffsq'
params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
fixed = [None, None, None, None, 1.156e-6, None, 27.179]
vvm.fitter(E76, params0, fixed, model=model, profiles='all', cf_key=cf_key)
assess_vvm_fit(E76)

model = '1'
cf_key = 'diffsq'
params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
fixed = [None, None, None, None, 1.156e-6, None, 27.179]
vvm.fitter(E77, params0, fixed, model=model, profiles='all', cf_key=cf_key)
assess_vvm_fit(E77)





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
