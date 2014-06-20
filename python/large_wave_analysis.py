# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 16:48:46 2014

@author: jc3e13
"""

import scipy.signal as sig
import scipy.optimize as op
import matplotlib.pyplot as plt
import pylab as pyl
import numpy as np
import emapex
import pickle


try:

    print("Floats {} and {}.".format(E76.floatID, E77.floatID))

except NameError:
    
    reload(emapex)

    E76 = emapex.EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4976)
    E76.apply_w_model('../../data/EM-APEX/4976_fix_p0k0M_fit_info.p')

    E77 = emapex.EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4977)
    E77.apply_w_model('../../data/EM-APEX/4977_fix_p0k0M_fit_info.p')

    for Float, fstr in zip([E76, E77], ['76', '77']):
        with open('../../data/EM-APEX/49'+fstr+'_N2_ref_300dbar.p', 'rb') as f:
            N2_ref = pickle.load(f)
            setattr(Float, 'N2_ref', N2_ref)
            setattr(Float, 'strain_z', (Float.N2 - N2_ref)/N2_ref)
            Float.update_profiles()


def plane_wave(x, A, k, phase, C):
    return A*np.cos(2*np.pi*(k*x + phase)) + C

dz = 5.
dk = 1./dz

E76_hpids = np.arange(27, 34)
E77_hpids = np.arange(23, 30)

var_names = ['rWw']  #['rU_abs', 'rV_abs', 'rWw']

for Float, hpids in zip([E76, E77], [E76_hpids, E77_hpids]):

    for pfl in Float.get_profiles(hpids):

        z = getattr(pfl, 'rz')

        for name in var_names:

            var = getattr(pfl, name)
            N = var.size

            # Try fitting plane wave.
            popt, __ = op.curve_fit(plane_wave, z, var, 
                                    p0=[0.1, 1.6e-3, 0., 0.])
            mfit = popt[1]
            
            plt.figure()
            plt.subplot(1, 2, 1)            
            
            plt.plot(var, z, plane_wave(z, *popt), z)
            plt.xticks(rotation=45)
            plt.xlabel('$W_w$ (m s$^{-1}$)')
            plt.ylabel('$z$ (m)')
            title = ("{}: Float {}, profile {}\nm_fit {:1.1e} m$^{{-1}}$"
                     "\nlambda_fit {:4.0f} m"
                     "").format(name, Float.floatID, pfl.hpid[0], mfit, 1./mfit)
            plt.title(title)

            plt.subplot(1, 2, 2)
            # Try calculating power spectral density.
            Pzz, ms = plt.psd(var, NFFT=N//2, Fs=dk,
                              noverlap=int(0.2*N//2),
                              detrend=pyl.detrend_linear)
            plt.gca().set_xscale('log')
            mmax = ms[Pzz.argmax()]
            title = ("{}: Float {}, profile {}\nm_max {:1.1e} m$^{{-1}}$"
                     "\nlambda_max {:4.0f} m"
                     "").format(name, Float.floatID, pfl.hpid, mmax, 1./mmax)
            ylim = plt.ylim()
            plt.plot(2*[mmax], ylim, 'r', 2*[mfit], ylim, 'g')
            plt.title(title)
            plt.xlabel('$m$ (m$^{-1}$)')

def plane_wave2(params, x):
    A, k, m, om = params
    return A*np.cos(2*np.pi*(k*x[:,0] + m*x[:,1] + om*x[:,2]))
    
def cost(params, data, func, y):
    return (func(params, data) - y).flatten()
    
res = []

E76_hpids = np.arange(27, 34) # np.arange(31, 33)
E77_hpids = np.arange(23, 30) # np.arange(26, 28)
    
for Float, hpids in zip([E76, E77], [E76_hpids, E77_hpids]):
    
    __, idxs = Float.get_profiles(hpids, ret_idxs=True)
    
    z = Float.rz[:,idxs].flatten(order='F')
    x = Float.rdist_ctd_data[:, idxs].flatten(order='F')*1000.
    t = Float.rUTC[:, idxs].flatten(order='F')*86400.
    W = Float.rWw[:, idxs].flatten(order='F')
    
    nans = np.isnan(z) | np.isnan(x) | np.isnan(t) | np.isnan(W) | (z > -200)
    data = np.array([x[~nans], z[~nans], t[~nans]]).T
    W = W[~nans]
    
    x0 = [0.15, 0.001, 0.003, 0.]
    fit = op.leastsq(cost, x0=x0, args=(data, plane_wave2, W))[0]
    print(fit)
    res.append(fit)
    
    Wm = plane_wave2(fit, data)
    Wm0 = plane_wave2(x0, data)
    
    plt.figure()
    plt.subplot(3, 1, 1)
    plt.plot(data[:,0], Wm, data[:,0], W, data[:,0], Wm0)
    plt.subplot(3, 1, 2)
    plt.plot(Wm, data[:,1], W, data[:,1], Wm0, data[:,1])
    plt.subplot(3, 1, 3)
    plt.plot(data[:,2], Wm, data[:,2], W, data[:,2], Wm0)
    