#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 14:41:56 2017

@author: jc3e13
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import scipy.optimize as opt
import gsw

import emapex


# Figure save path.
sdir = os.path.join('..', 'figures', 'w_analysis_complex')
if not os.path.exists(sdir):
    os.makedirs(sdir)
# Universal figure font size.
matplotlib.rc('font', **{'size': 7})


###############################################################################

try:
    print("Floats {} and {} exist!".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.load(4976, apply_w=False)
    E77 = emapex.load(4977, apply_w=False)


# %%
def float_volume(p, k, V0=2.62e-2, alpha_p=3.6e-6, p0=2000., alpha_k=1.2e-6,
                 k0=16.):
    return V0*(1 - alpha_p*(p - p0) + alpha_k*(k - k0))


def w_cost(w, wf, p, k, V0, alpha_p, p0, alpha_k, k0, M, dwfdt, g, rho, C_D, A):
    V = float_volume(p, k, V0, alpha_p, p0, alpha_k, k0)
    accel = M*dwfdt
    buoy = g*(M - rho*V)
    drag = rho*C_D*A*(wf - w)*np.abs(wf - w)
    return (accel - buoy + drag)**2


def calc_w(wf, p, k, V0, alpha_p, p0, alpha_k, k0, M, dwfdt, g, rho, C_D, A):
    wf_ = wf.flatten()
    p_ = p.flatten()
    k_ = k.flatten()
    dwfdt_ = dwfdt.flatten()
    g_ = g.flatten()
    rho_ = rho.flatten()
    w_ = np.full_like(p_, np.nan)
    ni = len(w_)

    for i in range(ni):
        if np.isnan(dwfdt_[i]):
            continue

        args = (wf_[i], p_[i], k_[i], V0, alpha_p, p0, alpha_k, k0, M,
                dwfdt_[i], g_[i], rho_[i], C_D, A)
        res = opt.minimize(w_cost, 0., args=args, method='Nelder-Mead')
        w_[i] = res.x

    return w_.reshape(p.shape)


def calc_w_2(wf, p, k, V0, alpha_p, p0, alpha_k, k0, M, dwfdt, g, rho, C_D, A):
    wf_ = wf.flatten()
    p_ = p.flatten()
    k_ = k.flatten()
    dwfdt_ = dwfdt.flatten()
    g_ = g.flatten()
    rho_ = rho.flatten()
    w_ = np.full_like(p_, np.nan)

    nnans = ~np.isnan(dwfdt_)

    wf_ = wf_[nnans]
    p_ = p_[nnans]
    k_ = k_[nnans]
    dwfdt_ = dwfdt_[nnans]
    g_ = g_[nnans]
    rho_ = rho_[nnans]

    args = (wf_, p_, k_, V0, alpha_p, p0, alpha_k, k0, M, dwfdt_, g_, rho_,
            C_D, A)
    x0 = np.zeros_like(wf_)
    res = opt.fsolve(w_cost, x0, args=args)
    w_[nnans] = res.x

    return w_.reshape(p.shape)


def calc_w_ud(wf, p, k, V0, alpha_p, p0, alpha_k, k0, M, dwfdt, g, rho, C_D_u, C_D_d, A, up):
    wf_ = wf.flatten()
    p_ = p.flatten()
    k_ = k.flatten()
    dwfdt_ = dwfdt.flatten()
    g_ = g.flatten()
    rho_ = rho.flatten()
    w_ = np.full_like(p_, np.nan)
    C_D = (up*C_D_u + ~up*C_D_d).flatten()
    ni = len(w_)

    for i in range(ni):
        if np.isnan(dwfdt_[i]):
            continue

        args = (wf_[i], p_[i], k_[i], V0, alpha_p, p0, alpha_k, k0, M,
                dwfdt_[i], g_[i], rho_[i], C_D[i], A)
        res = opt.minimize(w_cost, 0., args=args, method='Nelder-Mead')
        w_[i] = res.x

    return w_.reshape(p.shape)


def cost(params, wf, p, k, V0, p0, k0, M, dwfdt, g, rho, A):
    C_D, alpha_p, alpha_k = params
    w = calc_w(wf, p, k, V0, alpha_p, p0, alpha_k, k0, M, dwfdt, g, rho, C_D, A)
    return (w**2).sum()


def cost_ud(params, wf, p, k, V0, p0, k0, M, dwfdt, g, rho, A, up):
    C_D_u, C_D_d, alpha_p, alpha_k = params
    w = calc_w_ud(wf, p, k, V0, alpha_p, p0, alpha_k, k0, M, dwfdt, g, rho, C_D_u, C_D_d, A, up)
    return (w**2).sum()


def fit(Float, hpids, pmin, pmax, V0, p0, k0, M, A):

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)

    lat = Float.lat_start[idxs]
    p = Float.P[:, idxs]
    g = -gsw.grav(lat, p)
    rho = Float.rho[:, idxs]
    k = Float.ppos[:, idxs]
    wf = Float.Wz[:, idxs]
    t = Float.UTC[:, idxs]*86400.
    dwfdt = np.gradient(wf, axis=0)/np.gradient(t, axis=0)

    use = (p > pmin) & (p < pmax) & ~np.isnan(dwfdt)

    args = (wf[use], p[use], k[use], V0, p0, k0, M, dwfdt[use], g[use], rho[use], A)
    res = opt.minimize(cost, [1.5, 3.6e-6, 1.1e-6], args=args)

    return res


def fit_ud(Float, hpids, pmin, pmax, V0, p0, k0, M, A):

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)

    lat = Float.lat_start[idxs]
    p = Float.P[:, idxs]
    g = -gsw.grav(lat, p)
    rho = Float.rho[:, idxs]
    k = Float.ppos[:, idxs]
    wf = Float.Wz[:, idxs]
    t = Float.UTC[:, idxs]*86400.
    dwfdt = np.gradient(wf, axis=0)/np.gradient(t, axis=0)
    hpids = Float.hpid[idxs]
    up = np.full_like(p, True, dtype=bool)
    up[:, hpids % 2 != 0] = False

    use = (p > pmin) & (p < pmax) & ~np.isnan(dwfdt)

    args = (wf[use], p[use], k[use], V0, p0, k0, M, dwfdt[use], g[use], rho[use], A, up[use])
    res = opt.minimize(cost_ud, [1.5, 1.5, 3.6e-6, 1.1e-6], args=args)

    return res


def apply_w_fit(Float, hpids, res, V0, p0, k0, M, A):

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)

    lat = Float.lat_start[idxs]
    p = Float.P[:, idxs]
    g = -gsw.grav(lat, p)
    rho = Float.rho[:, idxs]
    k = Float.ppos[:, idxs]
    wf = Float.Wz[:, idxs]
    t = Float.UTC[:, idxs]*86400.
    dwfdt = np.gradient(wf, axis=0)/np.gradient(t, axis=0)

    C_D, alpha_p, alpha_k = res.x

    return calc_w(wf, p, k, V0, alpha_p, p0, alpha_k, k0, M, dwfdt, g, rho, C_D, A)


def apply_w_fit_ud(Float, hpids, res, V0, p0, k0, M, A):

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)

    lat = Float.lat_start[idxs]
    p = Float.P[:, idxs]
    g = -gsw.grav(lat, p)
    rho = Float.rho[:, idxs]
    k = Float.ppos[:, idxs]
    wf = Float.Wz[:, idxs]
    t = Float.UTC[:, idxs]*86400.
    dwfdt = np.gradient(wf, axis=0)/np.gradient(t, axis=0)
    hpids = Float.hpid[idxs]
    up = np.full_like(p, True, dtype=bool)
    up[:, hpids % 2 != 0] = False

    C_D_u, C_D_d, alpha_p, alpha_k = res.x

    return calc_w_ud(wf, p, k, V0, alpha_p, p0, alpha_k, k0, M, dwfdt, g, rho, C_D_u, C_D_d, A, up)

# %% up and down separate
Float = E76

pmin = 50.
pmax = 1400.

hpids_up = np.arange(50, 150, 2)
hpids_down = np.arange(51, 151, 2)
hpids_updown = np.arange(50, 150)

p0 = 2000.
k0 = 16.
A = 0.0214
V0 = 2.62e-2
M = 27.197

res_up = fit(Float, hpids_up, pmin, pmax, V0, p0, k0, M, A)
res_down = fit(Float, hpids_down, pmin, pmax, V0, p0, k0, M, A)
res_updown = fit(Float, hpids_updown, pmin, pmax, V0, p0, k0, M, A)

# %%
hpids_up = np.arange(2, 150, 2)
hpids_down = np.arange(1, 151, 2)

__, idxs_up = Float.get_profiles(hpids_up, ret_idxs=True)
__, idxs_down = Float.get_profiles(hpids_down, ret_idxs=True)

w_up = apply_w_fit(Float, hpids_up, res_up, V0, p0, k0, M, A)
w_down = apply_w_fit(Float, hpids_down, res_down, V0, p0, k0, M, A)

w_up_lin = Float.Ww[:, idxs_up]
w_down_lin = Float.Ww[:, idxs_down]


# %% up down seperated drag coef only
Float = E76

pmin = 50.
pmax = 1400.

hpids_updown = np.arange(50, 150)

p0 = 2000.
k0 = 16.
A = 0.0214
V0 = 2.62e-2
M = 27.197

res_ud = fit_ud(Float, hpids_updown, pmin, pmax, V0, p0, k0, M, A)





