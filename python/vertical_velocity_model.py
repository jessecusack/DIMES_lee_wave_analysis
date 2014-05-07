# -*- coding: utf-8 -*-
"""
Created on Tue Apr 08 14:33:17 2014

@author: jc3e13
"""

import pylab as pl
import emapex
import scipy.optimize as op


def still_water_model(params, data):
    """Calculates and returns the vertical velocity that the glider would have
    if it were moving in still water."""

    ppos, p, rho = data
    # Known parameters.
    g = 9.81  # Gravitational acceleration [m s-2].
    M = 27.179  # Float mass [kg].
    alpha_ppos = 1.156e-6  # Coeff. of expansion with piston position [m3].
    ppos_0 = 16.  # Piston position when neutrally buoyant [-].

    # Unknown parameters.
    V_0 = 0.  # Float volume when neutrally buoyant [m3]
    # Combination of drag and surface area.
    CA = 0.

    V_0, CA, alpha_p, p_0 = params

    alpha_p = 3.76e-6  # Coeff. of expansion with pressure [-].
    p_0 = 2000.  # Pressure when float is neutrally buoyant [dbar].

    float_volume = V_0 + alpha_ppos*(ppos - ppos_0)/V_0 - alpha_p*(p - p_0)
    effective_water_volume = M/rho

    # Need to keep track of whether the float is going up or down but also not
    # cause problems with the square root.
    is_going_down = float_volume < effective_water_volume

    w = pl.sqrt(g/(CA)*pl.absolute(float_volume - effective_water_volume))
    w[is_going_down] = -1.*w[is_going_down]

    return w


def still_water_model_2(params, data):
    """"""
    a, b, c, d = params
    ppos, p, rho = data

    w_sqrd = a + b*ppos + c*p + d*rho

    is_going_down = w_sqrd < 0.

    w = pl.sqrt(pl.absolute(w_sqrd))
    w[is_going_down] = -1.*w[is_going_down]

    return w


def cost_func(params, data, w_float_sqrd):
    """"""
    ws = still_water_model_2(params, data)

    cost = pl.mean(ws**2 - w_float_sqrd, axis=0)

    return cost


def fit_model(Float, hpids):
    """"""
    P_vals = pl.arange(100., 1400., 8)
    __, Pg, rhog = Float.get_interp_grid(hpids, P_vals, 'P', 'rho')
    __, __, pposg = Float.get_interp_grid(hpids, P_vals, 'P', 'ppos')
    __, __, w_floatg = Float.get_interp_grid(hpids, P_vals, 'P', 'Wpef')
    __, idxs = Float.get_profiles(hpids, ret_idxs=True)
    existing_hpids = Float.hpid[idxs]

    w_floatg_sqrd = w_floatg**2

    up_idxs = emapex.up_down_indices(existing_hpids, 'up')
    down_idxs = emapex.up_down_indices(existing_hpids, 'down')

#    # Why is this necessary?
#    w_floatg_sqrd[:, down_idxs] = -w_floatg_sqrd[:, down_idxs]
    data = (pposg, Pg, rhog)

    params0 = pl.array([1., 1., 1., 1.])
    params = op.leastsq(cost_func, params0,
                        args=(data, w_floatg_sqrd))[0]

#    params = results.x
#    print(results.message)
    print("The final parameter values:")
    print(params)


    return params

#    Wref(i) = sqrt(x(i,4)/dens0);
#    alpha(i) = -x(i,2)*dens0/x(i,4);
#    kappa(i) = x(i,3)*dens0/x(i,4);
#    k0(i) = -x(i,1)/x(i,2) - x(i,4)/(x(i,2)*dens0);
#
