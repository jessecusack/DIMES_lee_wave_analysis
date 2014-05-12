# -*- coding: utf-8 -*-
"""
Created on Tue Apr 08 14:33:17 2014

@author: jc3e13
"""

import emapex
import scipy.optimize as op
import numpy as np
import matplotlib.pyplot as plt


def fitter(Float, params0, model='1', hpids=np.arange(50, 151), profiles='all',
           cf_key='sqdiff', P_vals=np.arange(100., 1400., 8),
           data_names=['ppos', 'P', 'rho']):
    """This function takes an EM-APEX float, fits a vertical velocity model and
    applies it to the float using the given scope and saves the scope and
    parameters for future use."""

    if model == '1':
        still_water_model = still_water_model_1
    elif model == '2':
        still_water_model = still_water_model_2

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)
    existing_hpids = Float.hpid[idxs]

    if profiles == 'all':
        hpids = existing_hpids
    elif profiles == 'up':
        up_idxs = emapex.up_down_indices(existing_hpids, 'up')
        hpids = existing_hpids[up_idxs]
    elif profiles == 'down':
        down_idxs = emapex.up_down_indices(existing_hpids, 'down')
        hpids = existing_hpids[down_idxs]

    data = [Float.get_interp_grid(hpids, P_vals, 'P', data_name)[2]
            for data_name in data_names]

    __, __, w_floatg = Float.get_interp_grid(hpids, P_vals, 'P', 'Wz')

    cargs = (still_water_model, w_floatg, data, cf_key)

    params = op.leastsq(cost, params0, args=cargs)[0]

    print("The final parameter values:")
    print(params)

    data = [getattr(Float, data_name) for data_name in data_names]

    setattr(Float, 'rWs', still_water_model(params, data))
    setattr(Float, 'rWw', Float.rWz - Float.rWs)

    Float.update_profiles()

    t1, Ws = Float.get_timeseries(hpids, 'rWs')
    __, Wz = Float.get_timeseries(hpids, 'rWz')
    __, Wf = Float.get_timeseries(hpids, 'rWf')
    __, Ww = Float.get_timeseries(hpids, 'rWw')

    plt.figure()
    plt.plot(t1, Ws, 'r', t1, Wz, 'g', t1, Wf, 'r--', t1, Ww, 'b')


def still_water_model_1(params, data):
    """Calculates and returns the vertical velocity that the glider would have
    if it were moving in still water."""
    # TODO: introduce a way of fixing certain parameters...
    ppos, p, rho = data
    # Known parameters. I guess they should be put into 'data'.
    g = 9.81  # Gravitational acceleration [m s-2].
    M = 27.179  # Float mass [kg].
    alpha_ppos = 1.156e-6  # Coeff. of expansion with piston position [m3].
    ppos_0 = 16.  # Piston position when neutrally buoyant [-].

    # Unknown parameters.
    V_0 = 0.  # Float volume when neutrally buoyant [m3]
    CA = 0.  # Combination of drag and surface area.
    alpha_p = 3.76e-6  # Coeff. of expansion with pressure [-].
    p_0 = 2000.  # Pressure when float is neutrally buoyant [dbar].

    V_0, CA, alpha_p, p_0 = params

    float_volume = V_0*(1 + alpha_ppos*(ppos - ppos_0)/V_0 - alpha_p*(p - p_0))
    effective_water_volume = M/rho

    # Need to keep track of whether the float is going up or down but also not
    # cause problems with the square root.
    is_going_down = float_volume < effective_water_volume

    w = np.sqrt(g/(CA)*np.abs(float_volume - effective_water_volume))
    w[is_going_down] = -1.*w[is_going_down]

    return w


def still_water_model_2(params, data):
    """"""
    a, b, c, d = params
    ppos, p, rho = data

    w_sqrd = a + b*ppos + c*p + d*rho

    is_going_down = w_sqrd < 0.

    w = np.sqrt(np.abs(w_sqrd))
    w[is_going_down] = -1.*w[is_going_down]

    return w


def cost(params, model, wf, data, cf_key):
    """The cost function should be minimised when the model parameters are
    optimised.

      Parameters
      ----------
      params : 1-D numpy.ndarray.
          The profile ID numbers at for which to construct grid.
      model : function.
          A model of the float in still water.
      wf : numpy.ndarray
          The float absolute velocity.
      model_args : tuple.
          Additional arguments to model, (model(params, *model_args)).
      cf_key : string.
          Dictionary key to select cost function.

      Returns
      -------
      c : numpy.ndarray
          The cost calculated from cost_func.

      Notes
      -----
      Uses the Profile.interp function.


    """

    def square_diff(ws, wf):
        return ws**2 - wf**2

    cfd = {'sqdiff': square_diff}

    cost_func = cfd[cf_key]

    ws = model(params, data)

    c = cost_func(ws, wf)
    print c, c.shape

    return c




#    Wref(i) = sqrt(x(i,4)/dens0);
#    alpha(i) = -x(i,2)*dens0/x(i,4);
#    kappa(i) = x(i,3)*dens0/x(i,4);
#    k0(i) = -x(i,1)/x(i,2) - x(i,4)/(x(i,2)*dens0);
#
