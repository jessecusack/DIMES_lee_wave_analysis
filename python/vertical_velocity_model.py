# -*- coding: utf-8 -*-
"""
Created on Tue Apr 08 14:33:17 2014

@author: jc3e13
"""

import emapex
import scipy.optimize as op
import numpy as np


class vv_fit_info:
    """Storage class for information relating to the vertical velocity
    fitting."""
    def __init__(self, params0, fixed, model_func, hpids, profiles, cf_key,
                 P_vals, data_names, params, cov, info, mesg, ier):

        self.params0 = params0
        self.fixed = fixed
        self.model_func = model_func
        self.hpids = hpids
        self.profiles = profiles
        self.cf_key = cf_key
        self.P_vals = P_vals
        self.data_names = data_names
        self.params = params
        self.cov = cov
        self.info = info
        self.mesg = mesg
        self.ier = ier


def fitter(Float, params0, fixed, model='1', hpids=np.arange(50, 151),
           profiles='all', cf_key='sqdiff', P_vals=np.arange(60., 1500., 8),
           data_names=['ppos', 'P', 'rho']):
    """This function takes an EM-APEX float, fits a vertical velocity model and
    applies it to the float using the given scope and saves the scope and
    parameters for future use as an attribute of the float class."""

    if model == '1':
        still_water_model = still_water_model_1
    elif model == '2':
        still_water_model = still_water_model_2
    else:
        raise ValueError('Cannot find model.')

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
    else:
        raise ValueError('Cannot understand what type of profiles')

    data = [Float.get_interp_grid(hpids, P_vals, 'P', data_name)[2]
            for data_name in data_names]

    __, __, w_floatg = Float.get_interp_grid(hpids, P_vals, 'P', 'Wz')

    cargs = (fixed, still_water_model, w_floatg, data, cf_key)

    params, cov, info, mesg, ier = op.leastsq(cost, params0, args=cargs,
                                              full_output=True)

    print("The final parameter values:")
    print(params)

    data = [getattr(Float, 'r' + data_name) for data_name in data_names]

    setattr(Float, 'rWs', still_water_model(params, data, fixed))
    setattr(Float, 'rWw', Float.rWz - Float.rWs)

    vfi = vv_fit_info(params0, fixed, still_water_model, hpids, profiles,
                      cf_key, P_vals, data_names, params, cov, info, mesg, ier)

    setattr(Float, '__vfi', vfi)  # Don't want this added to Profile classes.

    Float.update_profiles()

#    t1, Ws = Float.get_timeseries(hpids, 'rWs')
#    __, Wz = Float.get_timeseries(hpids, 'rWz')
#    __, Wf = Float.get_timeseries(hpids, 'rWf')
#    __, Ww = Float.get_timeseries(hpids, 'rWw')
#
#    plt.figure()
#    plt.plot(t1, Ws, 'r', t1, Wz, 'g', t1, Wf, 'r--', t1, Ww, 'b')


def still_water_model_1(params, data, fixed):
    """Calculates and returns the vertical velocity that the glider would have
    if it were moving in still water.

    params:

    0: V_0 = 1.  # Float volume when neutrally buoyant [m3]
    1: CA = 1.  # Combination of drag and surface area.
    2: alpha_p = 3.76e-6  # Coeff. of expansion with pressure [-].
    3: p_0 = 2000.  # Pressure when float is neutrally buoyant [dbar].
    4: alpha_ppos = 1.156e-6  # Coeff. of expansion with piston position [m3].
    5: ppos_0 = 16.  # Piston position when neutrally buoyant [-].
    6: M = 27.179  # Float mass [kg].

    data:

    ppos, p, rho

    fixed:

    List of values to fix with the same numbering as parameters. Use None for
    varying parameters.

    Gravity is given a value of 9.8 m s-2.

    """

    ppos, p, rho = data

    g = 9.8  # Gravitational acceleration [m s-2].

    P = params

    for i, val in enumerate(fixed):
        if val is not None:
            P[i] = val

    float_volume = P[0]*(1 + P[4]*(ppos - P[5])/P[0] - P[2]*(p - P[3]))
    effective_water_volume = P[6]/rho

    # Need to keep track of whether the float is going up or down but also not
    # cause problems with the square root.
    is_going_down = float_volume < effective_water_volume

    w = np.sqrt(g/(P[1])*np.abs(float_volume - effective_water_volume))
    w[is_going_down] = -1.*w[is_going_down]

    return w


def still_water_model_2(params, data, fixed):
    """Currently cannot fix parameters."""
    a, b, c, d = params
    ppos, p, rho = data

    w_sqrd = a + b*ppos + c*p + d/rho

    # Not sure about this condition...
    is_going_down = w_sqrd < 0.

    w = np.sqrt(np.abs(w_sqrd))
    w[is_going_down] = -1.*w[is_going_down]

    return w

#    Wref(i) = sqrt(x(i,4)/dens0);
#    alpha(i) = -x(i,2)*dens0/x(i,4);
#    kappa(i) = x(i,3)*dens0/x(i,4);
#    k0(i) = -x(i,1)/x(i,2) - x(i,4)/(x(i,2)*dens0);
#


def cost(params, fixed, model, wf, data, cf_key):
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

    def diff_square(ws, wf):
        return (ws - wf)**2

    cfd = {'sqdiff': square_diff,
           'diffsq': diff_square}

    cost_func = cfd[cf_key]

    ws = model(params, data, fixed)

    c = cost_func(ws, wf).flatten()

    return c
