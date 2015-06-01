# -*- coding: utf-8 -*-
"""
Created on Tue Apr 08 14:33:17 2014

@author: jc3e13
"""

import emapex
import scipy.optimize as op
import numpy as np
from utils import Bunch


def fitter(Float, params0, fixed, model='1', hpids=None,
           profiles='all', cf_key='sqdiff', P_vals=None,
           data_names=['ppos', 'P', 'rho']):
    """This function takes an EM-APEX float, fits a vertical velocity model
    using the given arguments, estimates errors using bootstrapping technique,
    and then passes this information to the float.

    """
    # Argument checks.

    # Defaults
    if hpids is None:
        hpids = np.arange(50, 151)

    if P_vals is None:
        P_vals = np.arange(60., 1500., 8)

    print('Fitting model '+model+'.')

    if model == '1':
        still_water_model = still_water_model_1
    elif model == '2':
        still_water_model = still_water_model_2
    else:
        raise ValueError('Cannot find model.')

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)
    hpids = Float.hpid[idxs]

    if profiles == 'all':

        # Fitting.
        data = [Float.get_interp_grid(hpids, P_vals, 'P', data_name)[2]
                for data_name in data_names]

        __, __, w_f = Float.get_interp_grid(hpids, P_vals, 'P', 'Wz')

        cargs = (fixed, still_water_model, w_f, data, cf_key)

        p, pcov, info, mesg, ier = op.leastsq(cost, params0, args=cargs,
                                              full_output=True)

        # Bootstrapping.
        print('Starting bootstrap.')
        ps = []
        # 200 random data sets are generated and fitted
        for i in range(200):
            rand_idx = np.random.rand(*w_f.shape) < 0.25
            rand_data = [d[rand_idx] for d in data]
            cargs = (fixed, still_water_model, w_f[rand_idx], rand_data,
                     cf_key)
            rand_params, __ = op.leastsq(cost, params0, args=cargs)
            ps.append(rand_params)

        ps = np.array(ps)
        pcov = np.cov(ps.T)
        pmean = np.mean(ps, 0)
        pcorr = np.corrcoef(ps.T)

    elif profiles == 'updown':

        up_idxs = emapex.up_down_indices(hpids, 'up')
        up_hpids = hpids[up_idxs]
        down_idxs = emapex.up_down_indices(hpids, 'down')
        down_hpids = hpids[down_idxs]

        # We now contain variables in lists.
        p = 2*[0]
        ps = 2*[0]
        pmean = 2*[0]
        pcov = 2*[0]
        pcorr = 2*[0]
        info = 2*[0]
        mesg = 2*[0]
        ier = 2*[0]

        # This is a bit of hack since profiles gets changed here... bad.
        for i, _hpids in enumerate([up_hpids, down_hpids]):

            # Fitting.
            data = [Float.get_interp_grid(_hpids, P_vals, 'P', data_name)[2]
                    for data_name in data_names]

            __, __, w_f = Float.get_interp_grid(_hpids, P_vals, 'P', 'Wz')

            cargs = (fixed, still_water_model, w_f, data, cf_key)

            p[i], pcov[i], info[i], mesg[i], ier[i] = \
                op.leastsq(cost, params0, args=cargs, full_output=True)

            # Bootstrapping.
            print('Starting bootstrap.')

            bps = []
            # 100 random data sets are generated and fitted
            for __ in range(100):
                rand_idx = np.random.rand(*w_f.shape) < 0.25
                rand_data = [d[rand_idx] for d in data]
                cargs = (fixed, still_water_model, w_f[rand_idx], rand_data,
                         cf_key)
                rand_params, __ = op.leastsq(cost, params0, args=cargs)
                bps.append(rand_params)

            ps[i] = np.array(bps)
            pcov[i] = np.cov(ps[i].T)
            pmean[i] = np.mean(ps[i], 0)
            pcorr[i] = np.corrcoef(ps[i].T)

    else:
        raise ValueError("profiles can be 'all' or 'updown'")

    wfi = Bunch(params0=params0,
                fixed=fixed,
                model_func=still_water_model,
                hpids=hpids,
                profiles=profiles,
                cf_key=cf_key,
                P_vals=P_vals,
                data_names=data_names,
                p=p,
                ps=ps,
                pmean=pmean,
                pcov=pcov,
                pcorr=pcorr,
                info=info,
                mesg=mesg,
                ier=ier)

    return wfi


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
        Key to select cost function, either 'sqdiff' or 'diff_square'

    Returns
    -------
    c : numpy.ndarray
        The cost calculated from cost_func.

    Notes
    -----
    Uses the Profile.interp function.


    """

    ws = model(params, data, fixed)

    if cf_key == 'sqdiff':
        c = ws**2 - wf**2
    elif cf_key == 'diff_square':
        c = (ws - wf)**2
    else:
        raise ValueError('Incorrect cf_key')

    return c.flatten()
