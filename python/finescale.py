# -*- coding: utf-8 -*-
"""
Created on Tue May 20 15:45:36 2014

A place for finescale parameterisation functions.

@author: jc3e13
"""

import numpy as np
import scipy as sp
import gsw


def adiabatic_level(P, SA, T, N2, lat, P_bin_width=200.):
    """Generate smooth buoyancy frequency profile by applying the adiabatic
    levelling method of Bray and Fofonoff (1981).

    TODO: Add optional parameter ignore_above.

      Parameters
      ----------
      P : 1-D ndarray
          Pressure [dbar]
      SA : 1-D ndarray
          Absolute salinity [g/kg]
      T : 1-D ndarray
          Temperature [degrees C]
      N2 : 1-D ndarray
          Buoyancy frequency [s-2]
      lat : float
          The [-90...+90]
      p_bin_width : float, optional
          Pressure bin width [dbar]

      Returns
      -------
      N2_ref : 1-D ndarray of floats
          Reference buoyancy frequency [s-2]

      Raises
      ------

      Notes
      -----

      Examples
      --------

    """

    N2_ref = np.NaN*P.copy()


    for i in xrange(len(P)):
        P_min = np.maximum(P[i] - P_bin_width/2., np.nanmin(P))
        P_max = np.minimum(P[i] + P_bin_width/2., np.nanmax(P))
        in_bin = np.where((P >= P_min) & (P <= P_max))[0]
        if in_bin.size == 0:
            continue
        P_bar = np.nanmean(P[in_bin])
        T_bar = np.nanmean(T[in_bin])
        SA_bar = np.nanmean(SA[in_bin])
        rho_bar = gsw.pot_rho_t_exact(SA_bar, T_bar, P_bar, P_bar)
        sv = 1./gsw.pot_rho_t_exact(SA[in_bin], T[in_bin], P[in_bin], P_bar)

        deg = 1

        p = np.polyfit(P[in_bin], sv - np.nanmean(sv), deg)
        g = gsw.grav(lat, P_bar)
        N2_ref[i] = -1e-4 * rho_bar**2 * g**2 * p[0]

    return N2_ref
