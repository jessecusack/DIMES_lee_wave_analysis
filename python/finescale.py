# -*- coding: utf-8 -*-
"""
Created on Tue May 20 15:45:36 2014

A place for finescale parameterisation functions.

@author: jc3e13
"""

import numpy as np
import gsw
import scipy.signal as sig


def adiabatic_level(P, SA, T, lat, P_bin_width=200., deg=1):
    """Generate smooth buoyancy frequency profile by applying the adiabatic
    levelling method of Bray and Fofonoff (1981).

      Parameters
      ----------
      P : 1-D ndarray
          Pressure [dbar]
      SA : 1-D ndarray
          Absolute salinity [g/kg]
      T : 1-D ndarray
          Temperature [degrees C]
      lat : float
          Latitude [-90...+90]
      p_bin_width : float, optional
          Pressure bin width [dbar]
      deg : int, option
          Degree of polynomial fit.

      Returns
      -------
      N2_ref : 1-D ndarray
          Reference buoyancy frequency [s-2]

      Raises
      ------

      Notes
      -----
      Calls to the gibbs seawater toolbox are slow and therefore this function
      is quite slow.

      Examples
      --------

    """

    N2_ref = np.NaN*P.copy()
    nans = np.isnan(P) | np.isnan(SA) | np.isnan(T)

    # If there are nothing but NaN values don't waste time.
    if np.sum(nans) == nans.size:
        return N2_ref

    P = P[~nans]
    SA = SA[~nans]
    T = T[~nans]
    P_min, P_max = np.min(P), np.max(P)

    shape = (P.size, P.size)
    Pm = np.NaN*np.empty(shape)
    SAm = np.NaN*np.empty(shape)
    Tm = np.NaN*np.empty(shape)

    # Populate bins.
    for i in xrange(len(P)):

        P_bin_min = np.maximum(P[i] - P_bin_width/2., P_min)
        P_bin_max = np.minimum(P[i] + P_bin_width/2., P_max)
        in_bin = np.where((P >= P_bin_min) & (P <= P_bin_max))[0]

        Pm[in_bin, i] = P[in_bin]
        SAm[in_bin, i] = SA[in_bin]
        Tm[in_bin, i] = T[in_bin]

    P_bar = np.nanmean(Pm, axis=0)
    T_bar = np.nanmean(Tm, axis=0)
    SA_bar = np.nanmean(SAm, axis=0)

    # Perform thermodynamics once only...
    rho_bar = gsw.pot_rho_t_exact(SA_bar, T_bar, P_bar, P_bar)
    sv = 1./gsw.pot_rho_t_exact(SAm, Tm, Pm, P_bar)

    p = []
    for P_bin, sv_bin in zip(Pm.T, sv.T):
        bnans = np.isnan(P_bin)
        p.append(np.polyfit(P_bin[~bnans],
                            sv_bin[~bnans] - np.nanmean(sv_bin),
                            deg))

    p = np.array(p)

    g = gsw.grav(lat, P_bar)
    # The factor 1e-4 is needed for conversion from dbar to Pa.
    N2_ref[~nans] = -1e-4*rho_bar**2*g**2*p[:, 0]

    return N2_ref


def apply_strain(Float, P_bin_width=100.):
    """
    Takes a float and a pressure bin width then calculates N2_ref for all
    available profiles and sets N2_ref and strain_z as attributes of the
    float class.

    """

    Pg = getattr(Float, 'P')
    SAg = getattr(Float, 'SA')
    Tg = getattr(Float, 'T')
    lats = getattr(Float, 'lat_start')

    N2_ref = np.NaN*Pg.copy()

    for i, (P, SA, T, lat) in enumerate(zip(Pg.T, SAg.T, Tg.T, lats)):
        print("hpid: {}".format(Float.hpid[i]))
        N2_ref[:, i] = adiabatic_level(P, SA, T, lat, P_bin_width)

    N2 = getattr(Float, 'N2')
    strain_z = (N2 - N2_ref)/N2_ref

    setattr(Float, 'N2_ref', N2_ref)
    setattr(Float, 'strain_z', strain_z)


def h_gregg(R=3.):
    """Gregg 2003 implimentation."""
    return 3.*(R + 1)/(2.*R*np.sqrt(2*(R - 1)))


def h_whalen(R=3.):
    """Whalen 2012 implimentation based on Kunze 2006 implimentation (which
    is based on Gregg yet equations are different)."""
    return R*(R + 1)/(6.*np.sqrt(2*(R - 1)))


def L(f, N):
    f30 = 7.292115e-5  # rad s-1
    N0 = 5.2e-3  # rad s-1
    return f*np.arccosh(N/f)/(f30*np.arccosh(N0/f30))


def integrated_variance(x, fs, m_c, m_0):
    """Integrated variance of quantity x."""

    # Calculate the power spectral density.
    m, spec = sig.periodogram(x, fs=fs)


#def GM(omega, f, N):
#    """Garrett-Munk (power?) spectrum.
#
#    Model variant for the frequency spectrum of horizontal kinetic energy.
#    """
#
#    om2 = omega**2.
#    f2 = f**2.
#    N2 = N**2.
#
#    a = 0.5
#    b = 0.
#
#    continuum = f*(om2 + f2)*(N2 - om2)/(N2*om2**2.)
#    inertial = om2/(om2 - f2)**a
#    buoyancy = N2/(N2 - om2)**b
#
#    return continuum*inertial*buoyancy
#
#def GM2():
#    pass
# betastar=pi*jstar/b*sqrt(N2mean)/N0;
# phi_zeta = E*b^3*N0^2/(2*jstar*pi*N2mean) ./(1+kzax/betastar).^2; % power spectral density of vertical displacement
# phi_eta = kzax.^2.*phi_zeta; % power spectral density of vertical strain
