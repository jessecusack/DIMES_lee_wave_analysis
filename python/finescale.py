# -*- coding: utf-8 -*-
"""
Created on Tue May 20 15:45:36 2014

A place for finescale parameterisation functions.

@author: jc3e13
"""

import numpy as np
import gsw
import scipy.signal as sig
import matplotlib.pyplot as plt
import window as wdw


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

    Notes
    -----
    Calls to the gibbs seawater toolbox are slow and therefore this function
    is quite slow.

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


def coperiodogram(x, y, fs=1.0, window=None, nfft=None, detrend='linear',
                  scaling='density'):
    """
    Estimate co-power spectral density using periodogram method.

    Parameters
    ----------
    x : array_like
        Measurement values
    y : array_like
        Measurement values
    fs : float, optional
        Sampling frequency of `x` and `y`. e.g. If the measurements are a time
        series then `fs` in units of Hz. Defaults to 1.0.
    window : str or tuple or array_like, optional
        Desired window to use. See `get_window` for a list of windows and
        required parameters. If `window` is array_like it will be used
        directly as the window and its length will be used for nperseg.
        Defaults to `boxcar`. New window option `sin2taper` is available.
    nfft : int, optional
        Length of the FFT used, if a zero padded FFT is desired.  If None,
        the FFT length is `len(x)`. Defaults to None.
    detrend : str or function, optional
        Specifies how to detrend each segment. If `detrend` is a string,
        it is passed as the ``type`` argument to `detrend`. If it is a
        function, it takes a segment and returns a detrended segment.
        Defaults to 'linear'.
    scaling : { 'density', 'spectrum' }, optional
        Selects between computing the power spectral density ('density')
        where Pxx has units of V**2/Hz if x is measured in V and computing
        the power spectrum ('spectrum') where Pxx has units of V**2 if x is
        measured in V. Defaults to 'density'.

    Returns
    -------
    f : ndarray
        Array of sample frequencies.
    Pxx : ndarray
        Power spectral density or power spectrum of x.
    Pyy : ndarray
        Power spectral density or power spectrum of y.
    Pxy : ndarray (complex)
        Co-power spectral density or co-power spectrum of x.

    """

    x, y = np.asarray(x), np.asarray(y)

    if (x.size == 0) or (y.size == 0):
        raise ValueError('At least one of the inputs is empty.')

    if len(x) != len(y):
        raise ValueError('x and y must have the same length.')

    if nfft is None:
        nfft = len(x)

    if window is None:
        window = 'boxcar'
        win = sig.get_window(window, nfft)
    elif window == 'sin2taper':
        win = np.ones(nfft)
        idx10 = int(np.ceil(nfft/10.))
        idxs = np.arange(idx10)
        win[:idx10] = np.sin(np.pi*idxs/(2.*idx10))**2
        win[-idx10:] = np.cos(np.pi*(idxs - nfft)/(2.*idx10))**2
    else:
        win = sig.get_window(window, nfft)

    if scaling == 'density':
        scale = 1.0/(fs*(win*win).sum())
    elif scaling == 'spectrum':
        scale = 1.0/win.sum()**2
    else:
        raise ValueError('Unknown scaling: %r' % scaling)

    x_dt = sig.detrend(x, type=detrend)
    y_dt = sig.detrend(y, type=detrend)

    xft = np.fft.fft(win*x_dt, nfft)
    yft = np.fft.fft(win*y_dt, nfft)

    # Power spectral density in x, y and xy.
    Pxx = (xft*xft.conj()).real
    Pyy = (yft*yft.conj()).real
    Pxy = (xft*yft.conj())

    M = nfft/2 + 1

    # Chop spectrum in half.
    Pxx, Pyy, Pxy = Pxx[:M], Pyy[:M], Pxy[:M]

    # Multiply spectrum by 2 except for the Nyquist and constant elements to
    # account for the loss of negative frequencies.
    Pxx[..., 1:-1] *= 2*scale
    Pxx[..., (0,-1)] *= scale

    Pyy[..., 1:-1] *= 2*scale
    Pyy[..., (0,-1)] *= scale

    Pxy[..., 1:-1] *= 2*scale
    Pxy[..., (0,-1)] *= scale

    f = np.arange(Pxx.shape[-1])*(fs/nfft)

    return f, Pxx, Pyy, Pxy


def CW_ps(Pxx, Pyy, Pxy):
    """Clockwise power spectrum."""
    QS = -Pxy.imag
    return (Pxx + Pyy - 2.*QS)/2.


def CCW_ps(Pxx, Pyy, Pxy):
    """Counter clockwise power spectrum."""
    QS = -Pxy.imag
    return (Pxx + Pyy + 2*QS)/2.


def integrated_ps(m, P, m_c, m_0):
    """Integrates power spectrum P between wavenumbers m_c and m_0."""

    idxs = (m < m_c) & (m > m_0)
    return np.trapz(P[idxs], x=m[idxs])


def spectral_correction(m, use_range=True, use_diff=True, use_interp=True,
                        use_tilt=True, dzt=16., dzr=16., dzfd=16., dzg=20.,
                        ddash=9.):
    """Spectral corrections for LADCP data - see Polzin et. al. 2002.

    m is vertical wavenumber. Needs factor of 2 pi.

    See Sheen/Shuckburgh MATLAB code for comments.
    """

    # Range averaging.
    if use_range:
        C_range = np.sinc(m*dzt/(2.*np.pi))**2 * np.sinc(m*dzr/(2.*np.pi))**2
    else:
        C_range = 1.

    # First differencing.
    if use_diff:
        C_diff = np.sinc(m*dzfd/(2.*np.pi))**2
    else:
        C_diff = 1.

    # Interpolation.
    if use_interp:
        C_interp = np.sinc(m*dzr/(2.*np.pi))**4 * np.sinc(m*dzg/(2.*np.pi))**2
    else:
        C_interp = 1.

    # Tilting.
    if use_tilt:
        C_tilt = np.sinc(m*ddash/(2.*np.pi))**2
    else:
        C_tilt = 1.

    return C_range*C_diff*C_interp*C_tilt


def analyse_profile(Pfl, plot=False):

    # First remove NaN values and interpolate variables onto a regular grid.
    dz = 4.  # Depth [m]
    z = np.arange(np.nanmin(Pfl.z), 0., dz)
    U = Pfl.interp(z, 'zef', 'U_abs')
    V = Pfl.interp(z, 'zef', 'V_abs')
    dUdz = Pfl.interp(z, 'zef', 'dUdz')
    dVdz = Pfl.interp(z, 'zef', 'dVdz')
    strain = Pfl.interp(z, 'z', 'strain_z')
    N2_ref = Pfl.interp(z, 'z', 'N2_ref')

    if plot:
        X_names = ['U', 'V', 'dUdz', 'dVdz', 'strain', 'N2_ref']
        X = [U, V, dUdz, dVdz, strain, N2_ref]
        for x, name in zip(X, X_names):
            plt.figure()
            plt.plot(x, z)
            plt.title(name)

    # Split varables into overlapping window segments, bare in mind the last
    # window may not be full.
    width = 300.
    overlap = 200.
    WDW = [wdw.window(z, x, width=width, overlap=overlap) for x in X]

    N = WDW[0].shape[0]
    z_mean = np.empty(N)
    EK = np.empty(N)
    Rpol = np.empty(N)
    Rom = np.empty(N)

    # Next normalise shear by mean N and then get spectra for everything. The
    # for loop goes over each segment.
    for i, (wU, wV, wdUdz, wdVdz, wstrain, wN2_ref) in enumerate(zip(*WDW)):
        # Bad code I know but simplyfy what variables contain by removing the
        # z values.
        wz = wU[0]
        z_mean[i] = np.mean(wz)
        wU, wV, wdUdz, wdVdz, wstrain, wN2_ref = \
            wU[1], wV[1], wdUdz[1], wdVdz[1], wstrain[1], wN2_ref[1]

        # Normalise the shear by the mean buoyancy frequency.
        ndUdz = wdUdz/np.mean(np.sqrt(wN2_ref))
        ndVdz = wdVdz/np.mean(np.sqrt(wN2_ref))

        # Compute the (co)power spectral density.
        m, PdU, PdV, PdUdV = coperiodogram(ndUdz, ndVdz, fs=1./dz)
        __, PU, PV, __ = coperiodogram(wU, wV, fs=1./dz)
        __, Pstrain = sig.periodogram(wstrain, fs=1./dz, detrend='linear')
        PCW = CW_ps(PdU, PdV, PdUdV)
        PCCW = CCW_ps(PdU, PdV, PdUdV)
        Pshear = PdU + PdV

        # Plotting here generates a crazy number of plots.
#        if plot:
#            X_names = ['U', 'V', 'shear', 'strain', 'CW', 'CCW']
#            X = [PU, PV, Pshear, Pstrain, PCW, PCCW]
#            for x, name in zip(X, X_names):
#                plt.figure()
#                plt.loglog(m, x)
#                plt.title(name)

        # Integrate the spectra.
        m_0 = 1./200.
        m_c = 1./10.

        I = [integrated_ps(m, P, m_c, m_0) for P in [PU, PV, Pshear, Pstrain,
                                                     PCW, PCCW]]
        IU, IV, Ishear, Istrain, ICW, ICCW = I

        EK[i] = (IU + IV)/2.
        Rpol[i] = ICCW/ICW
        Rom[i] = Ishear/Istrain