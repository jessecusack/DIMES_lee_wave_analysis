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
from utils import Bunch


default_params = Bunch(
    window='sin2taper',
    dz=4.,
    bin_width=300.,
    bin_overlap=200.,
    plot=False,
    plot_dir='../figures/finescale',
    m_0=1./150.,
    m_c=1./15.
    )


def sin2taper(L):
    """A boxcar window that tapers the last 10% of points of both ends using a
    sin^2 function."""
    win = np.ones(L)
    idx10 = int(np.ceil(L/10.))
    idxs = np.arange(idx10)
    win[:idx10] = np.sin(np.pi*idxs/(2.*idx10))**2
    win[-idx10:] = np.cos(np.pi*(idxs + 1 - L)/(2.*idx10))**2
    return win


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
        win = sig.get_window('boxcar', nfft)
    elif window == 'sin2taper':
        win = sin2taper(nfft)
    else:
        win = sig.get_window(window, nfft)

    if scaling == 'density':
        scale = 1.0/(fs*(win*win).sum())
    elif scaling == 'spectrum':
        scale = 1.0/win.sum()**2
    elif scaling == 'waterman':
        scale = 1.0/(2.*np.pi*len(x))
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

    # Make sure the zero frequency is really zero.
    Pxx[0], Pyy[0], Pxy[0] = 0., 0., 0.

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
    # NOTE that in the MATLAB code QS = -Pxy.imag because they take the
    # conjugate transpose of Pxy first meaning that the imaginary parts are
    # multiplied by -1.
    QS = Pxy.imag
    return (Pxx + Pyy - 2.*QS)/2.


def CCW_ps(Pxx, Pyy, Pxy):
    """Counter clockwise power spectrum."""
    QS = Pxy.imag
    return (Pxx + Pyy + 2*QS)/2.


def integrated_ps(m, P, m_c, m_0):
    """Integrates power spectrum P between wavenumbers m_c and m_0."""

    idxs = (m < m_c) & (m > m_0)
    return np.trapz(P[idxs], x=m[idxs])


def spectral_correction(m, use_range=True, use_diff=True, use_interp=True,
                        use_tilt=True, use_bin=True, dzt=8., dzr=8.,
                        dzfd=8., dzg=8., ddash=5.4, dzs=8.):
    """Spectral corrections for LADCP data - see Polzin et. al. 2002.

    m is vertical wavenumber. Needs factor of 2 pi.

    See Sheen/Shuckburgh MATLAB code for comments.

    dzt: transmitted sound pulse length projected on the vertical
    dzr: receiver processing bin length
    dzfd: first-differencing interval (=dzr)
    dzg: interval of depth grid onto which single-ping piecewise-linear
        continuous profiles of vertical shear are binned
    dzs: superensemble pre-averaging interval, usually chosen to be dzg

    Notes
    -----

    A quadratic fit to the range-maxima (rmax) ? d? pairs given by Polzin et
    al. (2002) yields ddash = -1.2+0.0857r_{max} - 0.000136r_{max}^2 ,
    (5) max max
    which has an intercept near rmax = 14 m. It should be noted that
    expressions (4) and (5) are semi-empirical and apply strictly only to the
    data set of Polzin et al. (2002). Estimating rmax ? 255 m as the range at
    which 80% of all ensembles have valid velocities yields d? ? 11.8 m in case
    of this data set (i.e as in Thurherr 2011 NOT DIMES - nede to updaye!!).
    ddash is determined empirically by Polzin et al. (2002) and is dependent
    on range and the following assumptions:
       + small tilts (~ 3 deg)
       + instrument tilt and orientation are constant over measurement period
       + instrument tilt and orientation are independent
       + tilt attenuation is limited by bin-mapping capabilities of
       RDI (1996) processing

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

    # Binning
    if use_bin:
        C_bin = np.sinc(m*dzg/(2.*np.pi))**2 * np.sinc(m*dzs/(2.*np.pi))**2
    else:
        C_bin = 1.

    return C_range*C_diff*C_interp*C_tilt*C_bin


def window_ps(dz, U, V, dUdz, dVdz, strain, N2_ref, params=default_params):
    """Calculate the power spectra for a window of data."""

    window = params.window
    plot = params.plot

    # Normalise the shear by the mean buoyancy frequency.
    ndUdz = dUdz/np.mean(np.sqrt(N2_ref))
    ndVdz = dVdz/np.mean(np.sqrt(N2_ref))

    # Compute the (co)power spectral density.
    m, PdU, PdV, PdUdV = coperiodogram(ndUdz, ndVdz, fs=1./dz, window=window)
    # We only really want the cospectrum for shear so the next two lines are
    # something of a hack where we ignore unwanted output.
    __, PU, PV, __ = coperiodogram(U, V, fs=1./dz, window=window)
    __, Pstrain, __, __ = coperiodogram(strain, U, fs=1./dz, window=window)

    # Clockwise and counter clockwise spectra.
    PCW = CW_ps(PdU, PdV, PdUdV)
    PCCW = CCW_ps(PdU, PdV, PdUdV)
    # Shear spectra.
    Pshear = PdU + PdV
    # Kinetic energy spectra.
    PEK = (PU + PV)/2.

    # Plotting here generates a crazy number of plots.
    if plot:
        X_names = ['PEk', 'shear', 'strain', 'CW', 'CCW']
        X = [PEK, Pshear, Pstrain, PCW, PCCW]
        for x, name in zip(X, X_names):
            plt.figure()
            plt.loglog(m, x)
            plt.title(name)

    return m, Pshear, Pstrain, PCW, PCCW, PEK


def analyse_profile(Pfl, params=default_params):
    """"""

    plot = params.plot

    # First remove NaN values and interpolate variables onto a regular grid.
    dz = params.dz
    z = np.arange(np.nanmin(Pfl.z), 0., dz)
    U = Pfl.interp(z, 'zef', 'U_abs')
    V = Pfl.interp(z, 'zef', 'V_abs')
    dUdz = Pfl.interp(z, 'zef', 'dUdz')
    dVdz = Pfl.interp(z, 'zef', 'dVdz')
    strain = Pfl.interp(z, 'z', 'strain_z')
    N2_ref = Pfl.interp(z, 'z', 'N2_ref')

    X = [U, V, dUdz, dVdz, strain, N2_ref]

    if plot:
        X_names = ['U', 'V', 'dUdz', 'dVdz', 'strain', 'N2_ref']
        for x, name in zip(X, X_names):
            plt.figure()
            plt.plot(x, z)
            plt.title(name)

    # Split varables into overlapping window segments, bare in mind the last
    # window may not be full.
    width = params.bin_width
    overlap = params.bin_overlap
    wdws = [wdw.window(z, x, width=width, overlap=overlap) for x in X]

    N = wdws[0].shape[0]
    z_mean = np.empty(N)
    EK = np.empty(N)
    R_pol = np.empty(N)
    R_om = np.empty(N)

    for i, w in enumerate(zip(*wdws)):

        # This takes the z values from the horizontal velocity.
        wz = w[0][0]
        z_mean[i] = np.mean(wz)
        # (Terrible code) removes the z values from windowed variables.
        w = [var[1] for var in w]

        # Get the useful power spectra.
        m, PCW, PCCW, Pshear, Pstrain, PEK = window_ps(dz, *w, params=params)

        # Integrate the spectra.
        m_0 = params.m_0
        m_c = params.m_c

        I = [integrated_ps(m, P, m_c, m_0) for P in [Pshear, Pstrain,
                                                     PCW, PCCW, PEK]]
        Ishear, Istrain, ICW, ICCW, IEK = I

        EK[i] = IEK
        R_pol[i] = ICCW/ICW
        R_om[i] = Ishear/Istrain

    return z_mean, EK, R_pol, R_om


def analyse_float(Float, hpids):
    """"""
    # Nothing special for now.
    __, idxs = Float.get_profiles(hpids, ret_idxs=True)
    return [analyse_profile(Pfl) for Pfl in Float.Profiles[idxs]]
