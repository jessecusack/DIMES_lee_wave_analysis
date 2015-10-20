# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 11:40:47 2015

@author: jc3e13
"""

import numpy as np
import scipy as sp
import scipy.signal as sig
import matplotlib.pyplot as plt
import os
import sys

import gsw

lib_path = os.path.abspath('../modules')
if lib_path not in sys.path:
    sys.path.append(lib_path)

lib_path = os.path.abspath('../../ocean-tools')
if lib_path not in sys.path:
    sys.path.append(lib_path)

import finescale as fs
import window as wdw
import GM79

reload(fs)


default_corrections = {
    'use_range': False,
    'use_diff': False,
    'use_interp': True,
    'use_tilt': False,
    'use_bin': False,
    'use_volt': True,
    'dzt': 8.,
    'dzr': 8.,
    'dzfd': 8.,
    'dzg': 8.,
    'ddash': 5.4,
    'dzs': 8.,
    'vfi': 50.,
    'mfr': 0.12
    }


default_periodogram_params = {
    'window': 'sin2taper',
    'nfft': 256,
    'detrend': 'linear',
    'scaling': 'density',
    }


default_params = {
    'dz': 4.,
    'zmin': None,
    'zmax': -300.,
    'bin_width': 300.,
    'bin_overlap': 200.,
    'fine_grid_spectra': False,
    'print_diagnostics': False,
    'plot_profiles': False,
    'plot_spectra': False,
    'plot_results': False,
    'plot_dir': '../figures/finescale',
    'm_0': 1./150.,
    'm_c': 1./15.,
    'apply_corrections': False,
    'corrections': default_corrections,
    'periodogram_params': default_periodogram_params,
    'mixing_efficiency': 0.2
    }


def sin2taper(L):
    """A boxcar window that tapers the last 10% of points of both ends using a
    sin^2 function."""
    win = np.ones(L)
    idx10 = int(np.ceil(L/10.))
    idxs = np.arange(idx10)
    win[:idx10] = np.sin(np.pi*idxs/(2.*idx10))**2
    win[-idx10:] = np.cos(np.pi*(idxs + 1 - L)/(2.*idx10))**2
    return win


def h_gregg(R=3.):
    """Gregg 2003 implimentation."""
    return 3.*(R + 1)/(2.*R*np.sqrt(2*np.abs(R - 1)))


def h_whalen(R=3.):
    """Whalen 2012 implimentation based on Kunze 2006 implimentation (which
    is based on Gregg yet equations are different)."""
    return R*(R + 1)/(6.*np.sqrt(2*np.abs(R - 1)))


def L(f, N):
    f30 = 7.292115e-5  # rad s-1
    N0 = 5.2e-3  # rad s-1
    f = np.abs(f)
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
        win = sig.get_window('boxcar', len(x))
    elif window == 'sin2taper':
        win = sin2taper(len(x))
    else:
        win = sig.get_window(window, len(x))

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
    Pxy = xft*yft.conj()

    M = nfft/2 + 1

    # Chop spectrum in half.
    Pxx, Pyy, Pxy = Pxx[:M], Pyy[:M], Pxy[:M]

    # Make sure the zero frequency is really zero and not a very very small
    # non-zero number because that can mess up log plots.
    if detrend is not None:
        Pxx[..., 0], Pyy[..., 0], Pxy[..., 0] = 0., 0., 0.

    # Multiply spectrum by 2 except for the Nyquist and constant elements to
    # account for the loss of negative frequencies.
    Pxx[..., 1:-1] *= 2*scale
    Pxx[..., (0, -1)] *= scale

    Pyy[..., 1:-1] *= 2*scale
    Pyy[..., (0, -1)] *= scale

    Pxy[..., 1:-1] *= 2*scale
    Pxy[..., (0, -1)] *= scale

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

    m_int = np.logspace(np.log10(m_0), np.log10(m_c), 100)
    P_int = np.interp(m_int, m, P)

    return np.trapz(P_int, x=m_int)


def spectral_correction(m, use_range=True, use_diff=True, use_interp=True,
                        use_tilt=True, use_bin=True, use_volt=True, dzt=8.,
                        dzr=8., dzfd=8., dzg=8., ddash=5.4, dzs=8., vfi=50.,
                        mfr=0.12):
    """
    Calculates the appropriate transfer function to scale power spectra.

    Parameters
    ----------
    m : ndarray
        Vertical wavenumber. [rad s-1]
    use_range : boolean, optional (LADCP)
        Switch for range correction.
    use_diff : boolean, optional (LADCP)
        Switch for differencing correction.
    use_interp : boolean, optional (LADCP/EM-APEX)
        Switch for interpolation correction.
    use_tilt : boolean, optional (LADCP)
        Switch for tilt correction.
    use_bin : boolean, optional (LADCP)
        Switch for binning correction.
    use_volt : boolean, optional (EM-APEX)
        Switch for voltmeter noise correction.
    dzt : float, optional (LADCP)
        Transmitted sound pulse length projected on the vertical. [m]
    dzr : float, optional (LADCP/EM-APEX)
        Receiver processing bin length. [m]
    dzfd : float, optional (LADCP)
        First-differencing interval. [m]
    dzg : float, optional (LADCP/EM-APEX)
        Interval of depth grid onto which single-ping piecewise-linear
        continuous profiles of vertical shear are binned. [m]
    ddash : float, optional (LADCP)
        ?
    dzs : float, optional (LADCP)
        Superensemble pre-averaging interval, usually chosen to be dzg. [m]
    vfi : float, optional (EM-APEX)
        ? [s-1]
    mfr : float, optional (EM-APEX)
        ? [m s-1]

    Returns
    -------
    T : ndarray
        Transfer function, which is the product of all of the individual
        transfer functions for each separate spectral correction.


    Notes
    -----
    Spectral corrections for LADCP data - see Polzin et. al. 2002.

    There is another possible correction which isn't used.

    Notes from MATLAB code
    ----------------------

    A quadratic fit to the range maxima (r_max) pairs given by Polzin et al.
    (2002) yields.

    ddash = -1.2+0.0857r_max - 0.000136r_max^2 ,

    which has an intercept near r_max = 14 m. It should be noted that
    expressions (4) and (5) are semi-empirical and apply strictly only to the
    data set of Polzin et al. (2002). Estimating r_max ? 255 m as the range at
    which 80% of all ensembles have valid velocities yields d? ? 11.8 m in case
    of this data set (i.e as in Thurherr 2011 NOT DIMES - need to update!!).
    ddash is determined empirically by Polzin et al. (2002) and is dependent
    on range and the following assumptions:
        Small tilts (~ 3 deg).
        Instrument tilt and orientation are constant over measurement period.
        Instrument tilt and orientation are independent.
        Tilt attenuation is limited by bin-mapping capabilities of RDI (1996)
        processing.

    """

    pi2 = np.pi*2

    # Range averaging.
    if use_range:
        T_range = np.sinc(m*dzt/pi2)**2 * np.sinc(m*dzr/pi2)**2
    else:
        T_range = 1.

    # First differencing.
    if use_diff:
        T_diff = np.sinc(m*dzfd/pi2)**2
    else:
        T_diff = 1.

    # Interpolation.
    if use_interp:
        T_interp = np.sinc(m*dzr/pi2)**4 * np.sinc(m*dzg/pi2)**2
    else:
        T_interp = 1.

    # Tilting.
    if use_tilt:
        T_tilt = np.sinc(m*ddash/pi2)**2
    else:
        T_tilt = 1.

    # Binning
    if use_bin:
        T_bin = np.sinc(m*dzg/pi2)**2 * np.sinc(m*dzs/pi2)**2
    else:
        T_bin = 1.

    # Voltmeter
    if use_volt:
        T_volt = 1./np.sinc(m*vfi*mfr/pi2)
    else:
        T_volt = 1.

    T = T_range*T_diff*T_interp*T_tilt*T_bin*T_volt

    return T


def window_ps(dz, U, V, dUdz, dVdz, strain, N2_ref, params=default_params):
    """Calculate the power spectra for a window of data."""

    # Normalise the shear by the mean buoyancy frequency.
    ndUdz = dUdz/np.mean(np.sqrt(N2_ref))
    ndVdz = dVdz/np.mean(np.sqrt(N2_ref))

    # Compute the (co)power spectral density.
    m, PdU, PdV, PdUdV = coperiodogram(ndUdz, ndVdz, fs=1./dz,
                                       **params['periodogram_params'])
    # We only really want the cospectrum for shear so the next two lines are
    # something of a hack where we ignore unwanted output.
    __, PU, PV, __ = coperiodogram(U, V, fs=1./dz,
                                   **params['periodogram_params'])
    __, Pstrain, __, __ = coperiodogram(strain, U, fs=1./dz,
                                        **params['periodogram_params'])

    # Clockwise and counter clockwise spectra.
    PCW = CW_ps(PdU, PdV, PdUdV)
    PCCW = CCW_ps(PdU, PdV, PdUdV)
    # Shear spectra.
    Pshear = PdU + PdV

    if params['apply_corrections']:
        T = spectral_correction(m, **params['corrections'])
        PCW /= T
        PCCW /= T
        Pshear /= T
        PU /= T
        PV /= T

        if params['print_diagnostics']:
            print("T = {}".format(T))

    # Kinetic energy spectra.
    PEK = (PU + PV)/2.

    return m, Pshear, Pstrain, PCW, PCCW, PEK


def analyse(z, U, V, dUdz, dVdz, strain, N2_ref, lat, params=default_params):
    """ """

    X = [U, V, dUdz, dVdz, strain, N2_ref]

    if params['plot_profiles']:
        fig, axs = plt.subplots(1, 4, sharey=True)

        axs[0].set_ylabel('$z$ (m)')
        axs[0].plot(np.sqrt(N2_ref), z, 'k-', label='$N_{ref}$')
        axs[0].plot(np.sqrt(strain*N2_ref + N2_ref), z, 'k--', label='$N$')
        axs[0].set_xlabel('$N$ (rad s$^{-1}$)')
        axs[0].legend(loc=0)
        axs[0].set_xticklabels(axs[0].get_xticks(), rotation='vertical')
        axs[1].plot(U, z, 'k-', label='$U$')
        axs[1].plot(V, z, 'r-', label='$V$')
        axs[1].set_xlabel('$U$, $V$ (m s$^{-1}$)')
        axs[1].legend(loc=0)
        axs[1].set_xticklabels(axs[1].get_xticks(), rotation='vertical')
        axs[2].plot(dUdz, z, 'k-', label=r'$\frac{dU}{dz}$')
        axs[2].plot(dVdz, z, 'r-', label=r'$\frac{dV}{dz}$')
        axs[2].set_xlabel(r'$\frac{dU}{dz}$, $\frac{dV}{dz}$ (s$^{-1}$)')
        axs[2].legend(loc=0)
        axs[2].set_xticklabels(axs[2].get_xticks(), rotation='vertical')
        axs[3].plot(strain, z, 'k-')
        axs[3].set_xlabel(r'$\xi_z$ (-)')

    # Split varables into overlapping window segments, bare in mind the last
    # window may not be full.
    width = params['bin_width']
    overlap = params['bin_overlap']
    wdws = [wdw.window(z, x, width=width, overlap=overlap) for x in X]

    n = wdws[0].shape[0]
    z_mean = np.empty(n)
    EK = np.empty(n)
    R_pol = np.empty(n)
    R_om = np.empty(n)
    epsilon = np.empty(n)
    kappa = np.empty(n)

    for i, w in enumerate(zip(*wdws)):

        # This takes the z values from the horizontal velocity.
        wz = w[0][0]
        z_mean[i] = np.mean(wz)
        # This (poor code) removes the z values from windowed variables.
        w = [var[1] for var in w]
        N2_mean = np.mean(w[-1])
        N_mean = np.sqrt(N2_mean)

        # Get the useful power spectra.
        m, PCW, PCCW, Pshear, Pstrain, PEK = \
            window_ps(params['dz'], *w, params=params)

        # Integrate the spectra.
        I = [integrated_ps(m, P, params['m_c'], params['m_0'])
             for P in [Pshear, Pstrain, PCW, PCCW, PEK]]

        Ishear, Istrain, ICW, ICCW, IEK = I

        # Garrett-Munk shear power spectral density normalised.
        GMshear = GM79.E_she_z(2*np.pi*m, N_mean)/N_mean

        IGMshear = integrated_ps(m, GMshear, params['m_c'], params['m_0'])

        EK[i] = IEK
        R_pol[i] = ICCW/ICW
        R_om[i] = Ishear/Istrain
        epsilon[i] = GM79.epsilon_0*N2_mean/GM79.N_0**2*Ishear**2/IGMshear**2
        # Apply correcting factors

        epsilon[i] *= L(gsw.f(lat), N_mean)*h_gregg(R_om[i])

        kappa[i] = params['mixing_efficiency']*epsilon[i]/N2_mean

        if params['print_diagnostics']:
            print("Ishear = {}".format(Ishear))
            print("IGMshear = {}".format(IGMshear))
            print("lat = {}. f = {}.".format(lat, gsw.f(lat)))
            print("N_mean = {}".format(N_mean))
            print("R_om = {}".format(R_om[i]))
            print("L = {}".format(L(gsw.f(lat), N_mean)))
            print("h = {}".format(h_gregg(R_om[i])))

        # Plotting here generates a crazy number of plots.
        if params['plot_spectra']:

            GMstrain = GM79.E_str_z(2*np.pi*m, N_mean)
            GMvel = GM79.E_vel_z(2*np.pi*m, N_mean)

            fig, axs = plt.subplots(4, 1, sharex=True)

            axs[0].loglog(m, PEK, 'k-', label="$E_{KE}$")
            axs[0].loglog(m, GMvel, 'k--', label="GM $E_{KE}$")
            axs[0].set_title("height {:1.0f} m".format(z_mean[i]))
            axs[1].loglog(m, Pshear, 'k-', label="$V_z$")
            axs[1].loglog(m, GMshear, 'k--', label="GM $V_z$")
            axs[2].loglog(m, Pstrain, 'k', label=r"$\xi_z$")
            axs[2].loglog(m, GMstrain, 'k--', label=r"GM $\xi_z$")
            axs[3].loglog(m, PCW, 'r-', label="CW")
            axs[3].loglog(m, PCCW, 'k-', label="CCW")

            axs[-1].set_xlabel('$k_z$ (m$^{-1}$)')

            for ax in axs:
                ax.vlines(params['m_c'], *ax.get_ylim())
                ax.vlines(params['m_0'], *ax.get_ylim())
                ax.grid()
                ax.legend()

    if params['plot_results']:

        fig, axs = plt.subplots(1, 5, sharey=True)

        axs[0].plot(np.log10(EK), z_mean, 'k-o')
        axs[0].set_xlabel('$\log_{10}E_{KE}$ (m$^{2}$ s$^{-2}$)')
        axs[0].set_ylabel('$z$ (m)')
        axs[1].plot(np.log10(R_pol), z_mean, 'k-o')
        axs[1].set_xlabel('$\log_{10}R_{pol}$ (-)')
        axs[1].set_xlim(-1, 1)
        axs[2].plot(np.log10(R_om), z_mean, 'k-o')
        axs[2].set_xlabel('$\log_{10}R_{\omega}$ (-)')
        axs[2].set_xlim(-1, 1)
        axs[3].plot(np.log10(epsilon), z_mean, 'k-o')
        axs[3].set_xlabel('$\log_{10}\epsilon$ (W kg$^{-1}$)')
        axs[4].plot(np.log10(kappa), z_mean, 'k-o')
        axs[4].set_xlabel('$\log_{10}\kappa$ (m$^{2}$ s$^{-1}$)')

        for ax in axs:
            ax.grid()
            ax.set_xticklabels(ax.get_xticks(), rotation='vertical')

    return z_mean, EK, R_pol, R_om, epsilon, kappa


data_file = "/noc/users/jc3e13/storage/test/ctd_ladcp.mat"
data = sp.io.loadmat(data_file, squeeze_me=True, struct_as_record=True)
ctd = data['allCTD']
ladcp = data['allLADCP']

use_station = 96

idx = use_station - 1

# Basic function...

# 1. Receive profiles of strain, shear, depth, N2, f, alternative depth.
# 2. Have an option to clean those profiles of nans and regrid onto a finer
#    grid.
# 3. Split up into bins according to user specifications (width, slide up/down)
# 4. Separately calculate shear and strain spectra with options for calculating
#    energy if U and V also given and optional spectral corrections.
# 5. Integrate spectra and calculate dissipation etc. in that bin.
# 6. Return to step 4. for next bin and continue through the profile.
# 7. Include options to print graphs and diagnostics along the way.

# 1.
etaz = ctd['etaz'][idx]
N2_ref = ctd['n2_bray'][idx]
P_ctd = ctd['p_ave'][idx]

shear = ladcp['shear'][idx]
dudz = shear.real
dvdz = shear.imag
P_ladcp = ladcp['p'][idx]

lat = np.mean(ctd['lat'][1:])
f = gsw.f(lat)

# 2.
nans_ctd = np.isnan(etaz) | np.isnan(N2_ref) | np.isnan(P_ctd)
nans_ladcp = np.isnan(dudz) | np.isnan(dvdz) | np.isnan(P_ladcp)

etaz = etaz[~nans_ctd]
N2_ref = N2_ref[~nans_ctd]
P_ctd = P_ctd[~nans_ctd]

dudz = dudz[~nans_ladcp]
dvdz = dvdz[~nans_ladcp]
P_ladcp = P_ladcp[~nans_ladcp]

P_top = np.max((np.min(P_ctd), np.min(P_ladcp)))
P_bot = np.min((np.max(P_ctd), np.max(P_ladcp)))

# If the user chooses starting depth there should be a warning if exceeds
# bounds set by P_top and P_bot.

regrid = False
if regrid:
    # Choose new grid spacing.
    dP = 1.
    P = np.arange(P_top, P_bot + dP, dP)

    ETAi = np.interp(P, P_ctd, etaz)
    N2i = np.interp(P, P_ctd, N2_ref)

    Uzi = np.interp(P, P_ladcp, dudz)
    Vzi = np.interp(P, P_ladcp, dvdz)

# 3.
width = 240.
overlap = 180.
dudz_bins = wdw.window(P_ladcp, dudz, width, overlap, x_0=P_top)
dvdz_bins = wdw.window(P_ladcp, dvdz, width, overlap, x_0=P_top)
etaz_bins = wdw.window(P_ctd, etaz, width, overlap, x_0=P_top)
N2_ref_bins = wdw.window(P_ctd, N2_ref, width, overlap, x_0=P_top)

# Bare in mind that these could be different lengths -- a possible problem.
N_bins = len(dudz_bins)

P_mean = np.empty(N_bins)
R_pol = np.empty(N_bins)
R_om = np.empty(N_bins)
epsilon = np.empty(N_bins)
kappa = np.empty(N_bins)

print_diagnostics = True

for i in xrange(N_bins):

    N2_mean = np.mean(N2_ref_bins[i][1])
    N_mean = np.sqrt(N2_mean)

    P_shear = dudz_bins[i][0]
    nUz = dudz_bins[i][1]/N_mean
    nVz = dvdz_bins[i][1]/N_mean
    P_strain = netaz = etaz_bins[i][0]
    ietaz = etaz_bins[i][1]

    P_mean[i] = (P_shear[-1] + P_shear[0])/2.

    # This should always be true...
    dP_shear = P_shear[1] - P_shear[0]
    dP_strain = P_strain[1] - P_strain[0]

    # Compute the (co)power spectral density.
    m_shear, PdU, PdV, PdUdV = coperiodogram(nUz, nVz, fs=1./dP_shear,
                                             window='hanning', nfft=None,
                                             detrend='linear',
                                             scaling='density')
    # We only really want the cospectrum for shear so the next two lines are
    # something of a hack where we ignore unwanted output.

    m_strain, Pstrain, __, __ = coperiodogram(ietaz, ietaz, fs=1./dP_strain,
                                              window='hanning', nfft=None,
                                              detrend='linear',
                                              scaling='density')

    # Clockwise and counter clockwise spectra.
    PCW = CW_ps(PdU, PdV, PdUdV)
    PCCW = CCW_ps(PdU, PdV, PdUdV)
    # Shear spectra.
    Pshear = PdU + PdV

    # Also Garrett-Munk. Factor of 2 pi to convert to cyclical units
    GMstrain = 2*np.pi*GM79.E_str_z(2*np.pi*m_strain, N_mean)
    GMvel = 2*np.pi*GM79.E_vel_z(2*np.pi*m_shear, N_mean)
    GMshear = 2*np.pi*GM79.E_she_z(2*np.pi*m_shear, N_mean)/N_mean

    # Now apply the corrections.
    T_shear = spectral_correction(m_shear, use_range=True, use_diff=True,
                                  use_interp=True, use_tilt=True, use_bin=True,
                                  use_volt=False, dzt=10., dzr=10., dzfd=10.,
                                  dzg=10., ddash=5.4, dzs=5.)

    T_strain = spectral_correction(m_strain, use_range=False, use_diff=False,
                                   use_interp=True, use_tilt=False,
                                   use_bin=False, use_volt=False, dzt=10.,
                                   dzr=10., dzfd=8., dzg=10., ddash=5.4,
                                   dzs=5.)

    PCW /= T_shear
    PCCW /= T_shear
    Pshear /= T_shear

    Pstrain /= T_strain

    # Now calculate what we wanted in the first place.
    m_c = 1./90.
    m_0 = 1./180.

    I_strain = integrated_ps(m_strain, Pstrain, m_c, m_0)
    I_shear = integrated_ps(m_shear, Pshear, m_c, m_0)
    I_CW = integrated_ps(m_shear, PCW, m_c, m_0)
    I_CCW = integrated_ps(m_shear, PCCW, m_c, m_0)
    I_GMshear = integrated_ps(m_shear, GMshear, m_c, m_0)
    I_GMstrain = integrated_ps(m_strain, GMstrain, m_c, m_0)

    R_pol[i] = I_CCW/I_CW
    R_om[i] = I_shear/I_strain

    if R_om[i] < 3. or R_om[i] > 20.:
        R = 7.
    else:
        R = R_om[i]

    epsilon[i] = GM79.epsilon_0*N2_mean/GM79.N_0**2*I_shear/I_GMshear \
        * L(f, N_mean)* h_gregg(R)
    kappa[i] = 0.2*epsilon[i]/N2_mean

    # Now print diagnostics and plots to see if it makes any sense.

    if print_diagnostics:

    #    print("Ishear = {}".format(Ishear))
    #    print("IGMshear = {}".format(IGMshear))
    #    print("N_mean = {}".format(N_mean))
        print("R_om = {}".format(R_om[i]))
    #    print("L = {}".format(L(gsw.f(lat), N_mean)))
    #    print("h = {}".format(h_gregg(R_om[i])))

        fig, axs = plt.subplots(4, 1, sharex=True)

    #    axs[0].loglog(m, PEK, 'k-', label="$E_{KE}$")
        axs[0].loglog(m_shear, GMvel, 'k--', label="GM $E_{KE}$")
    #    axs[0].set_title("height {:1.0f} m".format(z_mean[i]))
        axs[1].loglog(m_shear, Pshear, 'k-', label="$V_z$")
        axs[1].loglog(m_shear, GMshear, 'k--', label="GM $V_z$")
        axs[2].loglog(m_strain, Pstrain, 'k', label=r"$\xi_z$")
        axs[2].loglog(m_strain, GMstrain, 'k--', label=r"GM $\xi_z$")
        axs[3].loglog(m_shear, PCW, 'r-', label="CW")
        axs[3].loglog(m_shear, PCCW, 'k-', label="CCW")

        axs[-1].set_xlabel('$k_z$ (m$^{-1}$)')

        for ax in axs:
            ax.vlines(m_c, *ax.get_ylim())
            ax.vlines(m_0, *ax.get_ylim())
            ax.grid()
            ax.legend()

z_mean = gsw.z_from_p(P_mean, lat)

fig, axs = plt.subplots(1, 4, sharey=True, figsize=(6.5, 3.5))

axs[0].set_ylabel('$z$ (m)')
axs[0].plot(np.log10(R_pol), z_mean, 'k-o')
axs[0].set_xlabel('$\log_{10}R_{pol}$ (-)')
axs[0].set_xlim(-1, 1)
axs[1].plot(np.log10(R_om), z_mean, 'k-o')
axs[1].set_xlabel('$\log_{10}R_{\omega}$ (-)')
axs[1].set_xlim(-1, 1)
axs[2].plot(np.log10(epsilon), z_mean, 'k-o')
axs[2].set_xlabel('$\log_{10}\epsilon$ (W kg$^{-1}$)')
axs[3].plot(np.log10(kappa), z_mean, 'k-o')
axs[3].set_xlabel('$\log_{10}\kappa$ (m$^{2}$ s$^{-1}$)')

