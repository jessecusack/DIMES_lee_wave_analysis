# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 16:31:11 2014

@author: jc3e13

This module contains functions for investigating internal gravity waves.

"""

import numpy as np


def plane_wave(amplitude, wavevector, position):
    """Simple equation for a plane wave.

    Parameters
    ----------
    amplitude : float
        The wave amplitude.
    wavevector : ndarray
        Shape (1, N) array, e.g. (k, l, m, omega).
    position : ndarray
        Shape (N, M) array, e.g. (x, y, z, t).

    Returns
    -------
    wave : ndarray
        Wave amplitude at given positions. Shape (1, M).

    """

    wave = amplitude*np.real(np.exp(1j*np.dot(wavevector, position)))
    return wave


def omega(N, k, m, l=None, f=None):
    """Dispersion relation for an internal gravity wave in a continuously
    stratified fluid.

    (Vallis 2012)

    Parameters
    ----------
    N : ndarray
        Buoyancy frequency [s-1]
    k : ndarray
        Horizontal wavenumber (x) [m-1]
    m : ndarray
        Vertical wavenumber (z) [m-1]
    l : ndarray, optional
        Horizontal wavenumber (y) [m-1]
    f : ndarray, optional
        Coriolis parameter [s-1]

    Returns
    -------
    omega : ndarray
        Frequency [s-1]

    Notes
    -----
    The appropriate equation will be used based on the function arguments.

    """

    N2 = N**2
    k2 = k**2
    m2 = m**2

    if l is None and f is None:
        # 2D case without rotation.
        return np.sqrt(N2*k2/(k2 + m2))
    elif l is not None and f is None:
        # 3D case without rotation.
        l2 = l**2
        return np.sqrt(N2*(k2 + l2)/(k2 + l2 + m2))
    elif l is None and f is not None:
        # 2D case with rotation.
        f2 = f**2
        return np.sqrt((f2*m2 + N2*k2)/(k2 + m2))
    elif l is not None and f is not None:
        # 3D case with rotation.
        f2 = f**2
        l2 = l**2
        return np.sqrt((f2*m2 + N2*(k2 + l2))/(k2 + l2 + m2))


def m_topo(k, N, U, f):
    """Relationship between vertical and horizontal wavenumbers for a
    topographically generated internal wave.

    (Vallis 2012)

    Parameters
    ----------
    N : ndarray
        Buoyancy frequency [s-1]
    k : ndarray
        Horizontal wavenumber [m-1]
    U : ndarray
        Bottom flow speed [m s-1]
    f : ndarray
        Coriolis parameter [s-1]

    Returns
    -------
    m : ndarray
        Vertical wavenumber [m-1]

    """

    k2 = k**2.
    N2 = N**2.
    U2 = U**2.
    f2 = f**2.

    return np.sqrt((k2*(N2 - U2*k2))/(U2*k2 - f2))


def witch_of_agnesi(x, a=1., h=1.):
    return h*a**2/(a**2 + x**2)


if __name__ == '__main__':
    pass