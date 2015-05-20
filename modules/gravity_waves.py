# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 16:31:11 2014

@author: jc3e13

This module contains functions for investigating internal gravity waves.

"""

import numpy as np


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


def phi(x, y, z, t, phi_0, k, l, m, om, U=0., phase_0=0.):
    """Pressure pertubation."""
    phase = 1j*(k*x + l*y + m*z - (om + k*U)*t + phase_0)
    return np.real(phi_0*np.exp(phase))


def b(x, y, z, t, phi_0, k, l, m, om, N, U=0., phase_0=0.):
    """Buoyancy pertubation. TODO... check weirdness happening."""
    amplitude = phi_0*1j*m*N**2/(N**2 - om**2)
    phase = 1j*(k*x + l*y + m*z - (om + k*U)*t + phase_0)
    return np.real(amplitude*np.exp(phase))


def u(x, y, z, t, phi_0, k, l, m, om, f=None, U=0., phase_0=0.):
    """Zonal velocity pertubation."""
    f = 0. if f is None else f
    amplitude = phi_0*(k*om + 1j*l*f)/(om**2 - f**2)
    phase = 1j*(k*x + l*y + m*z - (om + k*U)*t + phase_0)
    return np.real(amplitude*np.exp(phase))


def v(x, y, z, t, phi_0, k, l, m, om, f=None, U=0., phase_0=0.):
    """Meridional velocity pertubation."""
    f = 0. if f is None else f
    amplitude = phi_0*(l*om - 1j*k*f)/(om**2 - f**2)
    phase = 1j*(k*x + l*y + m*z - (om + k*U)*t + phase_0)
    return np.real(amplitude*np.exp(phase))


def w(x, y, z, t, phi_0, k, l, m, om, N, U=0., phase_0=0.):
    """Vertical velocity pertubation."""
    amplitude = -phi_0*m*om/(N**2 - om**2)
    phase = 1j*(k*x + l*y + m*z - (om + k*U)*t + phase_0)
    return np.real(amplitude*np.exp(phase))


def eta(x, y, z, t, phi_0, k, l, m, om, N, U=0., phase_0=0.):
    """Vertical displacement of isopycnals. TODO... check weirdness happening."""
    amplitude = phi_0*1j*m/(N**2 - om**2)
    phase = 1j*(k*x + l*y + m*z - (om + k*U)*t + phase_0)
    return np.real(amplitude*np.exp(phase))


def wave_vel(r, t, phi_0, U, k, l, m, om, N, f, phase_0=0.):
    x = r[..., 0]
    y = r[..., 1]
    z = r[..., 2]

    om2 = om**2
    f2 = f**2
    K2 = k**2 + l**2 + m**2

    u_0 = ((k*om + 1j*l*f)/(om2 - f2))*phi_0
    v_0 = ((l*om - 1j*k*f)/(om2 - f2))*phi_0
    w_0 = ((-om*K2)/((N**2 - f2)*m))*phi_0

    phase = 1j*(k*x + l*y + m*z - (om + k*U)*t + phase_0)

    u = np.real(u_0*np.exp(phase))
    v = np.real(v_0*np.exp(phase))
    w = np.real(w_0*np.exp(phase))

    return (np.vstack((u, v, w))).T


def buoy(r, t, phi_0, U, k, l, m, om, N, f, phase_0=0.):
    x = r[..., 0]
    y = r[..., 1]
    z = r[..., 2]

    om2 = om**2
    N2 = N**2

    b_0 = (1j*m*N2/(N2 - om2))*phi_0

    phase = 1j*(k*x + l*y + m*z - (om + k*U)*t + phase_0)

    b = np.real(b_0*np.exp(phase))

    return b


def U_0(phi_0, k, l, om, N, f):
    return np.abs(((k*om + 1j*l*f)/(om**2 - f**2))*phi_0)


def V_0(phi_0, k, l, om, N, f):
    return np.abs(((l*om - 1j*k*f)/(om**2 - f**2))*phi_0)


def W_0(phi_0, m, om, N):
    return np.abs((-om*m/(N**2 - om**2))*phi_0)


def B_0(phi_0, m, om, N):
    return np.abs((1j*m*N**2/(N**2 - om**2))*phi_0)


def ETA_0(phi_0, m, om, N):
    return np.abs(phi_0*1j*m/(N**2 - om**2))


def U(x, y, z, t, phi_0, k, m, N, l=None, om=None, f=None):
    """All components of velocity."""
    om = omega(N, k, m, l, f) if om is None else om

    if l is None:
        l = 0.
        U_x = u(x, y, z, t, phi_0, k, l, m, om, f)
        U_z = w(x, y, z, t, phi_0, k, l, m, om, N)
        return U_x, U_z
    else:
        U_x = u(x, y, z, t, phi_0, k, l, m, om, f)
        U_z = w(x, y, z, t, phi_0, k, l, m, om, N)
        U_y = v(x, y, z, t, phi_0, k, l, m, om, f)
        return U_x, U_y, U_z


def cgz(k, m, N, l=0., f=0.):
    num = -m*(N**2 - f**2)*(k**2 + l**2)
    den = (k**2 + l**2 + m**2)**1.5 * (f**2*m**2 + N**2*(k**2 + l**2))**0.5
    return num/den


def phip(k, l, m):
    return np.arctan2(m, np.sqrt(k**2 + l**2))


def Edens(w_0, k, l, m, rho_0=1025.):
    phi = phip(k, l, m)
    return 0.5*rho_0*(w_0/np.cos(phi))**2


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