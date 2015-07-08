# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 11:55:19 2015

@author: jc3e13
"""

import numpy as np
import scipy as sp
import sys
import os
import matplotlib.pyplot as plt

lib_path = os.path.abspath('../modules')
if lib_path not in sys.path:
    sys.path.append(lib_path)

import emapex
import gravity_waves as gw


try:
    print("Floats {} and {} exist!.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.load(4976)
    E77 = emapex.load(4977)


def dXdt(X, t, z, rho, p, k, V_0, M, alpha_p, alpha_k, p_0, k_0, CA, g):
    """The fully nonlinear equation of motion for a float."""

    zs = X[0]
    ws = X[1]

    # Should these be reference profiles rather than measured? Or use time as
    # the interpolant instead.
    rhoi = np.interp(zs, z, rho)
    pi = np.interp(zs, z, p)
    ki = np.interp(zs, z, k)

    F_grav = g
    F_buoy = -g*rhoi*V_0/M
    F_comp = g*rhoi*V_0/M*alpha_p*(pi - p_0)
    F_pist = -g*rhoi/M*alpha_k*(ki - k_0)
    F_drag = -rhoi/M*CA*ws**2

    dwsdt = F_grav + F_buoy + F_comp + F_pist + F_drag

    dzsdt = ws

    return np.array([dzsdt, dwsdt])

# %%

Float = E76

pfl = Float.get_profiles(186)

g = -9.8
V_0 = Float.__wfi.pmean[0]
CA = Float.__wfi.pmean[1]
alpha_p = Float.__wfi.pmean[2]
p_0 = Float.__wfi.pmean[3]
alpha_k = Float.__wfi.pmean[4]
k_0 = Float.__wfi.pmean[5]
M = Float.__wfi.pmean[6]

z = pfl.z[~np.isnan(pfl.z)]
rho = pfl.interp(z, 'z', 'rho')
p = pfl.interp(z, 'z', 'P')
k = pfl.interp(z, 'z', 'ppos')
Ws = pfl.interp(z, 'z', 'Ws')
Wz = pfl.interp(z, 'z', 'Wz')
UTC = pfl.UTC[~np.isnan(pfl.z)]*86400.
UTC -= UTC[0]

w_0 = Ws[0]
z_0 = z[0]
X_0 = np.array([z_0, w_0])

t_max = 11000.
dt = 1.
t = np.arange(0., t_max+dt, dt)

args = (z, rho, p, k, V_0, M, alpha_p, alpha_k, p_0, k_0, CA, g)

X = sp.integrate.odeint(dXdt, X_0, t, args)

plt.figure()
plt.plot(X[:, 1], X[:, 0])
plt.plot(Ws, z)

wsi = np.interp(UTC, t, X[:, 1])

plt.figure()
plt.plot(Ws - wsi, z)

plt.figure()
plt.plot(Wz - wsi, z)
plt.plot(Wz - Ws, z)

Wsi = np.interp(t, UTC, Ws)
tcor = np.arange(-t_max, t_max+1, dt)
xcor = np.correlate(Wsi, X[:, 1], 'full')
plt.figure()
plt.plot(tcor, xcor)
plt.grid()

# %% Impulse test.

z = np.arange(-1500., -1492.)
rho = 1030.*np.ones_like(z)
p = z.copy()
k = -15.*np.ones_like(z)
k[k.size/2:] = 5.

w_0 = 0.
z_0 = z[0]
X_0 = np.array([z_0, w_0])

t_max = 50.
dt = 1.
t = np.arange(0., t_max+dt, dt)

args = (z, rho, p, k, V_0, M, alpha_p, alpha_k, p_0, k_0, CA, g)

X = sp.integrate.odeint(dXdt, X_0, t, args)

w_stdy = Float.__wfi.model_func((V_0, CA, alpha_p, p_0, alpha_k, k_0, M),
                                (k, p, rho), 7*[None])

plt.figure()
plt.plot(X[:, 1], X[:, 0])
plt.plot(w_stdy, z)
plt.xlabel('$w_s$ (m s$^{-1}$)')
plt.ylabel('$z$ (m)')
plt.grid()
plt.figure()
plt.plot(t, X[:, 1])
w_stdyi = np.interp(X[:, 0], z, w_stdy)
plt.plot(t, w_stdyi)
plt.xlabel('$t$ (s)')
plt.ylabel('$w_s$ (m s$^{-1}$)')
plt.grid()

# %%


def dX2dt(X, t, z, rho, p, ppos, wave, V_0, M, alpha_p, alpha_k, p_0, k_0, CA, g):
    """The fully nonlinear equation of motion for a float."""

    zs = X[0]
    ws = X[1]

    k, l, m, om, phi_0, N = wave

    # Should these be reference profiles rather than measured? Or use time as
    # the interpolant instead.
    rhoi = np.interp(zs, z, rho)
    pi = np.interp(zs, z, p)
    pposi = np.interp(zs, z, ppos)
    wi = gw.w(0., 0., zs, t, phi_0, k, l, m, om, N)

    V = V_0*(1 + alpha_p*(pi + p_0)) + alpha_k*(pposi - k_0)

    F_buoy = g - g*rhoi*V/M
    F_drag = -(rhoi/M)*CA*np.abs(ws - wi)*(ws - wi)

    dwsdt = F_buoy + F_drag

    dzsdt = ws

    return np.array([dzsdt, dwsdt])


z = np.arange(-1500., -600.)
rho = 1030.*np.ones_like(z)
p = -z.copy()
ppos = 10.*np.ones_like(z)
#ppos[ppos.size/2:] = 50.

# Wave params

X = -2000.
Y = -2000.
Z = -2000.

N = 2e-3
phi_0 = 0.03
k = np.pi*2./X
l = np.pi*2./Y
m = np.pi*2./Z
om = gw.omega(N, k, m, l)

wave = (k, l, m, om, phi_0, N)

w_0 = 0.
z_0 = z[0]
X_0 = np.array([z_0, w_0])

t_max = 5000.
dt = 1.
t = np.arange(0., t_max+dt, dt)

M = 27.179
p_0 = 2000.
k_0 = 16.
V_0 = 0.0262
CA = 0.0362
g = -9.8
alpha_p = 3.50e-6
alpha_k = 1.58e-6

args = (z, rho, p, ppos, wave, V_0, M, alpha_p, alpha_k, p_0, k_0, CA, g)

X = sp.integrate.odeint(dX2dt, X_0, t, args)

w_stdy = -Float.__wfi.model_func((V_0, CA, alpha_p, p_0, alpha_k, k_0, M),
                                 (ppos, p, rho), 7*[None])

w_stdyi = np.interp(X[:, 0], z, w_stdy)

plt.figure()
plt.plot(X[:, 1], X[:, 0])
plt.plot(w_stdy, z)
plt.plot(X[:, 1] - w_stdyi, X[:, 0])
plt.plot(gw.w(0., 0., X[:, 0], t, phi_0, k, l, m, om, N), X[:, 0])
plt.xlabel('$w$ (m s$^{-1}$)')
plt.ylabel('$z$ (m)')
plt.grid()