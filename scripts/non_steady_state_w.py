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

Float = E76

pfl = Float.get_profiles(32)

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

wsi = np.interp(z, X[:, 0], X[:, 1])

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