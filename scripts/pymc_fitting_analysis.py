# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 12:58:56 2014

@author: jc3e13
"""

import numpy as np
import sys
import matplotlib.pyplot as plt
import os

import gsw
import triangle

lib_path = os.path.abspath('../modules')
if lib_path not in sys.path:
    sys.path.append(lib_path)

import pymc
import emapex
import float_advection_routines as far
import utils
import gravity_waves as gw


try:
    print("Floats {} and {} exist!.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.load(4976)
    E77 = emapex.load(4977)


# %% Using the integrated model


pfl26 = E77.get_profiles(26)

zmax = -650.
use = ~np.isnan(pfl26.z) & (pfl26.z < zmax)
zf = pfl26.z[use]
wf = pfl26.Ww[use]
uf = pfl26.U_abs[use]
vf = pfl26.V_abs[use]
bf = pfl26.b[use]
zmin = np.min(zf)

uf = utils.nan_detrend(zf, uf, 2)
vf = utils.nan_detrend(zf, vf, 2)

data_stack = np.hstack((uf, vf, wf, 250*bf))

plt.figure()
plt.plot(data_stack)


def model():

    # Priors.
    sig = pymc.Uniform('sig', 0.0, 5., value=0.01)
    X = pymc.Uniform('X', -2e4, 2e4, value=-2e3)
    Y = pymc.Uniform('Y', -1e5, 1e5, value=-4e3)
    Z = pymc.Uniform('Z', -2e4, 2e4, value=-2e3)
    phase = pymc.Uniform('phase', 0., np.pi*2, value=0.)

    @pymc.deterministic()
    def wave_model(zf=zf, X=X, Y=Y, Z=Z, phase=phase):
        return far.model_pymc(zf, X, Y, Z, phase)

    # Likelihood
    y = pymc.Normal('y', mu=wave_model, tau=1./sig**2, value=data_stack, observed=True)

    return locals()

M = pymc.MCMC(model(), db='pickle', dbname='trace.p')
samples = 100000
burn = 50000
thin = 5
M.sample(samples, burn, thin)
pymc.Matplot.plot(M, common_scale=False)

plt.figure()
plt.plot(data_stack)
plt.plot(far.model_pymc(zf, np.median(M.trace('X')[:]), np.median(M.trace('Y')[:]),
                      np.median(M.trace('Z')[:])))
plt.savefig('output.png')

plt.figure()
plt.plot(data_stack)
plt.plot(far.model_pymc(zf, 5000., 15000, 8000, 0.))

# %% Fitting to profiles


def w_model(params, data):

    phi_0, X, Y, Z, phase_0 = params

    time, x, y, z, U, V, W, B, N, f = data

    k = 2*np.pi/X
    l = 2*np.pi/Y
    m = 2*np.pi/Z

    om = gw.omega(N, k, m, l, f)

    w = gw.w(x, y, z, time, phi_0, k, l, m, om, N, phase_0=phase_0)

    return w


def u_model(params, data):

    phi_0, X, Y, Z, phase_0 = params

    time, x, y, z, U, V, W, B, N, f = data

    k = 2*np.pi/X
    l = 2*np.pi/Y
    m = 2*np.pi/Z

    om = gw.omega(N, k, m, l, f)

    u = gw.u(x, y, z, time, phi_0, k, l, m, om, phase_0=phase_0)

    return u


def v_model(params, data):

    phi_0, X, Y, Z, phase_0 = params

    time, x, y, z, U, V, W, B, N, f = data

    k = 2*np.pi/X
    l = 2*np.pi/Y
    m = 2*np.pi/Z

    om = gw.omega(N, k, m, l, f)

    v = gw.v(x, y, z, time, phi_0, k, l, m, om, phase_0=phase_0)

    return v


def b_model(params, data):

    phi_0, X, Y, Z, phase_0 = params

    time, x, y, z, U, V, W, B, N, f = data

    k = 2*np.pi/X
    l = 2*np.pi/Y
    m = 2*np.pi/Z

    om = gw.omega(N, k, m, l, f)

    b = gw.b(x, y, z, time, phi_0, k, l, m, om, N, phase_0=phase_0)

    return 250.*b


def full_model(params, data):
    return np.hstack((u_model(params, data),
                      v_model(params, data),
                      w_model(params, data),
                      b_model(params, data)))


# Previously this looked like E76.get_timeseries([31, 32], ) etc. and the below
# bits of code were uncommented.

# %% PROFILES 31 32 ###########################################################

time, z = E76.get_timeseries([31, 32], 'z')
timeef, U = E76.get_timeseries([31, 32], 'U')
__, V = E76.get_timeseries([31, 32], 'V')
__, W = E76.get_timeseries([31, 32], 'Ww')
__, B = E76.get_timeseries([31, 32], 'b')
__, N2 = E76.get_timeseries([31, 32], 'N2_ref')
__, x = E76.get_timeseries([31, 32], 'x_ctd')
__, y = E76.get_timeseries([31, 32], 'y_ctd')

t_split = E76.get_profiles(31).UTC_end

nope = z > -600.

time = time[~nope]
W = W[~nope]
B = B[~nope]
x = x[~nope]
y = y[~nope]
z = z[~nope]

N = np.nanmean(np.sqrt(N2))
f = gsw.f(-57.5)

Unope = np.isnan(U)

timeef = timeef[~Unope]
U = U[~Unope]
V = V[~Unope]

U = np.interp(time, timeef, U)
V = np.interp(time, timeef, V)

U[time > t_split] = utils.nan_detrend(z[time > t_split], U[time > t_split], 2)
U[time < t_split] = utils.nan_detrend(z[time < t_split], U[time < t_split], 2)
U[U > 0.3] = 0.

V[time > t_split] = utils.nan_detrend(z[time > t_split], V[time > t_split], 2)
V[time < t_split] = utils.nan_detrend(z[time < t_split], V[time < t_split], 2)

#U = utils.nan_detrend(depth, U, 2)
#V = utils.nan_detrend(depth, V, 2)

time *= 60.*60.*24
time -= np.min(time)

data = [time, x, y, z, U, V, W, B, N, f]

data_stack = np.hstack((U, V, W, 250.*B))


def model():

    # Priors.
#    sig = pymc.Uniform('sig', 0.0, 5., value=0.01)
    sig = 0.02
    phi_0 = pymc.Uniform('phi_0', 0, 10, value=0.3)
    X = pymc.Uniform('X', -100000., -500., value=-1500.)
    Y = pymc.Uniform('Y', -100000., -500., value=-1500.)
    Z = pymc.Uniform('Z', -50000., -500., value=-1000.)
    phase = pymc.Uniform('phase', 0., np.pi*2, value=1.)

    @pymc.deterministic()
    def wave_model(phi_0=phi_0, X=X, Y=Y, Z=Z, phase=phase):
        params = [phi_0, X, Y, Z, phase]
        return full_model(params, data)

    # Likelihood
    y = pymc.Normal('y', mu=wave_model, tau=1./sig**2, value=data_stack,
                    observed=True)

    return locals()

M = pymc.MCMC(model(), db='pickle', dbname='trace_31_32.p') # , db='pickle', dbname='trace.p'
samples = 12000000
burn = 11800000
thin = 10
M.sample(samples, burn, thin)
pymc.Matplot.plot(M, common_scale=False)

triangle.corner(np.transpose(np.asarray([M.trace('X')[:], M.trace('Y')[:],
                                         M.trace('Z')[:], M.trace('phi_0')[:],
                                         M.trace('phase')[:]])),
                labels=['$\lambda_x$', '$\lambda_y$', '$\lambda_z$',
                        '$\phi_0$', 'phase'])


plt.figure()
plt.plot(data_stack)

Np = (samples - burn)/thin

for i in xrange(0, Np, 20):
    params = [M.trace('phi_0')[i], M.trace('X')[i], M.trace('Y')[i],
              M.trace('Z')[i], M.trace('phase')[i]]
    plt.plot(full_model(params, data), 'k', alpha=0.01)

# %% PROFILE 26 ###############################################################

time, z = E77.get_timeseries([26], 'z')
timeef, U = E77.get_timeseries([26], 'U')
__, V = E77.get_timeseries([26], 'V')
__, W = E77.get_timeseries([26], 'Ww')
__, B = E77.get_timeseries([26], 'b')
__, N2 = E77.get_timeseries([26], 'N2_ref')
__, x = E77.get_timeseries([26], 'x_ctd')
__, y = E77.get_timeseries([26], 'y_ctd')

nope = z > -600.

time = time[~nope]
W = W[~nope]
B = B[~nope]
x = x[~nope]
y = y[~nope]
z = z[~nope]

N = np.nanmean(np.sqrt(N2))
f = gsw.f(-57.5)

Unope = np.isnan(U)

timeef = timeef[~Unope]
U = U[~Unope]
V = V[~Unope]

U = np.interp(time, timeef, U)
V = np.interp(time, timeef, V)

U = utils.nan_detrend(z, U, 2)
V = utils.nan_detrend(z, V, 2)

time *= 60.*60.*24
time -= np.min(time)

data = [time, x, y, z, U, V, W, B, N, f]

data_stack = np.hstack((U, V, W, 250.*B))


def model():

    # Priors.
#    sig = pymc.Uniform('sig', 0.0, 5., value=0.01)
    sig = 0.02
    phi_0 = pymc.Uniform('phi_0', 0, 10, value=0.1)
    X = pymc.Uniform('X', -100000., 100000., value=-15000.)
    Y = pymc.Uniform('Y', -100000., 100000, value=-15000.)
    Z = pymc.Uniform('Z', -100000., 100000, value=-25000.)
    phase = pymc.Uniform('phase', -1000., 1000., value=1.)

    @pymc.deterministic()
    def wave_model(phi_0=phi_0, X=X, Y=Y, Z=Z, phase=phase):
        params = [phi_0, X, Y, Z, phase]
        return full_model(params, data)

    # Likelihood
    y = pymc.Normal('y', mu=wave_model, tau=1./sig**2, value=data_stack,
                    observed=True)

    return locals()

M = pymc.MCMC(model(), db='pickle', dbname='trace_26.p')
samples = 500000
burn = 300000
thin = 10
M.sample(samples, burn, thin)
pymc.Matplot.plot(M, common_scale=False)

triangle.corner(np.transpose(np.asarray([M.trace('X')[:], M.trace('Y')[:],
                                         M.trace('Z')[:], M.trace('phi_0')[:],
                                         M.trace('phase')[:]])),
                labels=['$\lambda_x$', '$\lambda_y$', '$\lambda_z$',
                        '$\phi_0$', 'phase'])


plt.figure()
plt.plot(data_stack)

Np = (samples - burn)/thin

for i in xrange(0, Np, 20):
    params = [M.trace('phi_0')[i], M.trace('X')[i], M.trace('Y')[i],
              M.trace('Z')[i], M.trace('phase')[i]]
    plt.plot(full_model(params, data), 'k', alpha=0.01)
