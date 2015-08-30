# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 11:59:10 2015

@author: jc3e13
"""

import argparse
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt
import os

import gsw
import triangle

lib_path = os.path.abspath('../modules')
if lib_path not in sys.path:
    sys.path.append(lib_path)

import pymc
import emapex
import utils
import gravity_waves as gw
import plotting_functions as pf


# Figure save path.
sdir = '../figures/pymc_fitting_analysis'
if not os.path.exists(sdir):
    os.makedirs(sdir)

matplotlib.rc('font', **{'size': 8})

# Parse arguments...
parser = argparse.ArgumentParser(description='Run MCMC fitting on a profile.')
parser.add_argument('--floatID', type=int, help='EM-APEX float ID number')
parser.add_argument('--hpid', type=int, help='half profile ID number')
parser.add_argument('--zrange', nargs=2, type=float,
                    help='min max height range for fit')
parser.add_argument('--xyz', nargs=3, type=float, help='initial conditions for'
                    ' fit X, Y, Z')
args = parser.parse_args()

# Model
wscale = 1.5
bscale = 250.


def w_model(params, data):

    phi_0, X, Y, Z, phase_0 = params

    time, x, y, z, U, V, N, f = data

    k = 2*np.pi/X
    l = 2*np.pi/Y
    m = 2*np.pi/Z

    om = gw.omega(N, k, m, l, f)

    w = gw.w(x, y, z, time, phi_0, k, l, m, om, N, U=U, V=V, phase_0=phase_0)

    return wscale*w


def u_model(params, data):

    phi_0, X, Y, Z, phase_0 = params

    time, x, y, z, U, V, N, f = data

    k = 2*np.pi/X
    l = 2*np.pi/Y
    m = 2*np.pi/Z

    om = gw.omega(N, k, m, l, f)

    u = gw.u(x, y, z, time, phi_0, k, l, m, om, f=f, U=U, V=V, phase_0=phase_0)

    return u


def v_model(params, data):

    phi_0, X, Y, Z, phase_0 = params

    time, x, y, z, U, V, N, f = data

    k = 2*np.pi/X
    l = 2*np.pi/Y
    m = 2*np.pi/Z

    om = gw.omega(N, k, m, l, f)

    v = gw.v(x, y, z, time, phi_0, k, l, m, om, f=f, U=U, V=V, phase_0=phase_0)

    return v


def b_model(params, data):

    phi_0, X, Y, Z, phase_0 = params

    time, x, y, z, U, V, N, f = data

    k = 2*np.pi/X
    l = 2*np.pi/Y
    m = 2*np.pi/Z

    om = gw.omega(N, k, m, l, f)

    b = gw.b(x, y, z, time, phi_0, k, l, m, om, N, U=U, V=V, phase_0=phase_0)

    return bscale*b


def full_model(params, data):
    return np.hstack((u_model(params, data),
                      v_model(params, data),
                      w_model(params, data),
                      b_model(params, data)))


Float = emapex.load(args.floatID)

time, z = Float.get_timeseries([args.hpid], 'z')
timeef, U = Float.get_timeseries([args.hpid], 'U_abs')
__, V = Float.get_timeseries([args.hpid], 'V_abs')
__, W = Float.get_timeseries([args.hpid], 'Ww')
__, B = Float.get_timeseries([args.hpid], 'b')
__, N2 = Float.get_timeseries([args.hpid], 'N2_ref')
__, x = Float.get_timeseries([args.hpid], 'x_ctd')
__, y = Float.get_timeseries([args.hpid], 'y_ctd')

zmin, zmax = args.zrange
nope = (z > zmax) | (z < zmin)

time = time[~nope]
W = W[~nope]
B = B[~nope]
x = x[~nope]
y = y[~nope]
z = z[~nope]

N = np.nanmean(np.sqrt(N2))
f = gsw.f(Float.get_profiles(args.hpid).lat_start)

Unope = np.isnan(U)

timeef = timeef[~Unope]
U = U[~Unope]
V = V[~Unope]

U = np.interp(time, timeef, U)
V = np.interp(time, timeef, V)

Umean = np.mean(U)
Vmean = np.mean(V)

U = utils.nan_detrend(z, U, 2)
V = utils.nan_detrend(z, V, 2)

time *= 60.*60.*24
time -= np.min(time)

data = [time, x, y, z, Umean, Vmean, N, f]

data_stack = np.hstack((U, V, wscale*W, bscale*B))

X0, Y0, Z0 = args.xyz


def model():

    # Priors.
    sig = 0.02
    phi_0 = pymc.Uniform('phi_0', 0, 100, value=0.1)
    X = pymc.Uniform('X', -1000000., 1000000., value=X0)
    Y = pymc.Uniform('Y', -100000., 1000000., value=Y0)
    Z = pymc.Uniform('Z', -50000., 50000., value=Z0)
    phase = pymc.Uniform('phase', 0., np.pi*4, value=2.)

    @pymc.deterministic()
    def wave_model(phi_0=phi_0, X=X, Y=Y, Z=Z, phase=phase):
        params = [phi_0, X, Y, Z, phase]
        return full_model(params, data)

    # Likelihood
    y = pymc.Normal('y', mu=wave_model, tau=1./sig**2, value=data_stack,
                    observed=True)

    return locals()

save_string = str(args.floatID) + '_' + str(args.hpid) + '_X' + str(int(X0)) \
    + '_Y' + str(int(Y0)) + '_Z'+str(int(Z0))

dbname = '/noc/users/jc3e13/storage/processed/trace_' + save_string + '.p'
tfname = '/noc/users/jc3e13/storage/processed/results_' + save_string + '.txt'

M = pymc.MCMC(model(), db='pickle', dbname=dbname)
samples = 10000000
burn = 9800000
thin = 10
M.sample(samples, burn, thin)

# Analysis and plotting.
phi_0 = M.trace('phi_0')[:]
k = np.pi*2/M.trace('X')[:]
l = np.pi*2/M.trace('Y')[:]
m = np.pi*2/M.trace('Z')[:]

om = gw.omega(N, k, m, l, f)
w_0 = gw.W_0(phi_0, m, om, N)
Efluxz = gw.Efluxz(w_0, k, m, N, l, f)
Mfluxz = gw.Mfluxz(phi_0, k, l, m, om, N)

with open(tfname, 'w') as f:
    f.writelines("Mean frequency: {} +/- {}\n"
                 "Mean X: {} +/- {}\n"
                 "Mean Y: {} +/- {}\n"
                 "Mean Z: {} +/- {}\n"
                 "Mean phi_0: {} +/- {}\n"
                 "Mean vertical energy flux: {} +/- {}\n"
                 "Mean vertical momentum flux: {} +/- {}\n"
                 "".format(
                     np.mean(om), np.std(om),
                     np.mean(M.trace('X')[:]), np.std(M.trace('X')[:]),
                     np.mean(M.trace('Y')[:]), np.std(M.trace('Y')[:]),
                     np.mean(M.trace('Z')[:]), np.std(M.trace('Z')[:]),
                     np.mean(M.trace('phi_0')[:]), np.std(M.trace('phi_0')[:]),
                     np.mean(Efluxz), np.std(Efluxz),
                     np.mean(Mfluxz), np.std(Mfluxz)
                 ))

triangle.corner(np.transpose(np.asarray([M.trace('X')[:],
                                         M.trace('Y')[:],
                                         M.trace('Z')[:],
                                         M.trace('phi_0')[:]])),
                labels=['$\lambda_x$ (m)', '$\lambda_y$ (m)',
                        '$\lambda_z$ (m)', '$\phi_0$ (m$^2$ s$^{-2}$)'])
fig = plt.gcf()
pf.my_savefig(fig, save_string, 'MCMC_triangle', sdir, ftype='png',
              fsize='double_col')
plt.close()

# Plot fit comparison.
fig, axs = plt.subplots(1, 4, sharey=True, figsize=(6.5, 3))
axs[0].plot(100.*U, z, color='black')
axs[0].set_xlabel('$u$ (cm s$^{-1}$)')
axs[0].set_ylabel('$z$ (m)')
axs[1].plot(100.*V, z, color='black')
axs[1].set_xlabel('$v$ (cm s$^{-1}$)')
axs[2].plot(100.*W, z, color='black')
axs[2].set_xlabel('$w$ (cm s$^{-1}$)')
axs[3].plot(10000.*B, z, color='black')
axs[3].set_xlabel('$b$ ($10^{-4}$ m s$^{-2}$)')

Ns = (samples - burn)/thin

for i in xrange(0, Ns, 40):
    params = [M.trace('phi_0')[i], M.trace('X')[i], M.trace('Y')[i],
              M.trace('Z')[i], M.trace('phase')[i]]
    axs[0].plot(100.*u_model(params, data), z, color='red', alpha=0.03)
    axs[1].plot(100.*v_model(params, data), z, color='red', alpha=0.03)
    axs[2].plot(100.*w_model(params, data)/wscale, z, color='red', alpha=0.03)
    axs[3].plot(10000.*b_model(params, data)/bscale, z, color='red',
                alpha=0.03)

for ax in axs:
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=60)

pf.my_savefig(fig, save_string, 'MCMC_profiles', sdir, ftype='png',
              fsize='double_col')
plt.close()
