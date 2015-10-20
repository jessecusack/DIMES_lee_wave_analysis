# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 16:20:33 2014

@author: jc3e13
"""

import datetime
import os
import sys
import glob
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from mpl_toolkits.axes_grid1 import make_axes_locatable
# from matplotlib.colors import LogNorm
# from matplotlib.ticker import LogFormatterMathtext
import scipy.optimize as op
from scipy.integrate import cumtrapz

import gsw

lib_path = os.path.abspath('../modules')
if lib_path not in sys.path:
    sys.path.append(lib_path)

lib_path = os.path.abspath('../../ocean-tools')
if lib_path not in sys.path:
    sys.path.append(lib_path)
    
import sandwell
import utils
import emapex
import plotting_functions as pf
import float_advection_routines as far
import gravity_waves as gw

try:
    print("Floats {} and {} exist!.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.load(4976)
    E77 = emapex.load(4977)
    E76.calculate_pressure_perturbation()
    E76.update_profiles()
    E77.calculate_pressure_perturbation()
    E77.update_profiles()

# %% Script params.

# Bathymetry file path.
bf = os.path.abspath(glob.glob('/noc/users/jc3e13/storage/smith_sandwell/topo_*.img')[0])
# Figure save path.
sdir = '../figures/sinusoid_fitting_analysis'
if not os.path.exists(sdir):
    os.makedirs(sdir)
# Universal figure font size.
matplotlib.rc('font', **{'size': 9})


# %% Is the wave stationary?

bwr = plt.get_cmap('bwr')
E76_hpids = np.arange(27, 38) # np.arange(31, 33)
E77_hpids = np.arange(22, 33) # np.arange(26, 28)

min_depth = -1400
max_depth = -200
bin_size = 40
step = 500

depth_mins = np.arange(min_depth - bin_size/2, max_depth - bin_size/2, step)
depth_maxs = np.arange(min_depth + bin_size/2, max_depth + bin_size/2, step)

fig, axs = plt.subplots(1, len(depth_mins), figsize=(16, 8), sharey=True)

for Float, hpids in zip([E76, E77], [E76_hpids, E77_hpids]):

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)

    Ww = getattr(Float, 'Ww')[:, idxs].flatten(order='F')
    z = getattr(Float, 'z')[:, idxs].flatten(order='F')
    d = getattr(Float, 'dist_ctd')[:, idxs].flatten(order='F')
    t = getattr(Float, 'UTC')[:, idxs].flatten(order='F')

    # This block gets the bathymetry at the estimated position of each CTD measurement.
    # NaNs removed so array size is different to the above block.
    tgps = getattr(Float, 'UTC_start')[idxs]
    lon = getattr(Float, 'lon_start')[idxs]
    lat = getattr(Float, 'lat_start')[idxs]
    tctd = getattr(Float, 'UTC')[:, idxs].flatten(order='F')
    nans = np.isnan(d) | np.isnan(tctd)
    tctd = tctd[~nans]
    dctd = d[~nans]
    lonctd = np.interp(tctd, tgps, lon)
    latctd = np.interp(tctd, tgps, lat)
    bathy = sandwell.interp_track(lonctd, latctd, bf)

    # Zero the distances at the top of the sea mount.
    d -= dctd[bathy.argmax()]

    for ax, dmin, dmax in zip(axs, depth_mins, depth_maxs):

        in_range = (z > dmin) & (z < dmax)

        C = ax.scatter(d[in_range], utils.datenum_to_datetime(t[in_range]),
                       c=Ww[in_range], s=30, cmap=bwr, vmin=-0.08, vmax=0.08)
        ax.set_ylim(datetime.datetime(2011, 1, 2, 12),
                    datetime.datetime(2011, 1, 4, 6))
        ax.set_xlim(-20, 40)
        ax.set_xlabel('Distance (km)')
        ax.set_title("{} to {} m".format(dmin, dmax))

#divider2 = make_axes_locatable(axs[2])
#cax2 = divider2.append_axes("right", size="10%", pad=0.4)
#cbar = plt.colorbar(C, cax=cax2, extend='both')

fmt = mdates.DateFormatter('%j %Hh')
axs[0].yaxis.set_major_formatter(fmt)
axs[0].set_ylabel('Time')
cbar = plt.colorbar(C, extend='both')
cbar.set_label('$w$ (m s$^{-1}$)')
[ax.grid() for ax in axs];

#pf.my_savefig(fig, 'both', 'time-dist', sdir, fsize='double_col')

# %% Modelling float motion

params = far.default_params

def U_const(z):
    return 0.5

params['Ufunc'] = U_const

lx = -2050.
ly = -1700.
lz = -1880.
phi_0 = 0.038
phase_0 = 3.24

pfl = E77.get_profiles(26)

params['z_0'] = np.nanmin(pfl.z)
params['N'] = 2.0e-3

X = far.model_verbose(phi_0, lx, ly, lz, phase_0, params)
X.u[:, 0] -= X.U

fig, axs = plt.subplots(1, 6, sharey=True, figsize=(16,6))
axs[0].plot(1e2*utils.nan_detrend(pfl.zef, pfl.U_abs, 2), pfl.zef, 'red')
plt.setp(axs[0].xaxis.get_majorticklabels(), rotation=60)
axs[1].plot(1e2*utils.nan_detrend(pfl.zef, pfl.V_abs, 2), pfl.zef, 'red')
plt.setp(axs[1].xaxis.get_majorticklabels(), rotation=60)
axs[2].plot(1e2*pfl.Ww, pfl.z, color='red')
plt.setp(axs[2].xaxis.get_majorticklabels(), rotation=60)
axs[3].plot(1e4*pfl.b, pfl.z, 'red')
plt.setp(axs[3].xaxis.get_majorticklabels(), rotation=60)
axs[4].plot(pfl.Pprime, pfl.z, 'red')
plt.setp(axs[4].xaxis.get_majorticklabels(), rotation=60)
axs[5].plot((pfl.dist_ctd - np.nanmin(pfl.dist_ctd))*1000., pfl.z, 'red')
plt.setp(axs[5].xaxis.get_majorticklabels(), rotation=60)

#pf.my_savefig(fig, '4977', 'pfl26_UVWB', sdir, fsize='double_col')
use = X.r[:, 2] < -600.

dwdt = 0. # utils.finite_diff(X.t, X.u[:, 2])
bw = X.b + dwdt
bi = cumtrapz(bw, X.r[:, 2], initial=0.)
bii = cumtrapz(bi, X.r[:, 2], initial=0.)
phii = bi + (bii[0] - bii[-1])/(-X.r[0, 2])
phi = gw.phi(X.r[:, 0], X.r[:, 1], X.r[:, 2], X.t, X.phi_0, X.k, X.l, X.m,
             X.om, U=U_const(0.), phase_0=phase_0)

axs[0].set_ylabel('$z$')
axs[0].plot(1e2*X.u[use, 0], X.r[use, 2])
axs[0].set_xlabel('$u$ (cm s$^{-1}$)')
plt.setp(axs[0].xaxis.get_majorticklabels(), rotation=60)
axs[1].plot(1e2*X.u[use, 1], X.r[use, 2])
axs[1].set_xlabel('$v$ (cm s$^{-1}$)')
plt.setp(axs[1].xaxis.get_majorticklabels(), rotation=60)
axs[2].plot(1e2*X.u[use, 2], X.r[use, 2])
axs[2].set_xlabel('$w$ (cm s$^{-1}$)')
plt.setp(axs[2].xaxis.get_majorticklabels(), rotation=60)
axs[3].plot(1e4*X.b[use], X.r[use, 2])
axs[3].set_xlabel('$b$ ($10^{-4}$ m s$^{-2}$)')
plt.setp(axs[3].xaxis.get_majorticklabels(), rotation=60)
axs[4].plot(phi[use], X.r[use, 2])
#axs[4].plot(phii[use], X.r[use, 2])
axs[4].set_xlabel('$\phi$ (m$^2$ s$^{-2}$)')
plt.setp(axs[4].xaxis.get_majorticklabels(), rotation=60)
axs[5].plot(X.r[use,0], X.r[use, 2])
axs[5].set_xlabel('$x$ (m)')
plt.setp(axs[5].xaxis.get_majorticklabels(), rotation=60)
plt.ylim(X.z_0, 0)

#pf.my_savefig(fig, 'model_data_pfl26', 'comparison', sdir, fsize='double_col')

k = X.k
l = X.l
m = X.m
om = X.om
N = X.N
f = X.f
U = X.U
w_0 = X.w_0
phi_0 = X.phi_0
r = X.r
t = X.t

phi = gw.phi(r[:, 0], r[:, 1], r[:, 2], t, phi_0, k, l, m, om, U, phase_0)
eta = gw.buoy(r, t, phi_0, U, k, l, m, om, N, f, phase_0)/N**2

#sintheta2 = X.m**2/(X.m**2 + X.k**2 + X.l**2)
#theta = np.rad2deg(np.arcsin(np.sqrt(sintheta2)))
#
#rho_0 = 1025.
h0 = 750.
#Eflux = m*U*w_0**2/(2*k)
#Eflux2 = 0.5*rho_0*U*m*h0**2*(U**2*k**2 - f**2)/k
#
#cp = X.om/np.sqrt(k**2 + l**2 + m**2)
#
## Group velocity.
##om0 = om - U_depth*k
#om0 = om
#cg = np.array([k*(N**2-om0**2)**2/(om0*m**2*(N**2-f**2)),
#               l*(N**2-om0**2)**2/(om0*m**2*(N**2-f**2)),
#               -(om0**2-f**2)*(N**2-om0**2)/(om0*m**2*(N**2-f**2))])
cgz = gw.cgz(k, m, N, l, f)
phip = np.arctan2(m, np.sqrt(k**2 + l**2))
#sin2phi = m**2/(k**2+l**2+m**2)
#phip = np.arcsin(np.sqrt(sin2phi))
rho0 = 1025.

E = 0.5*rho0*(w_0/np.cos(phip))**2

#Fv = 0.5*rho0*U*m/k*h0**2*(U**2*k**2 - f**2)

##wphi = phi_0**2 *

Fz = E*cgz
Fv = Fz*k/om
print("cgz", cgz)
print("E", E)
print("Fz", Fz)
print("Fv", Fv)

# %% EXTRAS
# Domain

xg, zg = np.meshgrid(np.arange(-1500, np.max(r[:,0]), 50), np.arange(-1600, 25, 25))

om2 = om**2
f2 = f**2
K2 = k**2 + l**2 + m**2
N2 = N**2

for j, ts in enumerate(np.arange(0, X.t.max(), 500.)): #

    idx = t.searchsorted(ts)
    C = []

    u_xg = np.real(((k*om + 1j*l*f)/(om2 - f2))*phi_0*np.exp(1j*(k*xg + m*zg - (om + k*U)*ts + phase_0)))
    u_yg = np.real(((l*om - 1j*k*f)/(om2 - f2))*phi_0*np.exp(1j*(k*xg + m*zg - (om + k*U)*ts + phase_0)))
    u_zg = np.real(((-om*K2)/((N**2 - f2)*m))*phi_0*np.exp(1j*(k*xg + m*zg - (om + k*U)*ts + phase_0)))
    bg = np.real((1j*m*N2/(N2 - om2))*phi_0*np.exp(1j*(k*xg + m*zg - (om + k*U)*ts + phase_0)))

    fig, axs = plt.subplots(1, 4, sharey=True, figsize=(14,6))
    fig.suptitle('t = {:1.0f} s'.format(ts))
    axs[0].set_ylabel('$z$ (m)')
    C.append(axs[0].pcolormesh(xg, zg, 1e2*(u_xg + U), cmap=plt.get_cmap('bwr')))
    axs[0].set_title('$U$ (cm s$^{-1}$)')

    divider0 = make_axes_locatable(axs[0])
    cax0 = divider0.append_axes("right", size="10%", pad=0.05)
    plt.colorbar(C[0], cax=cax0)

    C.append(axs[1].pcolormesh(xg, zg, 1e2*u_yg, cmap=plt.get_cmap('bwr')))
    axs[1].set_title('$V$ (cm s$^{-1}$)')

    divider1 = make_axes_locatable(axs[1])
    cax1 = divider1.append_axes("right", size="10%", pad=0.05)
    plt.colorbar(C[1], cax=cax1)

    C.append(axs[2].pcolormesh(xg, zg, 1e2*u_zg, cmap=plt.get_cmap('bwr')))
    axs[2].set_title('$W$ (cm s$^{-1}$)')

    divider2 = make_axes_locatable(axs[2])
    cax2 = divider2.append_axes("right", size="10%", pad=0.05)
    plt.colorbar(C[2], cax=cax2)

    C.append(axs[3].pcolormesh(xg, zg, 1e4*bg, cmap=plt.get_cmap('bwr')))
    axs[3].set_title('$b$ ($10^{-4}$ m s$^{-2}$)')

    divider3 = make_axes_locatable(axs[3])
    cax3 = divider3.append_axes("right", size="10%", pad=0.05)
    plt.colorbar(C[3], cax=cax3)

    for i in xrange(4):
        axs[i].set_xlabel('$x$ (m)')
        axs[i].set_ylim(X.z_0, 0)
        axs[i].set_xlim(np.min(xg), np.max(xg))

        axs[i].plot(r[:idx, 0], r[:idx, 2], 'k--', linewidth=3)
        axs[i].plot(r[idx, 0], r[idx, 2], 'yo', linewidth=3)
        plt.setp(axs[i].xaxis.get_majorticklabels(), rotation=60)


    C[0].set_clim(-1e2*(X.u_0+U), 1e2*(X.u_0+U))
    C[1].set_clim(-1e2*X.v_0, 1e2*X.v_0)
    C[2].set_clim(-1e2*X.w_0, 1e2*X.w_0)
    C[3].set_clim(-1e4*X.b_0, 1e4*X.b_0)

    pf.my_savefig(fig, 'model_contour', 't{:1.0f}'.format(j), sdir,
                  fsize='double_col')
    plt.close()


# %%
# Model fitting using EM-APEX positions

def u_model(params, pfl, zlim, deg):

    phi_0, X, Z, phase_0 = params
    zmin, zmax = zlim
    nope = np.isnan(pfl.zef) | (pfl.zef < zmin) | (pfl.zef > zmax)

    t = 60*60*24*(pfl.UTCef - np.nanmin(pfl.UTCef))
    x = 1000.*(pfl.dist_ef - np.nanmin(pfl.dist_ef))
    k = 2*np.pi/X
    l = 0.
    m = 2*np.pi/Z
    f = gsw.f(pfl.lat_start)
    N = np.mean(np.sqrt(pfl.N2_ref[~nope]))

    om = gw.omega(N, k, m, l, f)

    u = gw.u(x, 0., pfl.zef, t, phi_0, k, l, m, om, phase_0=phase_0)

    return u[~nope] - utils.nan_detrend(pfl.zef[~nope], pfl.U[~nope], deg)

def w_model(params, pfl, zlim, deg):

    phi_0, X, Z, phase_0 = params
    zmin, zmax = zlim
    nope = np.isnan(pfl.z) | (pfl.z < zmin) | (pfl.z > zmax)

    t = 60*60*24*(pfl.UTC - np.nanmin(pfl.UTC))
    x = 1000.*(pfl.dist_ctd - np.nanmin(pfl.dist_ctd))

    k = 2*np.pi/X
    l = 0.
    m = 2*np.pi/Z
    f = gsw.f(pfl.lat_start)
    N = np.mean(np.sqrt(pfl.N2_ref[~nope]))
    om = gw.omega(N, k, m, l, f)

    w = gw.w(x, 0., pfl.z, t, phi_0, k, l, m, om, N, phase_0=phase_0)

    return w[~nope] - pfl.Ww[~nope]


def b_model(params, pfl, zlim, deg):

    phi_0, X, Z, phase_0 = params
    zmin, zmax = zlim
    nope = np.isnan(pfl.z) | (pfl.z < zmin) | (pfl.z > zmax)

    t = 60*60*24*(pfl.UTC - np.nanmin(pfl.UTC))
    x = 1000.*(pfl.dist_ctd - np.nanmin(pfl.dist_ctd))

    k = 2*np.pi/X
    l = 0.
    m = 2*np.pi/Z
    f = gsw.f(pfl.lat_start)
    N = np.mean(np.sqrt(pfl.N2_ref[~nope]))
    om = gw.omega(N, k, m, l, f)

    b = gw.b(x, 0., pfl.z, t, phi_0, k, l, m, om, N, phase_0=phase_0)

    resid = b[~nope] - pfl.b[~nope]

    return 250.*resid


def full_model(params, pfl, zlims, deg):
    return np.hstack((w_model(params, pfl, zlims, deg),
                      u_model(params, pfl, zlims, deg),
#                      v_model(params, pfl, zlims, deg),
                      b_model(params, pfl, zlims, deg)))


# The portions of the profiles that contain the wave. Found by eye.
zlims = {4976: {29: (-1600, -200),
                30: (-1000, -200),
                31: (-1600, -600),
                32: (-1600, -400)},
         4977: {24: (-1600, -200),
                25: (-1400, -600),
                26: (-1600, -600),
                27: (-1200, -200)}}

hpids_76 = np.array([29, 30, 31, 32])
hpids_77 = np.array([24, 25, 26, 27])


# Detrend degre
deg = 2


for Float, hpids in zip([E76, E77], [hpids_76, hpids_77]):
    for hpid in hpids:
        zlim = zlims[Float.floatID][hpid]
        pfl = Float.get_profiles(hpid)

        popt1, __, info, __, __ = op.leastsq(full_model, x0=[0.05, -2000., -2000., 0.],
                                             args=(pfl, zlim, deg),
                                             full_output=True)

        zmin, zmax = zlim
        nope = np.isnan(pfl.z) | (pfl.z < zmin) | (pfl.z > zmax)
        nope2 = np.isnan(pfl.zef) | (pfl.zef < zmin) | (pfl.zef > zmax)
        fig, axs = plt.subplots(1, 3, sharey=True)
        axs[0].plot(pfl.Ww, pfl.z, pfl.Ww[~nope] + w_model(popt1, pfl, zlim, deg),
                 pfl.z[~nope])
        axs[1].plot(utils.nan_detrend(pfl.zef[~nope2], pfl.U[~nope2], deg),
                    pfl.zef[~nope2],
                    utils.nan_detrend(pfl.zef[~nope2], pfl.U[~nope2], deg)
                    + u_model(popt1, pfl, zlim, deg),
                    pfl.zef[~nope2])
        axs[2].plot(250.*pfl.b, pfl.z, 250*pfl.b[~nope] + b_model(popt1, pfl, zlim, deg),
                 pfl.z[~nope])
        fig.suptitle("Float {}. hpid {:}.".format(pfl.floatID, pfl.hpid))

        print("Float {}. hpid {:}.".format(pfl.floatID, pfl.hpid))
        print("X = {:g}, Z = {:g}, phase = {:1.1f}".format(*popt1))

# %% EDIT TO ABOVE FITTING CODE:

def w_model(params, data):

    phi_0, X, Y, Z, phase_0 = params

    time, x, y, z, U, V, N, f, w_data, u_data, v_data, b_data = data

    k = 2*np.pi/X
    l = 2*np.pi/Y
    m = 2*np.pi/Z

    om = gw.omega(N, k, m, l, f)

    w = gw.w(x, y, z, time, phi_0, k, l, m, om, N, U=U, V=V, phase_0=phase_0)

    return w


def u_model(params, data):

    phi_0, X, Y, Z, phase_0 = params

    time, x, y, z, U, V, N, f, w_data, u_data, v_data, b_data = data

    k = 2*np.pi/X
    l = 2*np.pi/Y
    m = 2*np.pi/Z

    om = gw.omega(N, k, m, l, f)

    u = gw.u(x, y, z, time, phi_0, k, l, m, om, f=f, U=U, V=V, phase_0=phase_0)

    return u


def v_model(params, data):

    phi_0, X, Y, Z, phase_0 = params

    time, x, y, z, U, V, N, f, w_data, u_data, v_data, b_data = data

    k = 2*np.pi/X
    l = 2*np.pi/Y
    m = 2*np.pi/Z

    om = gw.omega(N, k, m, l, f)

    v = gw.v(x, y, z, time, phi_0, k, l, m, om, f=f, U=U, V=V, phase_0=phase_0)

    return v


def b_model(params, data):

    phi_0, X, Y, Z, phase_0 = params

    time, x, y, z, U, V, N, f, w_data, u_data, v_data, b_data = data

    k = 2*np.pi/X
    l = 2*np.pi/Y
    m = 2*np.pi/Z

    om = gw.omega(N, k, m, l, f)

    b = gw.b(x, y, z, time, phi_0, k, l, m, om, N, U=U, V=V, phase_0=phase_0)

    return 250.*b


def full_model(params, data):

    phi_0, X, Y, Z, phase_0 = params

    time, x, y, z, U, V, N, f, w_data, u_data, v_data, b_data = data

    k = 2*np.pi/X
    l = 2*np.pi/Y
    m = 2*np.pi/Z

    om = gw.omega(N, k, m, l, f)

    u = gw.u(x, y, z, time, phi_0, k, l, m, om, f=f, U=U, V=V, phase_0=phase_0)
    v = gw.v(x, y, z, time, phi_0, k, l, m, om, f=f, U=U, V=V, phase_0=phase_0)
    w = gw.w(x, y, z, time, phi_0, k, l, m, om, N, U=U, V=V, phase_0=phase_0)
    b = gw.b(x, y, z, time, phi_0, k, l, m, om, N, U=U, V=V, phase_0=phase_0)

    return np.hstack((w - w_data, u - u_data, v - v_data, b - b_data))


# Previously this looked like E76.get_timeseries([31, 32], ) etc. and the below
# bits of code were uncommented.

Float = E76
hpids = [26]
detrend = 1
zlims = (-1600, -600)

time, z = Float.get_timeseries(hpids, 'z')
__, dist = Float.get_timeseries(hpids, 'dist_ctd')
timeef, u_data = Float.get_timeseries(hpids, 'U_abs')
__, v_data = Float.get_timeseries(hpids, 'V_abs')
__, w_data = Float.get_timeseries(hpids, 'Ww')
__, b_data = Float.get_timeseries(hpids, 'b')
__, N2 = Float.get_timeseries(hpids, 'N2_ref')
__, x = Float.get_timeseries(hpids, 'x_ctd')
__, y = Float.get_timeseries(hpids, 'y_ctd')


nope = (z > zlims[1]) | (z < zlims[0])

time = time[~nope]
dist = dist[~nope]
z = z[~nope]
w_data = w_data[~nope]
b_data = b_data[~nope]
x = x[~nope]
y = y[~nope]

N = np.nanmean(np.sqrt(N2))
f = gsw.f(-57.5)

Unope = np.isnan(u_data)

timeef = timeef[~Unope]
u_data = u_data[~Unope]
v_data = v_data[~Unope]

u_data = np.interp(time, timeef, u_data)
v_data = np.interp(time, timeef, v_data)

U = np.nanmean(u_data)
V = np.nanmean(v_data)

if len(hpids) == 1:
    u_data = utils.nan_detrend(time, u_data, detrend)
    v_data = utils.nan_detrend(time, v_data, detrend)
else:
    t_split = Float.get_profiles(hpids[0]).UTC_end

    u_data[time > t_split] = utils.nan_detrend(z[time > t_split],
                                               u_data[time > t_split], detrend)
    u_data[time < t_split] = utils.nan_detrend(z[time < t_split],
                                               u_data[time < t_split], detrend)
    u_data[u_data > 0.3] = 0.

    v_data[time > t_split] = utils.nan_detrend(z[time > t_split],
                                               v_data[time > t_split], detrend)
    v_data[time < t_split] = utils.nan_detrend(z[time < t_split],
                                               v_data[time < t_split], detrend)


time *= 60.*60.*24
dist *= 1000.
time -= np.min(time)
dist -= np.min(dist)

data = [time, x, y, z, U, V, N, f, w_data, u_data, v_data, b_data]

L = 10000
popt = []

for i in xrange(L):

    print(i)

    x0 = np.hstack((np.abs(0.1*np.random.randn()), 5e3*np.random.randn(3),
                    1 + 0.5*np.random.randn()))
    popt1, __, info, __, __ = op.leastsq(full_model, x0=x0, args=(data),
                                         full_output=True)

    popt.append(popt1)

popt = np.asarray(popt)

fig, axs = plt.subplots(1, 4, figsize=(10, 6))
axs[0].hist(popt[:, 0], bins=np.arange(-40000, 40000, 1000),
            normed=False, log=False, alpha=0.8)
axs[1].hist(popt[:, 1], bins=np.arange(-40000, 40000, 1000),
            normed=False, log=False, alpha=0.8)
axs[2].hist(popt[:, 2], bins=np.arange(-40000, 40000, 1000),
            normed=False, log=False, alpha=0.8)
axs[3].hist(popt[:, 3], bins=np.arange(-40000, 40000, 1000),
            normed=False, log=False, alpha=0.8)

#popt1 = [-2000., -2000., -2000., 0.]

popt1 = np.median(popt, axis=0)

fig, axs = plt.subplots(4, 1)
axs[0].plot(time, u_data, time, u_model(popt1, data))
axs[1].plot(time, v_data, time, v_model(popt1, data))
axs[2].plot(time, w_data, time, w_model(popt1, data))
axs[3].plot(time, 250.*b_data, time, b_model(popt1, data))

fig, axs = plt.subplots(1, 4, sharey=True)
axs[0].plot(u_data, z, u_model(popt1, data), z)
axs[1].plot(v_data, z, v_model(popt1, data), z)
axs[2].plot(w_data, z, w_model(popt1, data), z)
axs[3].plot(250.*b_data, z, b_model(popt1, data), z)
