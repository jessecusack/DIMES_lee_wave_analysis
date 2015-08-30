# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 12:58:56 2014

@author: jc3e13
"""

import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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


try:
    print("Floats {} and {} exist!.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.load(4976)
    E77 = emapex.load(4977)

# %% Script params.

# Figure save path.
sdir = '../figures/pymc_fitting_analysis'
if not os.path.exists(sdir):
    os.makedirs(sdir)
# Universal figure font size.
matplotlib.rc('font', **{'size': 9})


# %% Fitting to profiles

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


# %% Combined plots.
# Rewrite this for new trace files!
M1 = pymc.database.pickle.load('/noc/users/jc3e13/storage/processed/trace_31_C.p')
M2 = pymc.database.pickle.load('/noc/users/jc3e13/storage/processed/trace_32_C.p')
M3 = pymc.database.pickle.load('/noc/users/jc3e13/storage/processed/trace_26_C.p')
M4 = pymc.database.pickle.load('/noc/users/jc3e13/storage/processed/trace_27_C.p')

# %% Post-load.

X = np.hstack((M1.trace('X')[:], M2.trace('X')[:], M3.trace('X')[:], M4.trace('X')))
Y = np.hstack((M1.trace('Y')[:], M2.trace('Y')[:], M3.trace('Y')[:], M4.trace('Y')))
Z = np.hstack((M1.trace('Z')[:], M2.trace('Z')[:], M3.trace('Z')[:], M4.trace('Z')))
phi_0 = np.hstack((M1.trace('phi_0')[:], M2.trace('phi_0')[:],
                   M3.trace('phi_0')[:], M4.trace('phi_0')))

triangle.corner(np.transpose(np.asarray([X, Y, Z, phi_0])),
                labels=['$\lambda_x$ (m)', '$\lambda_y$ (m)',
                        '$\lambda_z$ (m)', '$\phi_0$ (m$^2$ s$^{-2}$)'])

fig = plt.gcf()
fig.set_size_inches(6.5, 6.5)
pf.my_savefig(fig, 'both', 'fit_histograms', sdir, ftype='pdf',
              fsize='double_col')

# %%

Ms = [M1, M2, M3, M4]

fig = plt.figure(figsize=(6.5, 3))

gs = gridspec.GridSpec(1, 2, width_ratios=[3,1])
gs.update(wspace=0.3)
axs = [plt.subplot(gs[0]), plt.subplot(gs[1])]
axs[1].yaxis.tick_right()
axs[1].yaxis.set_ticks_position('both')
axs[1].yaxis.set_label_position('right')

colors = ['blue', 'green', 'red', 'purple']

for i, M in enumerate(Ms):
    colprops={'color':colors[i]}
    data = [np.pi*2./M.trace('X')[:], np.pi*2./M.trace('Y')[:],
            np.pi*2./M.trace('Z')[:]]
    axs[0].boxplot(data, boxprops=colprops,
                   whiskerprops=colprops, capprops=colprops,
                   medianprops=colprops, showfliers=False,
                   labels=['$k$', '$l$', '$m$'])
    axs[1].boxplot(M.trace('phi_0')[:], boxprops=colprops,
                   whiskerprops=colprops, capprops=colprops,
                   medianprops=colprops, showfliers=False, labels=['$\phi_0$'])

ax0t = axs[0].twinx()
ax0t.yaxis.tick_right()
ax0t.yaxis.set_label_position('right')
ax0t.set_ylabel('wavelength (m)')
ax0t.set_ylim(axs[0].get_ylim())
yticklabels = np.array([1000, 1500, 2000, 3000, 5000, 10000])
yticks = -np.pi*2./yticklabels
ax0t.set_yticks(yticks)
ax0t.set_yticklabels(yticklabels)
ax0t.grid()

axs[0].set_ylabel('wavenumber (rad m$^{-1}$)')
axs[1].set_ylabel('pressure perturbation (m$^2$ s$^{-2}$)')

pf.my_savefig(fig, 'both', 'wavenumber_boxplots', sdir, ftype='pdf',
              fsize='double_col')


# %%

E76_hpids = [31, 32]
E77_hpids = [26, 27]

pfls = np.hstack((E76.get_profiles(E76_hpids), E77.get_profiles(E77_hpids)))

fig, axm = plt.subplots(len(pfls), 4, sharey='row', sharex='col',
                        figsize=(6.5, 10))
fig.subplots_adjust(hspace=0.05, wspace=0.1)
rot = 'vertical'
col = 'black'
deg = 2

U_var = 'U'
V_var = 'V'

ylims = {27: (-1000, -200),
         26: (-1540, -600),
         31: (-1474, -600),
         32: (-1580, -400)}

Ms = {26: M3,
      27: M4,
      31: M1,
      32: M2}


for pfl, axs in zip(pfls, axm):

    axs[0].set_ylabel('$z$ (m)')
    axs[0].plot(100.*utils.nan_detrend(pfl.zef, getattr(pfl, U_var), deg), pfl.zef, col)
    axs[1].plot(100.*utils.nan_detrend(pfl.zef, getattr(pfl, V_var), deg), pfl.zef, col)
    axs[2].plot(100.*pfl.Ww, pfl.z, col)
    axs[3].plot(10000.*pfl.b, pfl.z, col)
    axs[2].annotate("EM {}\nP {}".format(pfl.floatID, pfl.hpid[0]),
                    (-29, -250))

    for ax in axs:
        ax.vlines(0., *ax.get_ylim())
        ax.grid()
        ax.axhspan(*ylims[pfl.hpid[0]], color='grey', alpha=0.5)
        ax.set_ylim(-1600., 0.)

    Ns = len(M1.trace('X')[:])
    M = Ms[pfl.hpid[0]]

    time = pfl.UTC
    z = pfl.z
    timeef = pfl.UTCef
    U = pfl.U_abs
    V = pfl.V_abs
    N2 = pfl.N2_ref
    x = pfl.x_ctd
    y = pfl.y_ctd

    zlims = ylims[pfl.hpid[0]]
    nope = (z < zlims[0]) | (z > zlims[1]) | np.isnan(z)

    time = time[~nope]
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

    Umean = np.mean(U)
    Vmean = np.mean(V)

    time *= 60.*60.*24
    time -= np.min(time)

    data = [time, x, y, z, Umean, Vmean, N, f]

    for i in xrange(0, Ns, 40):
        params = [M.trace('phi_0')[i], M.trace('X')[i], M.trace('Y')[i],
                  M.trace('Z')[i], M.trace('phase')[i]]
        axs[0].plot(100.*u_model(params, data), z, color='red', alpha=0.03)
        axs[1].plot(100.*v_model(params, data), z, color='red', alpha=0.03)
        axs[2].plot(100.*w_model(params, data)/wscale, z, color='red', alpha=0.03)
        axs[3].plot(10000.*b_model(params, data)/bscale, z, color='red', alpha=0.03)

for ax in axm[-1, :]:
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=rot)

axm[-1, 0].set_xlabel('$u$ (cm s$^{-1}$)')
axm[-1, 0].set_xlim(-30, 30)
axm[-1, 1].set_xlabel('$v$ (cm s$^{-1}$)')
axm[-1, 1].set_xlim(-30, 30)
axm[-1, 2].set_xlabel('$w$ (cm s$^{-1}$)')
axm[-1, 2].set_xlim(-30, 30)
axm[-1, 3].set_xlabel('$b$ ($10^{-4}$ m s$^{-2}$)')

pf.my_savefig(fig, 'both', 'profiles_fit', sdir, ftype='pdf',
              fsize='double_col')