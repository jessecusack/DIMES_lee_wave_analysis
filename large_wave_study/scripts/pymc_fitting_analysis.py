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
import glob
import pickle

import gsw
import triangle
import pymc

lib_path = os.path.abspath('../modules')
if lib_path not in sys.path:
    sys.path.append(lib_path)

lib_path = os.path.abspath('../../ocean-tools')
if lib_path not in sys.path:
    sys.path.append(lib_path)

import emapex
import utils
import gravity_waves as gw
import plotting_functions as pf


try:
    print("Floats {} and {} exist!.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.load(4976)
    E77 = emapex.load(4977)

E76.calculate_pressure_perturbation()
E77.calculate_pressure_perturbation()

# %% Script params.

# Figure save path.
sdir = '../figures/pymc_fitting_analysis'
if not os.path.exists(sdir):
    os.makedirs(sdir)
# Universal figure font size.
matplotlib.rc('font', **{'size': 8})


# %% Fitting to profiles

def w_model(params, data):

    phi_0, X, Y, Z, phase_0 = params

    time, x, y, z, U, V, N, f = data

    k = 2*np.pi/X
    l = 2*np.pi/Y
    m = 2*np.pi/Z

    om = gw.omega(N, k, m, l, f)

    w = gw.w(x, y, z, time, phi_0, k, l, m, om, N, U=U, V=V, phase_0=phase_0)

    return w


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

    return b


def p_model(params, data):

    phi_0, X, Y, Z, phase_0 = params

    t, x, y, z, U, V, N, f = data

    k = 2*np.pi/X
    l = 2*np.pi/Y
    m = 2*np.pi/Z

    om = gw.omega(N, k, m, l, f)

    p = gw.phi(x, y, z, t, phi_0, k, l, m, om, U=U, V=V, phase_0=phase_0)

    return p


def full_model(params, data):
    return np.hstack((u_model(params, data),
                      v_model(params, data),
                      w_model(params, data),
                      b_model(params, data),
                      p_model(params, data)))


# %% Combined plots.
# Rewrite this for new trace files!
data_dir = '/noc/users/jc3e13/storage/processed/'
#M1 = pymc.database.pickle.load(os.path.join(data_dir, 'trace_4976_31_X1000_Y1000_Z1000.p'))
M2 = pymc.database.pickle.load(os.path.join(data_dir, 'trace_4976_32_X-1000_Y-1000_Z-1000.p'))
M3 = pymc.database.pickle.load(os.path.join(data_dir, 'trace_4977_26_X-1000_Y-1000_Z-1000.p'))
#M4 = pymc.database.pickle.load(os.path.join(data_dir, 'trace_4977_27_X1000_Y-1000_Z-1000.p'))

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

Ms = [M2, M3]

fig = plt.figure(figsize=(6.5, 3))

gs = gridspec.GridSpec(1, 5, width_ratios=[3, 1, 1, 1, 1])
gs.update(wspace=0.9)
axs = [plt.subplot(gs[0]), plt.subplot(gs[1]), plt.subplot(gs[2]),
       plt.subplot(gs[3]), plt.subplot(gs[4])]
#for ax in axs[1:]:
#    ax.yaxis.tick_right()
#    ax.yaxis.set_ticks_position('both')
#    ax.yaxis.set_label_position('right')

colors = ['blue', 'green', 'grey']

N = 2.2e-3
f = gsw.f(-57.5)

for i, M in enumerate(Ms):

    phi_0 = M.trace('phi_0')[:]
    k = np.pi*2./M.trace('X')[:]
    l = np.pi*2./M.trace('Y')[:]
    m = np.pi*2./M.trace('Z')[:]
    print("complete wavelength = {}".format(np.mean(np.pi*2/np.sqrt(k**2 + l**2 + m**2))))
    om = gw.omega(N, k, m, l, f)
    print("om/N = {}".format(np.mean(om/N)))
    w_0 = gw.W_0(phi_0, m, om, N)
    Edens = gw.Edens(w_0, k, m, l)
    Efluxz = gw.Efluxz(w_0, k, m, N, l, f)
    Mfluxz = gw.Mfluxz(phi_0, k, l, m, om, N)

    colprops={'color': colors[i]}
    data = [k, l, m]
    axs[0].boxplot(data, boxprops=colprops,
                   whiskerprops=colprops, capprops=colprops,
                   medianprops=colprops, showfliers=False,
                   labels=['$k$', '$l$', '$m$'])
    axs[1].boxplot(om/N, boxprops=colprops,
                   whiskerprops=colprops, capprops=colprops,
                   medianprops=colprops, showfliers=False,
                   labels=[''])
    axs[2].boxplot(Edens, boxprops=colprops,
                   whiskerprops=colprops, capprops=colprops,
                   medianprops=colprops, showfliers=False,
                   labels=[''])
    axs[3].boxplot(Efluxz, boxprops=colprops,
                   whiskerprops=colprops, capprops=colprops,
                   medianprops=colprops, showfliers=False,
                   labels=[""])
    axs[4].boxplot(Mfluxz, boxprops=colprops,
                   whiskerprops=colprops, capprops=colprops,
                   medianprops=colprops, showfliers=False,
                   labels=[''])

axs[0].hlines(0., *axs[0].get_xlim())
axs[0].set_ylim(-0.01, 0.01)

#ax0t = axs[0].twinx()
#ax0t.yaxis.tick_right()
#ax0t.yaxis.set_label_position('right')
#ax0t.set_ylabel('wavelength (m)')
#ax0t.set_ylim(axs[0].get_ylim())
#yticklabels_1 = np.array([500, 750, 1000, 1500, 5000])
#yticklabels = np.hstack((yticklabels_1, np.flipud(yticklabels_1)))
#yticks = np.pi*2./yticklabels
#yticks[:yticks.size/2] *= -1.
#ax0t.set_yticks(yticks)
#ax0t.set_yticklabels(yticklabels)
#ax0t.grid()

# Frequency log scale
axs[1].set_yscale('log')
axs[1].set_ylim(0.01, 10.)
axs[1].hlines(np.abs(f)/N, *axs[1].get_xlim())
axs[1].annotate('$f$', xy=(0.5, np.abs(f)/N))
axs[1].hlines(1., *axs[1].get_xlim())
axs[1].annotate('$N$', xy=(0.5, 1))


# The legend
labels = ['model fit to 4976 P 32', 'model fit to 4977 P 26', 'combined observations']
yloc = [0.009, 0.008, 0.007, 0.003]
for i in xrange(3):
    colprops = {'color': colors[i]}
    axs[0].text(0.53, yloc[i], labels[i], fontdict=colprops)

labelpad = -2
axs[0].set_ylabel('wavenumber (rad m$^{-1}$)', labelpad=labelpad)
axs[0].grid(axis='y')
plt.setp(axs[0].get_xticklabels(), fontsize=12)
axs[1].set_ylabel('Frequency $\omega/N$ (-)', labelpad=labelpad)
axs[1].grid(axis='y')
axs[2].set_ylabel('Energy density $E$ (J m$^{-3}$)', labelpad=labelpad)
axs[2].set_ylim(0., 30.)
axs[2].grid(axis='y')
axs[3].set_ylabel(r"Vertical energy flux $\overline{p'w'}$ (W m$^{-2}$)", labelpad=labelpad)
axs[3].grid(axis='y')
axs[3].set_ylim(0., 3.)
axs[4].set_ylabel('Vertical momentum flux $F_M^{(z)}$ (N m$^{-2}$)', labelpad=labelpad)
axs[4].grid(axis='y')
axs[4].set_ylim(0., 12.)

# Now to add on the observations.
with open('/noc/users/jc3e13/storage/processed/Edens.p', 'rb') as f:
    obsEdens = pickle.load(f)
with open('/noc/users/jc3e13/storage/processed/Eflux.p', 'rb') as f:
    obsEflux = pickle.load(f)
with open('/noc/users/jc3e13/storage/processed/Mflux.p', 'rb') as f:
    obsMflux = pickle.load(f)

obsEdens = np.hstack((obsEdens[0, 1, :], obsEdens[1, 0, :]))
obsEflux = np.hstack((obsEflux[0, 1, :], obsEflux[1, 0, :]))
obsMflux = np.hstack((obsMflux[0, 1, :], obsMflux[1, 0, :]))

colprops = {'color':'grey', 'alpha':0.5}

axs[2].boxplot(obsEdens, boxprops=colprops,
               whiskerprops=colprops, capprops=colprops,
               medianprops=colprops, showfliers=False,
               labels=[''])
axs[3].boxplot(obsEflux, boxprops=colprops,
               whiskerprops=colprops, capprops=colprops,
               medianprops=colprops, showfliers=False,
               labels=[''])
axs[4].boxplot(obsMflux, boxprops=colprops,
               whiskerprops=colprops, capprops=colprops,
               medianprops=colprops, showfliers=False,
               labels=[''])

amps_31 = np.array([[0.17, 0.11, 0.17, 4e-4],
                    [-0.12, -0.08, -0.025, -4e-4],
                    [0.2, 0.045, 0.07, -3.6e-4]])
amps_32 = np.array([[0.1, 0.18, 0.19, 2e-4],
                    [-0.06, -0.2, -0.2, -5e-4],
                    [-0.05, 0.22, 0.1, 5e-4],
                    [-0.17, 0., -0.2, -6e-4]])
amps_26 = np.array([[-0.095, -0.075, 0.16, 4.5e-4],
                    [0.2, 0.16, -0.1, -4e-4],
                    [-0.1, -0.04, 0.17, 3e-4],
                    [0.11, 0.17, -0.21, -6e-4]])
amps_27 = np.array([[-0.13, 0., -0.08, -3e-4],
                    [0.1, 0., 0.12, 2.7e-4],
                    [-0.03, 0., -0.07, 2.3e-4]])

amps = np.vstack((amps_31, amps_32, amps_26, amps_27))

u_amp = amps[:, 0]
v_amp = amps[:, 1]
w_amp = amps[:, 2]

alpha2 = w_amp**2/(u_amp**2 + v_amp**2)
obsom = np.sqrt(alpha2*N**2/(1 + alpha2))

axs[1].boxplot(obsom/N, boxprops=colprops,
               whiskerprops=colprops, capprops=colprops,
               medianprops=colprops, showfliers=False,
               labels=[''])

pf.my_savefig(fig, 'both', 'wavenumber_boxplots', sdir, ftype=('png', 'pdf'),
              fsize='double_col')


# %%

E76_hpids = [31, 32]
E77_hpids = [26, 27]

pfls = np.hstack((E76.get_profiles(E76_hpids), E77.get_profiles(E77_hpids)))

fig, axm = plt.subplots(len(pfls), 4, sharey='row', sharex='col',
                        figsize=(6.5, 7))
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
        axs[2].plot(100.*w_model(params, data), z, color='red', alpha=0.03)
        axs[3].plot(10000.*b_model(params, data), z, color='red', alpha=0.03)

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

# %% Figure for a profile

M_files = glob.glob('/noc/users/jc3e13/storage/processed/trace_497*.p')

info = {32: (-1600, -400),
        26: (-1600, -650)}

samples = 10000000
burn = 9800000
thin = 10

for M_file in M_files:

    M = pymc.database.pickle.load(M_file)
    savename = os.path.basename(M_file)[6:-2]
    hpid = int(M_file.split('_')[2])
    zmin, zmax = info[hpid]

    if hpid == 32:
        Float = E76
    elif hpid == 26:
        Float = E77

    time, z = Float.get_timeseries([hpid], 'z')
    timeef, U = Float.get_timeseries([hpid], 'U_abs')
    __, V = Float.get_timeseries([hpid], 'V_abs')
    __, W = Float.get_timeseries([hpid], 'Ww')
    __, B = Float.get_timeseries([hpid], 'b')
    __, N2 = Float.get_timeseries([hpid], 'N2_ref')
    __, x = Float.get_timeseries([hpid], 'x_ctd')
    __, y = Float.get_timeseries([hpid], 'y_ctd')
    __, PP = Float.get_timeseries([hpid], 'Pprime')

    # Correct for sign error in pp for downward profiles with large aspect
    # ratio. w is proxy for aspect.
    if (np.max(np.abs(W)) > 0.05) and (hpid % 2 != 0.):
        PP *= -1.
    else:
        PP *= 1.5

    nope = (z > zmax) | (z < zmin)

    time = time[~nope]
    W = W[~nope]
    B = B[~nope]
    PP = PP[~nope]
    x = x[~nope]
    y = y[~nope]
    z = z[~nope]

    N = np.nanmean(np.sqrt(N2))
    f = gsw.f(Float.get_profiles(hpid).lat_start)

    Unope = np.isnan(U)

    timeef = timeef[~Unope]
    U = U[~Unope]
    V = V[~Unope]

    U = np.interp(time, timeef, U)
    V = np.interp(time, timeef, V)

    Umean = np.mean(U)
    Vmean = np.mean(V)

    U = utils.nan_detrend(z, U, 1)
    V = utils.nan_detrend(z, V, 1)
    PP = utils.nan_detrend(z, PP, 1.)

    time *= 60.*60.*24
    time -= np.min(time)

    data = [time, x, y, z, Umean, Vmean, N, f]

    # Analysis and plotting.
    phi_0 = M.trace('phi_0')[:]
    k = np.pi*2/M.trace('X')[:]
    l = np.pi*2/M.trace('Y')[:]
    m = np.pi*2/M.trace('Z')[:]

    om = gw.omega(N, k, m, l, f)
    om_eulerian = om + k*Umean + l*Vmean
    w_0 = gw.W_0(phi_0, m, om, N)
    Efluxz = gw.Efluxz(w_0, k, m, N, l, f)
    Mfluxz = gw.Mfluxz(phi_0, k, l, m, om, N)

    stdf = 3.

    print("Mean frequency: {} +/- {}\n"
          "Eulerian frequency: {} +/- {}\n"
          "Mean X: {} +/- {}\n"
          "Mean Y: {} +/- {}\n"
          "Mean Z: {} +/- {}\n"
          "Mean k: {} +/- {}\n"
          "Mean l: {} +/- {}\n"
          "Mean m: {} +/- {}\n"
          "Mean phi_0: {} +/- {}\n"
          "Mean vertical energy flux: {} +/- {}\n"
          "Mean vertical momentum flux: {} +/- {}\n"
          "".format(np.mean(om), stdf*np.std(om),
                    np.mean(om_eulerian), stdf*np.std(om_eulerian),
                    np.mean(M.trace('X')[:]), stdf*np.std(M.trace('X')[:]),
                    np.mean(M.trace('Y')[:]), stdf*np.std(M.trace('Y')[:]),
                    np.mean(M.trace('Z')[:]), stdf*np.std(M.trace('Z')[:]),
                    np.mean(k), stdf*np.std(k),
                    np.mean(l), stdf*np.std(l),
                    np.mean(m), stdf*np.std(m),
                    np.mean(M.trace('phi_0')[:]), stdf*np.std(M.trace('phi_0')[:]),
                    np.mean(Efluxz),stdf* np.std(Efluxz),
                    np.mean(Mfluxz), stdf*np.std(Mfluxz)))

    # Plot fit comparison.
    fig, axs = plt.subplots(1, 5, sharey=True, figsize=(3.125, 3))

    z = z.copy()/1000.

    axs[0].plot(100.*U, z, color='black')
    axs[0].set_xlabel('$u^\prime$\n(cm s$^{-1}$)')
    axs[0].set_ylabel('$z$ (km)')
    axs[1].plot(100.*V, z, color='black')
    axs[1].set_xlabel('$v^\prime$\n(cm s$^{-1}$)')
    axs[2].plot(100.*W, z, color='black')
    axs[2].set_xlabel('$w^\prime$\n(cm s$^{-1}$)')
    axs[3].plot(10000.*B, z, color='black')
    axs[3].set_xlabel('$b^\prime$ ($10^{-4}$\nm s$^{-2}$)')
    axs[4].plot(100.*PP, z, color='black')
    axs[4].set_xlabel('$p^\prime$ ($10^{-2}$\nm$^2$ s$^{-2}$)')

    Ns = (samples - burn)/thin

    for i in xrange(0, Ns, 200):
        params = [M.trace('phi_0')[i], M.trace('X')[i], M.trace('Y')[i],
                  M.trace('Z')[i], M.trace('phase')[i]]
        axs[0].plot(100.*u_model(params, data), z, color='red', alpha=0.03)
        axs[1].plot(100.*v_model(params, data), z, color='red', alpha=0.03)
        axs[2].plot(100.*w_model(params, data), z, color='red', alpha=0.03)
        axs[3].plot(10000.*b_model(params, data), z, color='red',
                    alpha=0.03)
        axs[4].plot(100.*p_model(params, data), z, color='red', alpha=0.03)

    for ax in axs:
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=60)
        ax.locator_params(axis='x', tight=True, nbins=6)

    pf.my_savefig(fig, savename, 'profiles_fit', sdir, ftype={'png', 'pdf'},
                  fsize='single_col')
    plt.close()

#    triangle.corner(np.transpose(np.asarray([M.trace('X')[:],
#                                             M.trace('Y')[:],
#                                             M.trace('Z')[:],
#                                             M.trace('phi_0')[:]])),
#                    labels=['$\lambda_x$ (m)', '$\lambda_y$ (m)',
#                            '$\lambda_z$ (m)', '$\phi_0$ (m$^2$ s$^{-2}$)'])

    print("Completed " + savename)

# %% Figure for a profile

Ms = [M2, M3]
hpids = [32, 26]

info = {32: (-1600, -400),
        26: (-1600, -650)}

samples = 10000000
burn = 9800000
thin = 10

fig, axs = plt.subplots(2, 5, sharey=True, sharex='col', figsize=(3.125, 3))
#fig.tight_layout()

for i in xrange(2):

    M = Ms[i]
    hpid = hpids[i]
    zmin, zmax = info[hpid]

    if hpid == 32:
        Float = E76
    elif hpid == 26:
        Float = E77

    time, z = Float.get_timeseries([hpid], 'z')
    timeef, U = Float.get_timeseries([hpid], 'U_abs')
    __, V = Float.get_timeseries([hpid], 'V_abs')
    __, W = Float.get_timeseries([hpid], 'Ww')
    __, B = Float.get_timeseries([hpid], 'b')
    __, N2 = Float.get_timeseries([hpid], 'N2_ref')
    __, x = Float.get_timeseries([hpid], 'x_ctd')
    __, y = Float.get_timeseries([hpid], 'y_ctd')
    __, PP = Float.get_timeseries([hpid], 'Pprime')

    # Correct for sign error in pp for downward profiles with large aspect
    # ratio. w is proxy for aspect.
    if (np.max(np.abs(W)) > 0.05) and (hpid % 2 != 0.):
        PP *= -1.
    else:
        PP *= 1.5

    nope = (z > zmax) | (z < zmin)

    time = time[~nope]
    W = W[~nope]
    B = B[~nope]
    PP = PP[~nope]
    x = x[~nope]
    y = y[~nope]
    z = z[~nope]

    N = np.nanmean(np.sqrt(N2))
    f = gsw.f(Float.get_profiles(hpid).lat_start)

    Unope = np.isnan(U)

    timeef = timeef[~Unope]
    U = U[~Unope]
    V = V[~Unope]

    U = np.interp(time, timeef, U)
    V = np.interp(time, timeef, V)

    Umean = np.mean(U)
    Vmean = np.mean(V)

    U = utils.nan_detrend(z, U, 1)
    V = utils.nan_detrend(z, V, 1)
    PP = utils.nan_detrend(z, PP, 1.)

    time *= 60.*60.*24
    time -= np.min(time)

    data = [time, x, y, z, Umean, Vmean, N, f]

    # Analysis and plotting.
    phi_0 = M.trace('phi_0')[:]
    k = np.pi*2/M.trace('X')[:]
    l = np.pi*2/M.trace('Y')[:]
    m = np.pi*2/M.trace('Z')[:]

    om = gw.omega(N, k, m, l, f)
    om_eulerian = om + k*Umean + l*Vmean
    w_0 = gw.W_0(phi_0, m, om, N)
    Efluxz = gw.Efluxz(w_0, k, m, N, l, f)
    Mfluxz = gw.Mfluxz(phi_0, k, l, m, om, N)

    stdf = 3.

    print("Mean frequency: {} +/- {}\n"
          "Eulerian frequency: {} +/- {}\n"
          "Mean X: {} +/- {}\n"
          "Mean Y: {} +/- {}\n"
          "Mean Z: {} +/- {}\n"
          "Mean k: {} +/- {}\n"
          "Mean l: {} +/- {}\n"
          "Mean m: {} +/- {}\n"
          "Mean phi_0: {} +/- {}\n"
          "Mean vertical energy flux: {} +/- {}\n"
          "Mean vertical momentum flux: {} +/- {}\n"
          "".format(np.mean(om), stdf*np.std(om),
                    np.mean(om_eulerian), stdf*np.std(om_eulerian),
                    np.mean(M.trace('X')[:]), stdf*np.std(M.trace('X')[:]),
                    np.mean(M.trace('Y')[:]), stdf*np.std(M.trace('Y')[:]),
                    np.mean(M.trace('Z')[:]), stdf*np.std(M.trace('Z')[:]),
                    np.mean(k), stdf*np.std(k),
                    np.mean(l), stdf*np.std(l),
                    np.mean(m), stdf*np.std(m),
                    np.mean(M.trace('phi_0')[:]), stdf*np.std(M.trace('phi_0')[:]),
                    np.mean(Efluxz),stdf* np.std(Efluxz),
                    np.mean(Mfluxz), stdf*np.std(Mfluxz)))

    # Plot fit comparison.
    z = z.copy()/1000.

    axs[i, 0].set_ylabel('$z$ (km)')
    axs[i, 0].set_title(str(Float.floatID) + ' P ' + str(hpid))

    Ns = (samples - burn)/thin

    for j in xrange(0, Ns, 200):
        params = [M.trace('phi_0')[j], M.trace('X')[j], M.trace('Y')[j],
                  M.trace('Z')[j], M.trace('phase')[j]]
        if j == 0:
            label = 'fit'
        else:
            label = None
        alpha = 0.9
        axs[i, 0].plot(100.*u_model(params, data), z, color='red', alpha=alpha, label=label)
        axs[i, 1].plot(100.*v_model(params, data), z, color='red', alpha=alpha, label=label)
        axs[i, 2].plot(100.*w_model(params, data), z, color='red', alpha=alpha, label=label)
        axs[i, 3].plot(10000.*b_model(params, data), z, color='red', alpha=alpha, label=label)
        axs[i, 4].plot(100.*p_model(params, data), z, color='red', alpha=alpha, label=label)

    label = 'obs'
    alpha = 0.9
    axs[i, 0].plot(100.*U, z, color='black', alpha=alpha, label=label)
    axs[i, 1].plot(100.*V, z, color='black', alpha=alpha, label=label)
    axs[i, 2].plot(100.*W, z, color='black', alpha=alpha, label=label)
    axs[i, 3].plot(10000.*B, z, color='black', alpha=alpha, label=label)
    axs[i, 4].plot(100.*PP, z, color='black', alpha=alpha, label=label)

axs[1, 0].set_xlabel('$u^\prime$\n(cm s$^{-1}$)')
axs[1, 0].set_xlim(-30., 30.)
axs[1, 0].set_xticks([-20., 0., 20.])
axs[1, 1].set_xlabel('$v^\prime$\n(cm s$^{-1}$)')
axs[1, 1].set_xlim(-30., 30.)
axs[1, 1].set_xticks([-20., 0., 20.])
axs[1, 2].set_xlabel('$w^\prime$\n(cm s$^{-1}$)')
axs[1, 2].set_xlim(-30., 30.)
axs[1, 2].set_xticks([-20., 0., 20.])
axs[1, 3].set_xlabel('$b^\prime$ ($10^{-4}$\nm s$^{-2}$)')
axs[1, 3].set_xlim(-10., 10.)
axs[1, 3].set_xticks([-6., 0., 6.])
axs[1, 4].set_xlabel('$p^\prime$ ($10^{-2}$\nm$^2$ s$^{-2}$)')
axs[1, 4].set_xlim(-6., 6.)
axs[1, 4].set_xticks([-4., 0., 4.])

axs[1, 4].legend(loc='upper right', bbox_to_anchor=(1, 1.22), fontsize=7,
                 ncol=2)

for ax in axs[1, :]:
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=60, fontsize=6)
    ax.locator_params(axis='x', tight=True)

pf.my_savefig(fig, '32_26', 'profiles_fit', sdir, ftype={'png', 'pdf'},
              fsize='single_col')