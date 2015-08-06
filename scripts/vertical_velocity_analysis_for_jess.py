# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 11:40:08 2015

@author: jc3e13
"""

# Import the necessary modules, rename as required:
import numpy as np
import scipy as sp
import scipy.optimize as op
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import sys
import datetime
from glob import glob

import gsw  # Thermodynamic equation of state for seawater.
import triangle  # A cool plotting function.
import pymc  # MCMC toolbox.

lib_path = os.path.abspath('../modules')
if lib_path not in sys.path:
    sys.path.append(lib_path)

import utils

# %% SCRIPT PARAMS

do_plots = False

# %% LOAD DATA

dec_dir = '/noc/users/jc3e13/storage/jess_silvester/4980b/dec/'
ctd_files = np.sort(glob(os.path.join(dec_dir, 'ema-4980b-*-ctd.mat')))
gps_files = np.sort(glob(os.path.join(dec_dir, 'ema-4980b-*-gps.mat')))
mis_files = np.sort(glob(os.path.join(dec_dir, 'ema-4980b-*-mis.mat')))

# First run through files, deduce data size and grab important info
nprofs = int(os.path.basename(ctd_files[-1]).split('-')[2])
hpids = np.NaN*np.zeros(nprofs)
nobss = np.NaN*np.zeros(nprofs)
ngpsobss = np.NaN*np.zeros(nprofs)

ctd_vars = ['hpid', 'nobs']

for ctd_file in ctd_files:
    ctd_info = sp.io.loadmat(ctd_file, variable_names=ctd_vars,
                             squeeze_me=True)
    hpid = ctd_info['hpid']
    hpids[hpid-1] = ctd_info['hpid']
    nobss[hpid-1] = ctd_info['nobs']

for gps_file in gps_files:
    gps_info = sp.io.loadmat(gps_file, variable_names=ctd_vars,
                             squeeze_me=True)
    hpid = gps_info['hpid']
    ngpsobss[hpid-1] = gps_info['nobs']

# Initialise arrays.
mobs = np.nanmax(nobss)
mgpsobs = np.nanmax(ngpsobss)

P = np.NaN*np.zeros((mobs, nprofs))
T = np.NaN*np.zeros((mobs, nprofs))
S = np.NaN*np.zeros((mobs, nprofs))
pc = np.NaN*np.zeros((mobs, nprofs))
UXT = np.NaN*np.zeros((mobs, nprofs))

LAT = np.NaN*np.zeros((mgpsobs, nprofs))
LON = np.NaN*np.zeros((mgpsobs, nprofs))
STAT = np.NaN*np.zeros((mgpsobs, nprofs))
UXT_GPS = np.NaN*np.zeros((mgpsobs, nprofs))
UXT_APF9 = np.NaN*np.zeros((mgpsobs, nprofs))

FloatAlpha = np.NaN*np.zeros(nprofs)
FloatBeta = np.NaN*np.zeros(nprofs)
FloatMass = np.NaN*np.zeros(nprofs)
PistonCountsPerCC = np.NaN*np.zeros(nprofs)
PistonParkPosition = np.NaN*np.zeros(nprofs)
PressureBallastPoint = np.NaN*np.zeros(nprofs)

# Second run through files and grab data.
gps_vars = ['hpid', 'LAT', 'LON', 'STAT', 'UXT_GPS', 'UXT_APF9']
mis_vars = ['hpid', 'FloatAlpha', 'FloatBeta', 'FloatMass',
            'PistonCountsPerCC', 'PistonParkPosition',
            'PressureBallastPoint']
# Load CTD data.
for ctd_file in ctd_files:
    ctd_data = sp.io.loadmat(ctd_file, squeeze_me=True)
    hpid = ctd_data['hpid']
    idx = hpid - 1

    P[:nobss[idx], idx] = ctd_data['P']
    T[:nobss[idx], idx] = ctd_data['T']
    S[:nobss[idx], idx] = ctd_data['S']
    pc[:nobss[idx], idx] = ctd_data['pc']
    UXT[:nobss[idx], idx] = ctd_data['UXT']
# Load GPS data.
for i, gps_file in enumerate(gps_files):
    gps_data = sp.io.loadmat(gps_file, variable_names=gps_vars,
                             squeeze_me=True)
    hpid = gps_data['hpid']
    idx = hpid - 1

    LAT[:ngpsobss[idx], idx] = gps_data['LAT']
    LON[:ngpsobss[idx], idx] = gps_data['LON']
    STAT[:ngpsobss[idx], idx] = gps_data['STAT']
    UXT_GPS[:ngpsobss[idx], idx] = gps_data['UXT_GPS']
    UXT_APF9[:ngpsobss[idx], idx] = gps_data['UXT_APF9']
# Load mission data.
for i, mis_file in enumerate(mis_files):
    mis_info = sp.io.loadmat(mis_file, variable_names=mis_vars,
                             squeeze_me=True)
    hpid = mis_info['hpid']
    idx = hpid - 1

    FloatAlpha[idx] = mis_info['FloatAlpha']
    FloatBeta[idx] = mis_info['FloatBeta']
    FloatMass[idx] = mis_info['FloatMass']
    PistonCountsPerCC[idx] = mis_info['PistonCountsPerCC']
    PistonParkPosition[idx] = mis_info['PistonParkPosition']
    PressureBallastPoint[idx] = mis_info['PressureBallastPoint']

# %% ADDITIONAL PROCESSING
# Time of first measurement for each half profile.
UXT_START = UXT[0, :]
DT_START = []
for i in xrange(nprofs):
    try:
        DT_START.append(datetime.datetime.fromtimestamp(UXT_START[i]))
    except ValueError:
        DT_START.append(np.NaN)

DT_START = np.array(DT_START)


def make_timeseries(t, v):
    times = t[:, :].flatten(order='F')
    vals = v[:, :].flatten(order='F')
    nnans = ~np.isnan(times) & ~np.isnan(vals)
    return times[nnans], vals[nnans]

# There must be some offset between GPS time and the float time. Use the float
# time when interpolating float data.
tsUXT_GPS, tsLON = make_timeseries(UXT_APF9, LON)
__, tsLAT = make_timeseries(UXT_APF9, LAT)

LON_START = utils.nan_interp(UXT_START, tsUXT_GPS, tsLON)
LAT_START = utils.nan_interp(UXT_START, tsUXT_GPS, tsLAT)

if do_plots:
    plt.figure()
    plt.plot(LON_START, LAT_START, 'kx', markersize=2)
    plt.plot(LON, LAT)

# Use the TOES-10 toolbox to calculate height, absolute salinity and density.
z = gsw.z_from_p(P, LAT_START)
SA = gsw.SA_from_SP(S, P, LON_START, LAT_START)
CT = gsw.CT_from_t(SA, T, P)
rho = gsw.rho_t_exact(SA, T, P)
rho_0 = gsw.pot_rho_t_exact(SA, T, P, p_ref=0)
N2_mid, __ = gsw.Nsquared(SA, CT, P, LAT_START)
UXT_mid = (UXT[1:, :] + UXT[:-1, :])/2.

# Speed of the float in m s-1.
w_float = np.NaN*np.zeros((mobs, nprofs))
N2 = np.NaN*np.zeros((mobs, nprofs))
for i in xrange(nprofs):
    try:
        w_float[:, i] = utils.finite_diff(UXT[:, i], z[:, i])
        N2[:, i] = utils.nan_interp(UXT[:, i], UXT_mid[:, i], N2_mid[:, i])
    except ValueError:
        continue

if do_plots:
    plt.figure()
    plt.plot(w_float, z, 'k', alpha=0.1)

# %% OPTIMISE STEADY MODEL


def steady_model(params, data, fixed):
    """Calculates and returns the vertical velocity that the float would have
    if it were moving in still water.

    params:

    0: V_0 = 1.  # Float volume when neutrally buoyant [m3].
    1: CA = 1.  # Combination of drag and surface area [m2].
    2: alpha_p = 3.76e-6  # Coeff. of expansion with pressure [-].
    3: p_0 = 2000.  # Pressure when float is neutrally buoyant [dbar].
    4: alpha_ppos = 1.156e-6  # Coeff. of expansion with piston position [m3].
    5: ppos_0 = 16.  # Piston position when neutrally buoyant [-].
    6: M = 27.179  # Float mass [kg].

    data:

    ppos, p, rho

    fixed:

    List of values to fix with the same numbering as parameters. Use None for
    varying parameters.

    Gravity is given a value of 9.8 m s-2.

    """

    ppos, p, rho = data

    g = -9.8  # Gravitational acceleration [m s-2].

    # Apply fixed value for those parameters that are fixed.
    for i, val in enumerate(fixed):
        if val is not None:
            params[i] = val

    V_0, CA, alpha_p, p_0, alpha_ppos, ppos_0, M = params

    # Float volume
    V = V_0*(1 + alpha_p*(p - p_0)) + alpha_ppos*(ppos - ppos_0)

    return np.sign(rho*V - M)*np.sqrt(np.abs(g*(M - rho*V))/(rho*CA))


def cost(params, data, fixed, wf):
    ws = steady_model(params, data, fixed)
    return (ws - wf)**2


def fitter(ppos, P, rho, w_f, params0, fixed, Plims=(50., 1500.)):
    """This function takes an EM-APEX float, fits a vertical velocity model
    using the given arguments, estimates errors using bootstrapping technique.

    Parameters
    ----------
    Float : EMApexFloat object
        The float to fit.
    hpids : array
        Half profiles to optimise model for.
    params0 : array
        The variable values that need regridding.

    Returns
    -------
    wfi : Bunch object
        Class containing all fitting information.

    Notes
    -----

    """

    still_water_model = steady_model

    Pmin, Pmax = Plims

    use = (P > Pmin) & (P < Pmax)

    data = [ppos[use], P[use], rho[use]]
    w_f = w_f[use]

    cargs = (data, fixed, w_f)

    p, pcov, info, mesg, ier = op.leastsq(cost, params0, args=cargs,
                                          full_output=True)

    # Bootstrapping.
    ps = []
    # 200 random data sets are generated and fitted
    for i in range(300):
        rand_idx = np.random.rand(*w_f.shape) < 0.25
        rand_data = [d[rand_idx] for d in data]
        cargs = (rand_data, fixed, w_f[rand_idx])
        rand_params, __ = op.leastsq(cost, params0, args=cargs)
        ps.append(rand_params)

    ps = np.array(ps)
    pcov = np.cov(ps.T)
    pmean = np.mean(ps, 0)
    pcorr = np.corrcoef(ps.T)

    wfi = utils.Bunch(params0=params0,
                      fixed=fixed,
                      model_func=still_water_model,
                      Plims=Plims,
                      p=p,
                      ps=ps,
                      pmean=pmean,
                      pcov=pcov,
                      pcorr=pcorr,
                      info=info,
                      mesg=mesg,
                      ier=ier)

    return wfi

# Initial parameter estimates from mission data and buest guesses.
V_0 = 2.67e-2  # m3
CA = 3.412e-2  # Drag multiplied by cross sectional area m2
alpha_p = -3.401e-6  # np.nanmean(FloatAlpha)  # Coefficient of expansion w/ pressure.
p_0 = 1500. #  np.nanmean(PressureBallastPoint)
alpha_ppos = 1e-6/np.nanmean(PistonCountsPerCC)  # Convert units to m3
ppos_0 = 9.  # np.nanmean(PistonParkPosition)
M = np.nanmean(FloatMass)/1000.  # Convert to kg.

# Perform the optimisation on given range of hpids.
use_hpids = np.arange(110, 210, 1)
idxs = [i for i in xrange(len(hpids)) if hpids[i] in use_hpids]
idxs = np.array(idxs)

time, pc_fit = make_timeseries(UXT[:, idxs], pc[:, idxs])
__, P_fit = make_timeseries(UXT[:, idxs], P[:, idxs])
__, rho_fit = make_timeseries(UXT[:, idxs], rho[:, idxs])
__, w_f = make_timeseries(UXT[:, idxs], w_float[:, idxs])
params0 = np.array([V_0, CA, alpha_p, p_0, alpha_ppos, ppos_0, M])
fixed = [None, CA, alpha_p, p_0, alpha_ppos, ppos_0, M]
Plims = (350, 850)

# Do the fit!
fit_info = fitter(pc_fit, P_fit, rho_fit, w_f, params0, fixed, Plims)

w_steady = steady_model(fit_info.p, [pc, P, rho], fit_info.fixed)

# This little correction is necessary because of array index problems later.
# It isn't a good fix but it'll do for now.
w_steady[0, 148] = np.NaN

w = w_float - w_steady

# %% FIT ASSESSMENT


def assess_w_fit(floatID, t, Ww, Ws, Wz, wfi):
    """ """

    font = {'family': 'normal',
            'weight': 'normal',
            'size': 8}
    matplotlib.rc('font', **font)

    # Histogram of vertical water velocity.

    nans = np.isnan(Ww)

    Ww_mean = np.nanmean(Ww)
    Ww_std = np.nanstd(Ww)

    plt.figure(figsize=(3, 3))
    bins = np.arange(-0.15, 0.155, 0.005)
    Ww_hist, bins, patches = plt.hist(Ww[~nans], bins=bins,
                                      histtype='stepfilled')
    plt.setp(patches, 'facecolor', 'b', 'alpha', 0.75)
    plt.xlim(np.min(bins), np.max(bins))
    plt.xlabel('$W_w$ (m s$^{-1}$)')
    plt.xticks(rotation=45)
    title_str = ("Float {}\nmean = {:1.2e} m s$^{{-1}}$\nstd = {:1.2e} "
                 "m s$^{{-1}}$").format(floatID, Ww_mean, Ww_std)
    plt.title(title_str)
#    name = save_id + ud_id + '_ww_histogram.pdf'
#    fname = os.path.join(save_dir, name)
#    plt.savefig(fname, format='pdf', bbox_inches='tight')

    # Time series of different velocity measures.
    plt.figure(figsize=(6, 3))
    plt.plot(t, Ww)
    plt.plot(t, Wz)
    plt.plot(t, Ws)
    plt.ylabel('$W_w$, $W_f$, $W_s$ (m s$^{-1}$)')
    plt.xlabel('Time')
    plt.xticks(rotation=45)
    title_str = "Float {}".format(floatID)
    plt.title(title_str)
    plt.legend(['$W_w$', '$W_f$', '$W_s$'])
#    name = save_id + ud_id + '_ww_wf_w0_timeseries.pdf'
#    fname = os.path.join(save_dir, name)
#    plt.savefig(fname, format='pdf', bbox_inches='tight')

    # Parameter estimates and correlations.
    pnames = np.array(['$V_0$', '$CA$', r'$\alpha_p$', '$p_0$', r'$\alpha_k$',
                       '$k_0$', '$M$'])
    N = len(pnames)
    ticks = np.arange(0.5, N, 1)

    plt.figure(figsize=(3, 3))
    plt.pcolormesh(np.flipud(wfi.pcorr), cmap=plt.get_cmap('PiYG'))
    cbar = plt.colorbar()
    cbar.set_label('Correlation')
    plt.clim(-1, 1)
    plt.xticks(ticks, pnames)
    plt.yticks(ticks, pnames[::-1])
    title_str = ("Float {}").format(floatID)
    plt.title(title_str)
#    name = save_id + '_param_corr.pdf'
#    fname = os.path.join(save_dir, name)
#    plt.savefig(fname, format='pdf', bbox_inches='tight')

    not_fixed = np.array([(p is None) for p in wfi.fixed])
    ps = wfi.ps[:, not_fixed]
    p = wfi.p[not_fixed]
    params0 = wfi.params0[not_fixed]
    triangle.corner(ps, labels=pnames[not_fixed])
    f = plt.gcf()
    f.set_size_inches(7, 7)
    axs = f.axes
    N = np.shape(ps)[1]

    formatter = ticker.ScalarFormatter(useOffset=False)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 2))

    for i in xrange(N):
        for j in xrange(N):
            idx = i*N + j
            if i == N - 1:
                axs[idx].xaxis.set_major_formatter(formatter)
            if (j == 0) and (i > 0):
                axs[idx].yaxis.set_major_formatter(formatter)
            if i == j:
                axs[idx].vlines(p[i], *axs[idx].get_ylim(), color='r')
                axs[idx].vlines(params0[i], *axs[idx].get_ylim(),
                                color='g')
#
#    name = save_id + '_param_matrix_scatter.pdf'
#    fname = os.path.join(save_dir, name)
#    plt.savefig(fname, format='pdf', bbox_inches='tight')

    plt.draw()
    plt.show()

#    name = save_id + ud_id + '_fit_info.p'
#    fname = os.path.join(save_dir, name)
#    with open(fname, 'wb') as f:
#        pickle.dump(wfi, f)


def assess_profile(use_hpids, hpids, DT_START, z, S, T, pc, rho_0, N2,
                   w_steady, w_float, w, save=False, save_dir=None):

    idxs = [i for i in xrange(len(hpids)) if hpids[i] in use_hpids]
    idxs = np.array(idxs)

    if save and save_dir is None:
        save_dir = '/noc/users/jc3e13/storage/jess_silvester/figures'

    for idx in idxs:

        fig1, axs1 = plt.subplots(1, 7, figsize=(15, 3.5), sharey=True)

        axs1[0].set_title('Profile: {:1.0f}'.format(hpids[idx]))
        if type(DT_START[idx]) is datetime.datetime:
            axs1[-1].set_title(DT_START[idx].strftime("%Y-%m-%d %H:%M"))
        axs1[0].set_ylim(-1000., 0.)

        axs1[0].plot(100.*w[:, idx], z[:, idx], color='black')
        axs1[0].vlines(0., *axs1[0].get_ylim(), color='black')
        axs1[1].plot(100.*w_steady[:, idx], z[:, idx], color='black')
        axs1[1].plot(100.*w_float[:, idx], z[:, idx], linestyle=':',
                     color='black')
        axs1[2].plot(pc[:, idx], z[:, idx], color='black')
        axs1[3].plot(S[:, idx], z[:, idx], color='black')
        axs1[4].plot(T[:, idx], z[:, idx], color='black')
        axs1[4].vlines(0., *axs1[0].get_ylim(), color='black')
        axs1[5].plot(rho_0[:, idx]-1000., z[:, idx], color='black')
        axs1[6].plot(1e6*N2[:, idx], z[:, idx], color='black')
        axs1[6].vlines(0., *axs1[0].get_ylim(), color='black')

        axs1[0].set_ylabel('Height (m)')
        axs1[0].set_xlabel('Vertical velocity (cm s$^{-1}$)')
        axs1[0].set_xlim(-10., 10.)
        axs1[1].set_xlabel('Float vertical velocity (cm s$^{-1}$)')
        axs1[2].set_xlabel('Piston counts (-)')
        axs1[3].set_xlabel('Practical salinity (-)')
        axs1[3].set_xlim(34.1, 34.8)
        axs1[4].set_xlabel('In-situ temperature ($^\circ$C)')
        axs1[4].set_xlim(-0.2, 1.2)
        axs1[5].set_xlabel('Potential density (kg m$^{-3}$)')
        axs1[5].set_xlim(27.35, 27.85)
        axs1[6].set_xlabel('N$^{2}$ ($10^{-6}$ rad$^{2}$ s$^{-2}$)')
        axs1[6].set_xlim(-10., 40.)

        for ax in axs1:
            ax.grid()
            plt.setp(ax.xaxis.get_majorticklabels(), rotation=60)

        if save:
            save_name = "profile_{:1.0f}.png".format(hpids[idx])
            save_path = os.path.join(save_dir, save_name)
            plt.savefig(save_path, bbox_inches='tight')
            plt.close()


wmax = np.NaN*np.zeros(nprofs)
for i in xrange(nprofs):
    try:
        use = (P[:, i] < 800) & (P[:, i] > 200)
        wmax[i] = np.nanmax(w[use, i])
    except ValueError:
        continue


time, Ww = make_timeseries(UXT[:, idxs], w[:, idxs])
__, Ws = make_timeseries(UXT[:, idxs], w_steady[:, idxs])
__, Wz = make_timeseries(UXT[:, idxs], w_float[:, idxs])

assess_w_fit('4980b', time, Ww, Ws, Wz, fit_info)

# %% MCMC


def pymc_fitter(ppos, P, rho, w_f, Plims=(350, 850.), samples=10000):

    Pmin, Pmax = Plims
    use = (P > Pmin) & (P < Pmax)

    data = [ppos[use], P[use], rho[use]]

    w_f = w_f[use]

    def model():

        # Priors.
        V_0 = pymc.Uniform('V_0', 0., 1, value=2.7e-2)
        CA = pymc.Uniform('CA', 0., 1., value=3.412e-2)
        alpha_p = pymc.Uniform('alpha_p', -1, 0., value=-3.76e-6)
        p_0 = pymc.Uniform('p_0', 0., 4000., value=1500.)
        alpha_ppos = 1.156e-6
        ppos_0 = pymc.Uniform('ppos_0', -100, 300., value=16.)
        M = 27.75
        fixed = 7*[None]
        sig = pymc.Uniform('sig', 0., 0.1, value=0.02)

        @pymc.deterministic()
        def float_model(V_0=V_0, CA=CA, alpha_p=alpha_p, p_0=p_0,
                        alpha_ppos=alpha_ppos, ppos_0=ppos_0, M=M):
            params = [V_0, CA, alpha_p, p_0, alpha_ppos, ppos_0, M]
            return steady_model(params, data, fixed)

        # Likelihood
        y = pymc.Normal('y', mu=float_model, tau=1./sig**2, value=w_f,
                        observed=True)

        return locals()

    M = pymc.MCMC(model(), db='pickle',
                  dbname='/noc/users/jc3e13/storage/jess_silvester/pymc_trace.p')
    burn = samples/2
    thin = 5
    M.sample(samples, burn, thin)
    pymc.Matplot.plot(M, common_scale=False)

    chain = np.asarray([M.trace('V_0')[:],
                        M.trace('CA')[:],
                        M.trace('alpha_p')[:],
                        M.trace('p_0')[:],
#                        M.trace('alpha_ppos')[:],
                        M.trace('ppos_0')[:]])
    labels = [r'$V_0$', r'$C_D^*$', r'$\alpha_p$', r'$p_0$', #r'$\alpha_k$',
              r'$k_0$']
    triangle.corner(np.transpose(chain), labels=labels)

    return M

# Perform the optimisation on given range of hpids.
use_hpids = np.arange(110, 210, 1)
idxs = [i for i in xrange(len(hpids)) if hpids[i] in use_hpids]
idxs = np.array(idxs)

time, pc_fit = make_timeseries(UXT[:, idxs], pc[:, idxs])
__, P_fit = make_timeseries(UXT[:, idxs], P[:, idxs])
__, rho_fit = make_timeseries(UXT[:, idxs], rho[:, idxs])
__, w_f = make_timeseries(UXT[:, idxs], w_float[:, idxs])

Plims = (350, 850)

# Do the fit!
M = pymc_fitter(pc_fit, P_fit, rho_fit, w_f, Plims=Plims, samples=1000000)

#w_steady = steady_model(fit_info.p, [pc, P, rho], fit_info.fixed)

# This little correction is necessary because of array index problems later.
# It isn't a good fix but it'll do for now.
#w_steady[0, 148] = np.NaN