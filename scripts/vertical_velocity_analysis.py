# -*- coding: utf-8 -*-
"""
Created on Thu May 15 12:07:00 2014

@author: jc3e13
"""

import emapex
import vertical_velocity_model as vvm
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as mcolors
import pandas as pd
import os
import pickle

reload(emapex)
reload(vvm)


def assess_w_fit(Float, save_id=''):
    """ """

    font = {'family': 'normal',
            'weight': 'normal',
            'size': 8}
    matplotlib.rc('font', **font)

    wfi = Float.__wfi
    hpids = wfi.hpids
    floatID = Float.floatID
    if wfi.profiles == 'updown':
        updown = True
        ud_id = '_updown'
    else:
        updown = False
        ud_id = ''

    s = os.sep
    save_dir = '..'+s+'figures'+s+'vertical_velocity_analysis'

    # Histogram of vertical water velocity.
    Ww = Float.rWw.flatten(order='F')

    Ww_mean = np.nanmean(Ww)
    Ww_std = np.nanstd(Ww)

    plt.figure(figsize=(3, 3))
    bins = np.arange(-0.15, 0.155, 0.005)
    Ww_hist, bins, patches = plt.hist(Ww, bins=bins, histtype='stepfilled')
    plt.setp(patches, 'facecolor', 'b', 'alpha', 0.75)
    plt.xlim(np.min(bins), np.max(bins))
    plt.xlabel('$W_w$ (m s$^{-1}$)')
    plt.xticks(rotation=45)
    title_str = ("Float {}\nmean = {:1.2e} m s$^{{-1}}$\nstd = {:1.2e} "
                 "m s$^{{-1}}$").format(floatID, Ww_mean, Ww_std)
    plt.title(title_str)
    name = save_id + ud_id + '_ww_histogram.pdf'
    fname = os.path.join(save_dir, name)
    plt.savefig(fname, format='pdf', bbox_inches='tight')

    # Time series of different velocity measures.
    hpid1 = hpids[0]

    plt.figure(figsize=(6, 3))
    N = 4
    time, Ww = Float.get_timeseries(np.arange(hpid1, hpid1+N), 'rWw')
    __, Wz = Float.get_timeseries(np.arange(hpid1, hpid1+N), 'rWz')
    __, Ws = Float.get_timeseries(np.arange(hpid1, hpid1+N), 'rWs')
    plt.plot(time, Ww)
    plt.plot(time, Wz)
    plt.plot(time, Ws)
    plt.ylabel('$W_w$, $W_f$, $W_s$ (m s$^{-1}$)')
    plt.xlabel('Time')
    plt.xticks(rotation=45)
    title_str = ("Float {}, half profiles {}").format(floatID, hpids[0:N])
    plt.title(title_str)
    plt.legend(['$W_w$', '$W_f$', '$W_s$'])
    name = save_id + ud_id + '_ww_wf_w0_timeseries.pdf'
    fname = os.path.join(save_dir, name)
    plt.savefig(fname, format='pdf', bbox_inches='tight')

    # Scatter section of water velocity.

    cmap = mcolors.LinearSegmentedColormap.from_list(name='red_white_blue',
                                                     colors=[(0, 0, 1),
                                                             (1, 1, 1),
                                                             (1, 0, 0)],
                                                     N=40,
                                                     )

    z_vals = np.arange(-1500., 0., 10.)
    __, z, Ww = Float.get_interp_grid(np.arange(1, 600), z_vals, 'z', 'rWw')
    z = z.flatten(order='F')
    Ww = Ww.flatten(order='F')
    __, __, d = Float.get_interp_grid(np.arange(1, 600), z_vals, 'z',
                                      'dist_ctd')
    d = d.flatten(order='F')
    plt.figure(figsize=(6, 3))
    plt.scatter(d, z, c=Ww, edgecolor='none', cmap=cmap)
    plt.ylim(np.min(z), np.max(z))
    plt.xlim(np.min(d), np.max(d))
    plt.xlabel('Distance (km)')
    plt.ylabel('Depth (m)')
    plt.xlim(np.min(Float.dist), np.max(Float.dist))
    title_str = ("Float {}").format(floatID)
    plt.title(title_str)
    cbar = plt.colorbar(orientation='horizontal', extend='both')
    cbar.set_label('$W_w$ (m s$^{-1}$)')
    plt.clim(-0.15, 0.15)
    name = save_id + ud_id + '_ww_scatter_section.pdf'
    fname = os.path.join(save_dir, name)
    plt.savefig(fname, format='pdf', bbox_inches='tight')

#    # Varience of Ww2 with N2
#    Ww2 = Float.rWw.flatten(order='F')**2
#    N2 = Float.rN2.flatten(order='F')
#    Ww2 = Ww2[N2 > 0.]
#    N = np.sqrt(N2[N2 > 0.])
#    bins = np.linspace(0., np.max(N), 20)
#    idxs = np.digitize(N, bins)
#    Wmean = []
#    Wstd = []
#    Nm = (bins[1:] + bins[:-1])/2.
#    N0 = 1e-4
#    for i in xrange(len(bins) - 1):
#        Wmean.append(np.nanmean(Ww2[idxs == i+1]))
#        Wstd.append(np.nanstd(Ww2[idxs == i+1]))
#    plt.errorbar(Nm, Wmean, Wstd)
#    plt.plot(Nm, 0.25*N0/Nm)

    # Parameter estimates and correlations.
    pnames = ['$V_0$', '$CA$', r'$\alpha_p$', '$p_0$', r'$\alpha_k$', '$k_0$',
              '$M$']
    N = len(pnames)
    ticks = np.arange(0.5, N, 1)

    if updown:
        for n, ud in enumerate(['up', 'down']):
            plt.figure(figsize=(3, 3))
            plt.pcolormesh(np.flipud(wfi.pcorr[n]), cmap=plt.get_cmap('PiYG'))
            cbar = plt.colorbar()
            cbar.set_label('Correlation')
            plt.clim(-1, 1)
            plt.xticks(ticks, pnames)
            plt.yticks(ticks, pnames[::-1])
            title_str = ('Float {}' + ud).format(floatID)
            plt.title(title_str)
            name = save_id + '_' + ud + '_param_corr.pdf'
            fname = os.path.join(save_dir, name)
            plt.savefig(fname, format='pdf', bbox_inches='tight')

            pps = pd.DataFrame(wfi.ps[n], columns=pnames)
            axs = pd.tools.plotting.scatter_matrix(pps, hist_kwds={'bins': 12})
            f = plt.gcf()
            f.set_size_inches(7, 7)
            formatter = ticker.ScalarFormatter(useOffset=False)
            formatter.set_scientific(True)
            formatter.set_powerlimits((-1, 2))

            for i in xrange(N):
                for j in xrange(N):

                    axs[i, j].xaxis.set_major_formatter(formatter)
                    axs[i, j].yaxis.set_major_formatter(formatter)

                    if i == j:
                        y = np.array(axs[i, j].get_ylim())
                        x = np.array([wfi.p[n][i], wfi.p[n][i]])
                        axs[i, j].plot(x, y, 'r-')
                        x = np.array([wfi.params0[i], wfi.params0[i]])
                        axs[i, j].plot(x, y, 'g-')

            name = save_id + '_' + ud + '_param_matrix_scatter.pdf'
            fname = os.path.join(save_dir, name)
            plt.savefig(fname, format='pdf', bbox_inches='tight')

    else:

        plt.figure(figsize=(3, 3))
        plt.pcolormesh(np.flipud(wfi.pcorr), cmap=plt.get_cmap('PiYG'))
        cbar = plt.colorbar()
        cbar.set_label('Correlation')
        plt.clim(-1, 1)
        plt.xticks(ticks, pnames)
        plt.yticks(ticks, pnames[::-1])
        title_str = ("Float {}").format(floatID)
        plt.title(title_str)
        name = save_id + '_param_corr.pdf'
        fname = os.path.join(save_dir, name)
        plt.savefig(fname, format='pdf', bbox_inches='tight')

        pps = pd.DataFrame(wfi.ps, columns=pnames)
        axs = pd.tools.plotting.scatter_matrix(pps, hist_kwds={'bins': 12})
        f = plt.gcf()
        f.set_size_inches(7, 7)
        formatter = ticker.ScalarFormatter(useOffset=False)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-1, 2))

        for i in xrange(N):
            for j in xrange(N):

                axs[i, j].xaxis.set_major_formatter(formatter)
                axs[i, j].yaxis.set_major_formatter(formatter)

                if i == j:
                    y = np.array(axs[i, j].get_ylim())
                    x = np.array([wfi.p[i], wfi.p[i]])
                    axs[i, j].plot(x, y, 'r-')
                    x = np.array([wfi.params0[i], wfi.params0[i]])
                    axs[i, j].plot(x, y, 'g-')

        name = save_id + '_param_matrix_scatter.pdf'
        fname = os.path.join(save_dir, name)
        plt.savefig(fname, format='pdf', bbox_inches='tight')

    plt.draw()
    plt.show()

    name = save_id + ud_id + '_fit_info.p'
    fname = os.path.join(save_dir, name)
    with open(fname, 'wb') as f:
        pickle.dump(wfi, f)

###############################################################################

try:
    print("Floats {} and {}.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4976)
    E77 = emapex.EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4977)


# model = '1'
# cf_key = 'diffsq'
# params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
# fixed = [None, None, None, None, None, None, None]
# wfi = vvm.fitter(E76, params0, fixed, model=model, profiles='all',
#                 cf_key=cf_key)
# E76.apply_w_model(wfi)
# assess_w_fit(E76, str(E76.floatID))
# print(E76.__wfi.p)
#
# model = '1'
# cf_key = 'diffsq'
# params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
# fixed = [None, None, None, None, None, None, None]
# wfi = vvm.fitter(E77, params0, fixed, model=model, profiles='all',
#                 cf_key=cf_key)
# E77.apply_w_model(wfi)
# assess_w_fit(E77, str(E77.floatID))
# print(E77.__wfi.p)

###############################################################################

model = '1'
cf_key = 'diffsq'
params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
fixed = [None, None, None, None, None, None, 27.179]
wfi = vvm.fitter(E76, params0, fixed, model=model, profiles='all',
                 cf_key=cf_key)
E76.apply_w_model(wfi)
assess_w_fit(E76, str(E76.floatID)+'_fix_M')
print(E76.__wfi.p)

model = '1'
cf_key = 'diffsq'
params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
fixed = [None, None, None, None, None, None, 27.179]
wfi = vvm.fitter(E77, params0, fixed, model=model, profiles='all',
                 cf_key=cf_key)
E77.apply_w_model(wfi)
assess_w_fit(E77, str(E77.floatID)+'_fix_M')
print(E77.__wfi.p)

###############################################################################

model = '1'
cf_key = 'diffsq'
params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
fixed = [None, None, None, 2000., None, 16, 27.179]
wfi = vvm.fitter(E76, params0, fixed, model=model, profiles='all',
                 cf_key=cf_key)
E76.apply_w_model(wfi)
assess_w_fit(E76, str(E76.floatID)+'_fix_p0k0M')
print(E76.__wfi.p)

model = '1'
cf_key = 'diffsq'
params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
fixed = [None, None, None, 2000., None, 16, 27.179]
wfi = vvm.fitter(E77, params0, fixed, model=model, profiles='all',
                 cf_key=cf_key)
E77.apply_w_model(wfi)
assess_w_fit(E77, str(E77.floatID)+'_fix_p0k0M')
print(E77.__wfi.p)

###############################################################################

model = '1'
cf_key = 'diffsq'
params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
fixed = [None, None, None, None, 1.156e-6, None, 27.179]
wfi = vvm.fitter(E76, params0, fixed, model=model, profiles='all',
                 cf_key=cf_key)
E76.apply_w_model(wfi)
assess_w_fit(E76, str(E76.floatID)+'_fix_alphakM')
print(E76.__wfi.p)

model = '1'
cf_key = 'diffsq'
params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
fixed = [None, None, None, None, 1.156e-6, None, 27.179]
wfi = vvm.fitter(E77, params0, fixed, model=model, profiles='all',
                 cf_key=cf_key)
E77.apply_w_model(wfi)
assess_w_fit(E77, str(E77.floatID)+'_fix_alphakM')
print(E77.__wfi.p)

###############################################################################
#
# model = '1'
# cf_key = 'diffsq'
# params0 = np.array([3e-2, 5e-2, 3e-6, 2e+3, 1e-6, 16., 27.179])
# fixed = [None, None, None, 2e+3, None, None, 27.179]
# wfi = vvm.fitter(E76, params0, fixed, model=model, profiles='all',
#                 cf_key=cf_key)
# E76.apply_w_model(wfi)
# assess_w_fit(E76, str(E76.floatID)+'_fix_p0M')
# print(E76.__wfi.p)
#
# model = '1'
# cf_key = 'diffsq'
# params0 = np.array([3e-2, 5e-2, 2e-6, 2e+3, 1e-6, 16., 27.179])
# fixed = [None, None, None, 2e+3, None, None, 27.179]
# wfi = vvm.fitter(E77, params0, fixed, model=model, profiles='all',
#                 cf_key=cf_key)
# E77.apply_w_model(wfi)
# assess_w_fit(E77, str(E77.floatID)+'_fix_p0M')
# print(E77.__wfi.p)
#
###############################################################################
#
# model = '1'
# cf_key = 'diffsq'
# params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
# fixed = [None, None, None, None, None, None, None]
# wfi = vvm.fitter(E76, params0, fixed, model=model, profiles='updown',
#                 cf_key=cf_key)
# E76.apply_w_model(wfi)
# assess_w_fit(E76, str(E76.floatID))
# print(E76.__wfi.p)
#
# model = '1'
# cf_key = 'diffsq'
# params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
# fixed = [None, None, None, None, None, None, None]
# wfi = vvm.fitter(E77, params0, fixed, model=model, profiles='updown',
#                 cf_key=cf_key)
# E77.apply_w_model(wfi)
# assess_w_fit(E77, str(E77.floatID))
# print(E77.__wfi.p)
#
###############################################################################
#
# model = '1'
# cf_key = 'diffsq'
# params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
# fixed = [None, None, None, None, 1.156e-6, None, 27.179]
# wfi = vvm.fitter(E76, params0, fixed, model=model, profiles='updown',
#                 cf_key=cf_key)
# E76.apply_w_model(wfi)
# assess_w_fit(E76, str(E76.floatID)+'_fix_alphakM')
# print(E76.__wfi.p)
#
# model = '1'
# cf_key = 'diffsq'
# params0 = np.array([3e-2, 5e-2, 3e-6, 4e+2, 1e-6, 16., 27.179])
# fixed = [None, None, None, None, 1.156e-6, None, 27.179]
# wfi = vvm.fitter(E77, params0, fixed, model=model, profiles='updown',
#                 cf_key=cf_key)
# E77.apply_w_model(wfi)
# assess_w_fit(E77, str(E77.floatID)+'_fix_alphakM')
# print(E77.__wfi.p)