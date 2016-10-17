# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 14:53:08 2016

@author: jc3e13
"""

import os
import numpy as np
from scipy.stats import binned_statistic
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import emapex
import utils
from my_savefig import my_savefig

try:
    print("Floats {} and {} exist!.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.load(4976)
#    E76.generate_regular_grids(zmin=zmin, dz=dz)
    E77 = emapex.load(4977)

# %% Script params.

# Figure save path.
sdir = '../figures/inertial_energy'
if not os.path.exists(sdir):
    os.makedirs(sdir)
# Universal figure font size.
matplotlib.rc('font', **{'size': 8})


# %%
min_hpid = 3
max_hpid = 100
hpid = np.arange(min_hpid, max_hpid+1)
down_hpid = np.arange(min_hpid, max_hpid, 2)
up_hpid = np.arange(min_hpid+1, max_hpid+2, 2)

zmin = -1500.
zmax = -50.
dzbin = 10.
zbins = np.arange(zmin, zmax + dzbin, dzbin)
zmid = (zbins[:-1] + zbins[1:])/2.
Nbins = len(zbins) - 1
rho0 = 1025.

fig, axs = plt.subplots(2, 1, sharex='col', figsize=(6.5, 5))

for i, Float in enumerate([E76, E77]):
    __, idxs = Float.get_profiles(hpid, ret_idxs=True)
    __, idxs_d = Float.get_profiles(down_hpid, ret_idxs=True)
    __, idxs_u = Float.get_profiles(up_hpid, ret_idxs=True)

    U = Float.U_abs[:, idxs]
    V = Float.V_abs[:, idxs]
    z = Float.zef[:, idxs]
    t = Float.UTC_start[idxs]
    tp = (t[:-2] + t[2:])/2
    tp = utils.datenum_to_datetime(tp)

    Np = len(idxs)

    Ubin = np.NaN*np.zeros((Nbins, Np))
    Vbin = np.NaN*np.zeros((Nbins, Np))

    for j in xrange(len(idxs)):
        Ubin[:, j], __, __ = binned_statistic(z[:, j], U[:, j],
                                              statistic=np.nanmean, bins=zbins)
        Vbin[:, j], __, __ = binned_statistic(z[:, j], V[:, j],
                                              statistic=np.nanmean, bins=zbins)

    Uin = np.NaN*np.zeros((Nbins, Np-2))
    Vin = np.NaN*np.zeros((Nbins, Np-2))
    Uin[:, :] = (Ubin[:, :-2] - Ubin[:, 2:])/2.
    Vin[:, :] = (Vbin[:, :-2] - Vbin[:, 2:])/2.
    KEh = rho0*(Uin**2. + Vin**2)/2

    Uinm = np.ma.masked_invalid(Uin)
    Vinm = np.ma.masked_invalid(Vin)
    KEhm = np.ma.masked_invalid(KEh)

#    axs[0].contourf(tp, zmid, Uinm, 30, cmap=plt.get_cmap('bwr'), vmin=-0.25, vmax=0.25)
#    axs[0+i].contourf(tp, zmid, Vinm, 30, cmap=plt.get_cmap('bwr'), vmin=-0.25, vmax=0.25)
    vmax = 20.
    c = axs[0+i].pcolormesh(tp, zmid, KEhm, cmap=plt.get_cmap('viridis'), vmin=0., vmax=vmax)

fig.subplots_adjust(right=0.82)
cbar_ax = fig.add_axes([0.85, 0.2, 0.02, 0.6])
cbar = fig.colorbar(c, cax=cbar_ax, extend='max')
cbar.set_clim(0., vmax)
cbar.set_label('Horiz. Inertial KE (J m$^{-3}$)')
myFmt = mdates.DateFormatter('%j')
axs[-1].xaxis.set_major_formatter(myFmt)
axs[-1].set_xlabel('Julien day')

my_savefig(fig, 'KE_comparison', sdir, ftype=('png', 'pdf'), fsize='double_col')