#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 15:21:22 2017

@author: jc3e13
"""

import glob
import os
import numpy as np
import matplotlib
import emapex
import pickle
import matplotlib.pyplot as plt
from tabulate import tabulate
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import sandwell
import cmocean
import utils

psdir = '../processed_data/'
fsd = '../figures/w_sections'
bf = os.path.abspath(glob.glob('/noc/users/jc3e13/storage/smith_sandwell/topo_*.img')[0])

# Universal figure font size.
matplotlib.rc('font', **{'size': 8})

save_figs = True


fits = {}
for fid in emapex.FIDS_DIMES:
    filename = os.path.join(psdir, "{}_wfi.p".format(fid))
    with open(filename, 'rb') as f:
        fits[fid] = pickle.load(f)

# %%

Nf = len(emapex.FIDS_DIMES)

headers = ['FID', 'V0', 'CDu', 'CDd', 'alpha_p', 'alpha_k']

data_stack = np.empty((Nf, len(headers)))

for i in range(Nf):
    fid = emapex.FIDS_DIMES[i]
    data_stack[i, 0] = fid
    data_stack[i, 1:] = fits[fid]['p'][[0, 1, 2, 3, 5]]

print(tabulate(data_stack, headers, tablefmt='latex', floatfmt='.3g'))

# %%
N = 6
IDS = [4976, 4977, 6478, 6480, 6625, 6626]
FLOATS = [emapex.EMApexFloat(emapex.find_file(ID), ID) for ID in IDS]

for i in range(N):
    fid = IDS[i]
    filename = os.path.join(psdir, "{}_wfi.p".format(fid))
    with open(filename, 'rb') as f:
        wfi = pickle.load(f)
        FLOATS[i].apply_w_model(wfi)

E4976, E4977, E6478, E6480, E6625, E6626 = FLOATS

# %%
N = 2
IDS_ = [4596, 4814]
FLOATS_ALT = [emapex.EMApexFloat(emapex.find_file(ID), ID) for ID in IDS_]

for i in range(N):
    fid = IDS_[i]
    filename = os.path.join(psdir, "{}_wfi.p".format(fid))
    with open(filename, 'rb') as f:
        wfi = pickle.load(f)
        FLOATS_ALT[i].apply_w_model(wfi)

E4596, E4814 = FLOATS_ALT

# %% map
lon_lat = [-70, -30, -65, -45]
proj = ccrs.PlateCarree()

###############################################################################
lons, lats, bathy = sandwell.read_grid(lon_lat, bf)
bathy = -1*np.ma.masked_where(bathy > 0, bathy)

fig = plt.figure(figsize=(6.5, 5))
ax = plt.subplot(111, projection=proj)
ax.set_extent(lon_lat, ccrs.PlateCarree())
cax = fig.add_axes([0.6, 0.45, 0.2, 0.01])

ax.contour(lons, lats, bathy.data, [0.], colors='k')

C = ax.pcolormesh(lons[::10, ::10], lats[::10, ::10], bathy[::10, ::10],
                  vmin=0., vmax=6000., cmap=cmocean.cm.deep, rasterized=True)
cbp = plt.colorbar(C, cax=cax, orientation='horizontal')
cbp.set_ticks([0., 3000., 6000])
cbp.set_label('Depth (m)')

for Float in FLOATS + [FLOATS_ALT[1]]:

    lon = Float.lon_start
    lat = Float.lat_start

    ax.plot(lon, lat, label=str(Float.floatID))

ax.legend(loc=0, framealpha=1.)

ax.set_xticks([-70, -65, -60, -55, -50, -45, -40, -35, -30], crs=proj)
ax.set_yticks([-65, -60, -55, -50, -45], crs=proj)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)

if save_figs:
    filename = os.path.join(fsd, 'floats_map.pdf')
    fig.savefig(filename, bbox_inches='tight', pad_inches=0)

# %% vertical velocity sections
step = 20

###############################################################################

for Float in FLOATS + FLOATS_ALT:

    __, z = Float.get_timeseries(Float.hpid, 'z')
    t, Ww = Float.get_timeseries(Float.hpid, 'Ww')
    __, d = Float.get_timeseries(Float.hpid, 'dist_ctd')
    __, Wf = Float.get_timeseries(Float.hpid, 'Wz')

    use = t > 7000.
    z = z[use]
    Ww = Ww[use]
    Wf = Wf[use]
    t = t[use]
    d = d[use]

    t = utils.datenum_to_datetime(t)

    # change to x = d or x = t. Don't forget to change axis label.
    x = d

    Ww[np.abs(Wf) < 0.08] = np.nan

    fig, ax = plt.subplots(1, 1, figsize=(6.5, 3))
    cax = fig.add_axes([0.75, 0.25, 0.1, 0.02])

    C = ax.scatter(x[::step], z[::step], s=5, c=Ww[::step], edgecolor='none',
                cmap=plt.get_cmap('bwr'), vmin=-.05, vmax=.05, rasterized=True)
    ax.set_ylim(np.min(z), np.max(z))
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Height (m)')

#    title_str = ("Float {}").format(Float.floatID)
#    ax.set_title(title_str)
    cbp = plt.colorbar(C, cax=cax, orientation='horizontal', extend='both')
    cbp.set_label('$w$ (m s$^{-1}$)')
    cbp.set_ticks([-0.05, 0, 0.05])
    cbp.set_ticklabels(['-0.05', '0', '0.05'])

    if save_figs:
        filename = os.path.join(fsd, '{}_w_section.pdf'.format(Float.floatID))
        fig.savefig(filename, bbox_inches='tight', pad_inches=0)

# %% Float 4596 topo section
Float = E4596

b = sandwell.interp_track(Float.lon_start, Float.lat_start, bf)

__, z = Float.get_timeseries(Float.hpid, 'z')
t, Ww = Float.get_timeseries(Float.hpid, 'Ww')
__, Wf = Float.get_timeseries(Float.hpid, 'Wz')
__, d = Float.get_timeseries(Float.hpid, 'dist_ctd')

use = t > 7000.
z = z[use]
Ww = Ww[use]
Wf = Wf[use]
d = d[use]

Ww[np.abs(Wf) < 0.08] = np.nan

fig, ax = plt.subplots(1, 1, figsize=(6.5, 3))
cax = fig.add_axes([0.75, 0.25, 0.1, 0.02])

ax.plot(Float.dist, b, 'k')

C = ax.scatter(d[::step], z[::step], s=5, c=Ww[::step], edgecolor='none',
            cmap=plt.get_cmap('bwr'), vmin=-.05, vmax=.05, rasterized=True)
#ax.set_ylim(np.min(z), np.max(z))
ax.set_xlim(np.min(d), np.max(d))
ax.set_xlabel('Distance (km)')
ax.set_ylabel('Height (m)')

title_str = ("Float {}").format(Float.floatID)
ax.set_title(title_str)
cbp = plt.colorbar(C, cax=cax, orientation='horizontal', extend='both')
cbp.set_label('$w$ (m s$^{-1}$)')
cbp.set_ticks([-0.05, 0, 0.05])