# -*- coding: utf-8 -*-
"""
Created on Thu Jul 24 11:00:03 2014

@author: jc3e13
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.optimize as op
import gsw
import os
import sys
from scipy.interpolate import griddata
from scipy.integrate import cumtrapz
import scipy.signal as sig
import mpl_toolkits.basemap as bm

lib_path = os.path.abspath('../modules')
if lib_path not in sys.path:
    sys.path.append(lib_path)

import emapex
import sandwell
import plotting_functions as pf


try:
    print("Floats {} and {} exist!.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.load(4976)
    E77 = emapex.load(4977)

plot_dir = '../figures/poster_plots'
ftype = 'png'
bwr = plt.get_cmap('bwr')
bwr1 = mcolors.LinearSegmentedColormap.from_list(name='red_white_blue',
                                                 colors=[(0, 0, 1),
                                                         (1, 1, 1),
                                                         (1, 0, 0)],
                                                 N=40)


def my_savefig(fid, fname):
    fname = str(fid) + '_' + fname
    fname = os.path.join(plot_dir, fname) + '.' + ftype
    plt.savefig(fname, bbox_inches='tight')

# %% Mountain centered coords

bwr = plt.get_cmap('bwr')
E76_hpids = np.arange(20, 40) # np.arange(31, 33)
E77_hpids = np.arange(15, 35) # np.arange(26, 28)
bathy_file = '/noc/users/jc3e13/storage/smith_sandwell/topo_17.1.img'
vars = ['Ww']#, 'U_abs', 'V_abs']
zvars = ['z']#, 'zef', 'zef']
dvars = ['dist_ctd']#, 'dist_ef', 'dist_ef']
texvars = ['$W_w$']#, '$U$', '$V$']
clims = [(-10., 10.)]#, (-100., 100.), (-100, 100.)]

var_1_vals = np.linspace(-40., 40., 80)
var_2_vals = np.linspace(-1500, 0, 500)
Xg, Zg = np.meshgrid(var_1_vals, var_2_vals)

Wgs = []
ds = []
zs = []

for Float, hpids in zip([E76, E77], [E76_hpids, E77_hpids]):

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)

    for var, zvar, dvar, texvar, clim in zip(vars, zvars, dvars, texvars,
                                             clims):

        V = getattr(Float, var)[:, idxs].flatten(order='F')
        z = getattr(Float, zvar)[:, idxs].flatten(order='F')
        d = getattr(Float, dvar)[:, idxs].flatten(order='F')
        zs.append(z.copy())

        tgps = getattr(Float, 'UTC_start')[idxs]
        lon = getattr(Float, 'lon_start')[idxs]
        lat = getattr(Float, 'lat_start')[idxs]
        tctd = getattr(Float, 'UTC')[:, idxs].flatten(order='F')
        nans = np.isnan(d) | np.isnan(tctd)
        tctd = tctd[~nans]
        dctd = d[~nans]
        lonctd = np.interp(tctd, tgps, lon)
        latctd = np.interp(tctd, tgps, lat)
        bathy = sandwell.interp_track(lonctd, latctd, bathy_file)

        d -= dctd[bathy.argmax()]
        ds.append(d.copy())

        nans = np.isnan(d) | np.isnan(z) | np.isnan(V)

        Wg = griddata((d[~nans], z[~nans]), V[~nans], (Xg, Zg), method='linear')
        Wgs.append(Wg.copy())

        dctd -= dctd[bathy.argmax()]
        # Spectral analysis of bathymetry.
        dx = np.mean(np.diff(dctd))
        dk = 1./dx


        plt.figure(figsize=(7, 5))
        plt.scatter(d, z, s=50, c=V*100., edgecolor='none', cmap=bwr)
        cbar = plt.colorbar(orientation='horizontal', extend='both')
        cbar.set_label(texvar+' (cm s$^{-1}$)')
        plt.clim(*clim)

        plt.xlim(-10., np.nanmax(d))
        plt.xlabel('Distance (km)')
        plt.ylabel('Depth (m)')
        title_str = ("Float {}").format(Float.floatID)
        plt.title(title_str)

        plt.plot(dctd, bathy, 'k', linewidth=2)

        plt.ylim(np.nanmin(bathy), np.nanmax(z))

        plt.grid()

        my_savefig(Float.floatID, 'mountain_closeup')

# %% Track on bathymetry.

hpids = np.arange(1, 400)


llcrnrlon = -75
llcrnrlat = -60
urcrnrlon = -55
urcrnrlat = -54

lon_lat = np.array([llcrnrlon - 10, urcrnrlon + 10, llcrnrlat - 5, urcrnrlat + 5])

lon_grid, lat_grid, bathy_grid = sandwell.read_grid(lon_lat, bathy_file)
bathy_grid[bathy_grid > 0] = 0

m = bm.Basemap(projection='tmerc', llcrnrlon=llcrnrlon,
               llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon,
               urcrnrlat=urcrnrlat, lon_0=0.5*(llcrnrlon+urcrnrlon),
               lat_0=0.5*(llcrnrlat+urcrnrlat), resolution='f')

plt.figure(figsize=(10, 7))
x, y = m(lon_grid, lat_grid)
m.pcolormesh(x, y, bathy_grid, cmap=plt.get_cmap('binary_r'))
r = np.abs((urcrnrlon-llcrnrlon)/(urcrnrlat-llcrnrlat))

if r > 1.:
    Nm = 8
    Np = max(3, np.round(Nm/r))
    orientation = 'horizontal'
elif r < 1.:
    Np = 8
    Nm = max(3, np.round(Nm/r))
    orientation = 'vertical'

cbar = plt.colorbar(orientation=orientation)
cbar.set_label('Depth (m)')

parallels = np.array([-58, -56])
m.drawparallels(parallels, labels=[1, 0, 0, 0])
meridians = np.array([-65, -55])
m.drawmeridians(meridians, labels=[0, 0, 0, 1])

m.fillcontinents()
m.drawcoastlines()

colours = ['b', 'r']

for Float, c in zip([E76, E77], colours):

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)
    lons = Float.lon_start[idxs]
    lats = Float.lat_start[idxs]

    x, y = m(lons, lats)
    m.plot(x, y, linewidth=4, color=c)

my_savefig(2, 'tracks')

# %% Zoomed in variables near the wave.

hpids = np.arange(10, 50)
bathy_file = '../../data/sandwell_bathymetry/topo_17.1.img'
vars = ['Ww', 'U_abs', 'V_abs']
zvars = ['z', 'zef', 'zef']
dvars = ['dist_ctd', 'dist_ef', 'dist_ef']
texvars = ['$W_w$', '$U$', '$V$']
clims = [(-10., 10.), (-100., 100.), (-100, 100.)]

for Float in [E76, E77]:

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)

    for var, zvar, dvar, texvar, clim in zip(vars, zvars, dvars, texvars,
                                             clims):

        V = getattr(Float, var)[:, idxs].flatten(order='F')
        z = getattr(Float, zvar)[:, idxs].flatten(order='F')
        d = getattr(Float, dvar)[:, idxs].flatten(order='F')

        lons = Float.lon_start[idxs]
        lats = Float.lat_start[idxs]
        dist = Float.dist[idxs]
        bathy = sandwell.interp_track(lons, lats, bathy_file)

        idx = bathy.argmax()

        d -= dist[idx]
        dist -= dist[idx]

        plt.figure(figsize=(7, 5))
        plt.scatter(d, z, s=50, c=V*100., edgecolor='none', cmap=bwr)
        cbar = plt.colorbar(orientation='horizontal', extend='both')
        cbar.set_label(texvar+' (cm s$^{-1}$)')
        plt.clim(*clim)

        plt.xlim(np.nanmin(d), np.nanmax(d))
        plt.xlabel('Distance (km)')
        plt.ylabel('Depth (m)')
        title_str = ("Float {}").format(Float.floatID)
        plt.title(title_str)

        plt.plot(dist, bathy, 'k', linewidth=3)

        plt.ylim(np.nanmin(bathy), np.nanmax(z))

        my_savefig(Float.floatID, var + '_mountain_centered')
