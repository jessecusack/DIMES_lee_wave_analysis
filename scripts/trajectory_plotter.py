# -*- coding: utf-8 -*-
"""
Created on Tue Oct  7 11:34:21 2014

@author: jc3e13
"""


import numpy as np
from glob import glob
from scipy import io
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import os
import sys


lib_path = os.path.abspath('../modules')
if lib_path not in sys.path:
    sys.path.append(lib_path)

import utils
import emapex
import sandwell


# Fine all float IDs
files = glob('/noc/users/jc3e13/storage/DIMES/EM-APEX/allprof*.mat')
FIDs = np.array([])  # Stores all float IDs.
in_file = []
var_keys = ['flid', 'lon_gps', 'lat_gps', 'utc_dep']

for f in files:

    data = io.loadmat(f, squeeze_me=True, variable_names=var_keys)
    fids = data['flid']
    ufids = np.unique(fids[~np.isnan(fids)])
    FIDs = np.hstack((FIDs, ufids))

    for fid in ufids:

        in_file.append(f)

# Load all floats.
Floats = [emapex.EMApexFloat(f, FID, False, False)
          for f, FID in zip(in_file, FIDs)]

# Generate map.
dt = 24.
fstr = 'all'

t_mins, t_maxs = [], []
for Float in Floats:
    t_mins.append(np.min(Float.UTC_start))
    t_maxs.append(np.max(Float.UTC_start))

t_min = np.min(t_mins)
t_max = np.max(t_maxs)
ti = np.arange(t_min, t_max, dt/24.)
lons = np.empty((ti.size, len(Floats)))
lats = np.empty((ti.size, len(Floats)))

for i, Float in enumerate(Floats):
    lon = Float.lon_start
    lat = Float.lat_start
    t = Float.UTC_start
    lons[:, i] = np.interp(ti, t, lon, left=np.nan, right=np.nan)
    lats[:, i] = np.interp(ti, t, lat, left=np.nan, right=np.nan)

llcrnrlon = np.floor(np.nanmin(lons)) - 1.
llcrnrlat = np.floor(np.nanmin(lats)) - 1.
urcrnrlon = np.ceil(np.nanmax(lons)) + 1.
urcrnrlat = np.ceil(np.nanmax(lats)) + 1.

lon_lat = np.array([llcrnrlon, -30., -70, urcrnrlat])

lon_grid, lat_grid, bathy_grid = sandwell.read_grid(lon_lat)
bathy_grid[bathy_grid > 0] = 0

m = bm.Basemap(projection='cyl', llcrnrlon=llcrnrlon,
               llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon,
               urcrnrlat=urcrnrlat, lon_0=0.5*(llcrnrlon+urcrnrlon),
               lat_0=0.5*(llcrnrlat+urcrnrlat), resolution='l')

plt.figure(figsize=(10, 10))
x, y = m(lon_grid, lat_grid)
m.contour(x, y, bathy_grid, cmap=plt.get_cmap('binary_r'),
          levels=[-4000., -2000.])
m.fillcontinents()
m.drawcoastlines()

r = np.abs((urcrnrlon-llcrnrlon)/(urcrnrlat-llcrnrlat))

if r > 1.:
    Nm = 8
    Np = max(3, np.round(Nm/r))

elif r < 1.:
    Np = 8
    Nm = max(3, np.round(Nm/r))

parallels = np.round(np.linspace(llcrnrlat, urcrnrlat, Np), 1)
m.drawparallels(parallels, labels=[1, 0, 0, 0])
meridians = np.round(np.linspace(llcrnrlon, urcrnrlon, Nm), 1)
m.drawmeridians(meridians, labels=[0, 0, 0, 1])

for i, (time, lonr, latr) in enumerate(zip(ti, lons, lats)):

    print('Creating... {i:03d}'.format(i=i))
    plt.title(utils.datenum_to_datetime(time).strftime('%Y-%m-%d %H:%M'))

    for j, (lon, lat) in enumerate(zip(lonr, latr)):
        x, y = m(lon, lat)
        m.plot(x, y, '.', color='r')

    save_name = '../figures/animated_tracks/{}{i:03d}.png'.format(fstr,
                                                                  i=i)
    plt.savefig(save_name, bbox_inches='tight')

print('Finished.')
