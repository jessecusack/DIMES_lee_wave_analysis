# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 10:22:36 2014

@author: jc3e13
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumtrapz
import sys
import os

lib_path = os.path.abspath('../modules')
if lib_path not in sys.path:
    sys.path.append(lib_path)

import utils
import emapex

try:
    print("Floats {} and {} exist!.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.load(4976)
    E77 = emapex.load(4977)

# want an odd-even pair
Float = E76
pfls, idxs = Float.get_profiles(np.array([115, 116]), ret_idxs=True)

plt.figure()
plt.plot(Float.dist_ctd[:, idxs], Float.z[:, idxs])
plt.xlabel('Distance (km)')
plt.ylabel('z (m)')

z = Float.zef[:,idxs]
U = Float.U_abs[:, idxs] - Float.U[:, idxs]# - 0.1
V = Float.V_abs[:, idxs] - Float.V[:, idxs]# + 0.03
T = Float.UTCef[:, idxs]*86400

fig, axs = plt.subplots(1, 2, sharey=True)
axs[0].plot(U, z)
axs[0].set_ylabel('z (m)')
axs[0].set_xlabel('U (m/s)')
axs[1].plot(V, z)
axs[1].set_xlabel('V (m/s)')

Uf = U.flatten(order='F')
Vf = V.flatten(order='F')
Tf = T.flatten(order='F')

nans = np.isnan(Uf) | np.isnan(Vf) | np.isnan(Tf)

#D = np.nanmax(Float.dist_ef[:,idxs], axis=0) - np.nanmin(Float.dist_ef[:,idxs], axis=0)

x = cumtrapz(Uf[~nans], Tf[~nans], axis=0, initial=0.)/1e3
y = cumtrapz(Vf[~nans], Tf[~nans], axis=0, initial=0.)/1e3

lons, lats = utils.distll(Float.lon_start[idxs[0]], Float.lat_start[idxs[0]], x, y)

plt.figure()
plt.plot(lons, lats, 'k--')
ax = plt.gca()
ax.ticklabel_format(useOffset=False)

# Fudge factors.
fx = (Float.lon_end[idxs[1]] - Float.lon_start[idxs[0]])/(lons[-1] - lons[0])
fy = (Float.lat_end[idxs[1]] - Float.lat_start[idxs[0]])/(lats[-1] - lats[0])

x *= fx
y *= fy

lons, lats = utils.distll(Float.lon_start[idxs[0]], Float.lat_start[idxs[0]], x, y)

lon_ef = np.empty_like(T)
lat_ef = np.empty_like(T)

for i in range(2):
    lon_ef[:, i] = utils.nan_interp(T[:, i], Tf[~nans], lons)
    lat_ef[:, i] = utils.nan_interp(T[:, i], Tf[~nans], lats)

plt.plot(lon_ef, lat_ef)
plt.plot(Float.lon_start[idxs], Float.lat_start[idxs], 'yo', Float.lon_gps[[idxs[0]-1, idxs[-1]]],
         Float.lat_gps[[idxs[0]-1, idxs[-1]]], 'r^', Float.lon_end[idxs], Float.lat_end[idxs], 'y*')

plt.title("Scaling x: {:1.2f}, y: {:1.2f}".format(fx, fy))
plt.xlabel('Longitude')
plt.ylabel('Latitude')
