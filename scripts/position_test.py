# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 10:22:36 2014

@author: jc3e13
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumtrapz, trapz
import sys
import os

from geopy.distance import vincenty

lib_path = os.path.abspath('../modules')
if lib_path not in sys.path:
    sys.path.append(lib_path)

import emapex

try:
    print("Floats {} and {} exist!.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.load(4976)
    E77 = emapex.load(4977)

# %% want an odd-even pair
Float = E76
pfls, idxs = Float.get_profiles(np.array([305, 306]), ret_idxs=True)

#plt.figure()
#plt.plot(Float.dist_ctd[:, idxs], Float.z[:, idxs])
#plt.xlabel('Distance (km)')
#plt.ylabel('z (m)')

z = Float.zef[:, idxs]
U = Float.U[:, idxs]
V = Float.V[:, idxs]
Ua = Float.U_abs[:, idxs]
Va = Float.V_abs[:, idxs]
T = Float.UTCef[:, idxs]*86400

zf = z.flatten(order='F')
Uf = U.flatten(order='F')
Vf = V.flatten(order='F')
Uaf = Ua.flatten(order='F')
Vaf = Va.flatten(order='F')
Tf = T.flatten(order='F')
dT = np.nanmax(Tf) - np.nanmin(Tf)

nans = np.isnan(Uf) | np.isnan(Vf) | np.isnan(Tf)

# A half profile pair is bounded by a box defined by the lon-lat positions at
# its corners.
lon1 = Float.lon_start[idxs[0]]
lon2 = Float.lon_end[idxs[1]]
lat1 = Float.lat_start[idxs[0]]
lat2 = Float.lat_end[idxs[1]]

fX = -1 if lon1 > lon2 else 1.
fY = -1 if lat1 > lat2 else 1.

lonl = np.min((lon1, lon2))
lonr = np.max((lon1, lon2))
latb = np.min((lat1, lat2))
latt = np.max((lat1, lat2))

X = fX*vincenty((latb, lonl), (latb, lonr)).m
Y = fY*vincenty((latb, lonl), (latt, lonl)).m
D = vincenty((lat1, lon1), (lat2, lon2)).m
dD = np.sqrt(np.abs(X**2 + Y**2 - D**2))

x1 = cumtrapz(Uf[~nans], Tf[~nans], axis=0, initial=0.)
y1 = cumtrapz(Vf[~nans], Tf[~nans], axis=0, initial=0.)

U_abs = Uf[~nans] + X/dT - trapz(Uf[~nans], Tf[~nans], axis=0)/dT
V_abs = Vf[~nans] + Y/dT - trapz(Vf[~nans], Tf[~nans], axis=0)/dT

x2 = cumtrapz(U_abs, Tf[~nans], axis=0, initial=0.)
y2 = cumtrapz(V_abs, Tf[~nans], axis=0, initial=0.)

x3 = cumtrapz(Uaf[~nans], Tf[~nans], axis=0, initial=0.)
y3 = cumtrapz(Vaf[~nans], Tf[~nans], axis=0, initial=0.)

fig, axs = plt.subplots(1, 2, sharey=True)
axs[0].plot(U, z)
axs[0].set_ylabel('z (m)')
axs[0].set_xlabel('U (m/s)')
axs[1].plot(V, z)
axs[1].set_xlabel('V (m/s)')

axs[0].plot(U_abs, zf[~nans], 'k--')
axs[1].plot(V_abs, zf[~nans], 'k--')

plt.figure()
#plt.plot(lons, lats, 'k--')
plt.plot(x1, y1, 'y')
plt.plot(x2, y2, 'r')
plt.plot(x3, y3, 'k--')
plt.plot([X, X], [0., Y], 'k')
plt.plot([0., X], [0., 0.], 'k')
plt.plot([0., X], [0., Y], 'k')
ax = plt.gca()
ax.ticklabel_format(useOffset=False)

# %%

## Fudge factors.
#fx = (Float.lon_end[idxs[1]] - Float.lon_start[idxs[0]])/(lons[-1] - lons[0])
#fy = (Float.lat_end[idxs[1]] - Float.lat_start[idxs[0]])/(lats[-1] - lats[0])
#
#x *= fx
#y *= fy
#
#lons, lats = utils.distll(Float.lon_start[idxs[0]], Float.lat_start[idxs[0]], x, y)
#
#lon_ef = np.empty_like(T)
#lat_ef = np.empty_like(T)
#
#for i in range(2):
#    lon_ef[:, i] = utils.nan_interp(T[:, i], Tf[~nans], lons)
#    lat_ef[:, i] = utils.nan_interp(T[:, i], Tf[~nans], lats)
#
#plt.plot(lon_ef, lat_ef)
#plt.plot(Float.lon_start[idxs], Float.lat_start[idxs], 'yo', Float.lon_gps[[idxs[0]-1, idxs[-1]]],
#         Float.lat_gps[[idxs[0]-1, idxs[-1]]], 'r^', Float.lon_end[idxs], Float.lat_end[idxs], 'y*')
#
#plt.title("Scaling x: {:1.2f}, y: {:1.2f}".format(fx, fy))
#plt.xlabel('Longitude')
#plt.ylabel('Latitude')
