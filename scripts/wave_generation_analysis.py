# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 16:24:05 2014

@author: jc3e13
"""

# %% Loads and imports.

import numpy as np
#import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits import basemap as bm
import gsw
import os
import sys
import glob

lib_path = os.path.abspath('../python')
if lib_path not in sys.path:
    sys.path.append(lib_path)

import sandwell
import utils
import emapex
import plotting_functions as pf

try:
    print("Floats {} and {} exist!".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4976)
    E76.apply_w_model('../../data/EM-APEX/4976_fix_p0k0M_fit_info.p')
    E76.apply_strain('../../data/EM-APEX/4976_N2_ref_300dbar.p')
    E76.apply_isopycnal_displacement('../../data/EM-APEX/srho_4976_100mbin.p')
    E77 = emapex.EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4977)
    E77.apply_w_model('../../data/EM-APEX/4977_fix_p0k0M_fit_info.p')
    E77.apply_strain('../../data/EM-APEX/4977_N2_ref_300dbar.p')
    E77.apply_isopycnal_displacement('../../data/EM-APEX/srho_4977_100mbin.p')


# %% Script params.

# Bathymetry file path.
bf = os.path.abspath(glob.glob('../../data/sandwell_bathymetry/topo_*.img')[0])
# Figure save path.
sdir = '../figures/wave_generation_analysis'
if not os.path.exists(sdir):
    os.makedirs(sdir)
# Universal figure font size.
matplotlib.rc('font', **{'size': 9})

################
# START SCRIPT #
################

# %% Float trajectories FULL
# ------------------

lons = np.hstack([Float.lon_start for Float in [E76, E77]])
lats = np.hstack([Float.lat_start for Float in [E76, E77]])

llcrnrlon = np.floor(np.min(lons)) - 1.
llcrnrlat = np.floor(np.min(lats)) - 1.
urcrnrlon = np.ceil(np.max(lons)) + 1.
urcrnrlat = np.ceil(np.max(lats)) + 1.

lon_lat = np.array([llcrnrlon, urcrnrlon+5, llcrnrlat-5, urcrnrlat])

lon_grid, lat_grid, bathy_grid = sandwell.read_grid(lon_lat, bf)
bathy_grid[bathy_grid > 0] = 0

m = bm.Basemap(projection='tmerc', llcrnrlon=llcrnrlon,
               llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon,
               urcrnrlat=urcrnrlat, lon_0=0.5*(llcrnrlon+urcrnrlon),
               lat_0=0.5*(llcrnrlat+urcrnrlat), resolution='f')

fig = plt.figure()
x, y = m(lon_grid, lat_grid)
m.pcolormesh(x, y, bathy_grid, cmap=plt.get_cmap('binary_r'))

for Float, colour in zip([E76, E77], ['b', 'g']):
    x, y = m(Float.lon_start, Float.lat_start)
    m.plot(x, y, colour, linewidth=3, label=Float.floatID)

plt.legend()

m.fillcontinents()
m.drawcoastlines()

r = np.abs((urcrnrlon-llcrnrlon)/(urcrnrlat-llcrnrlat))

if r > 1.:
    Nm = 8
    Np = max(4, np.round(Nm/r))
    orientation = 'horizontal'
elif r < 1.:
    Np = 8
    Nm = max(4, np.round(Nm/r))
    orientation = 'vertical'
else:
    Np = 4
    Nm = 4
    orientation = 'horizontal'

parallels = np.round(np.linspace(llcrnrlat, urcrnrlat, Np), 1)
m.drawparallels(parallels, labels=[1, 0, 0, 0])
meridians = np.round(np.linspace(llcrnrlon, urcrnrlon, Nm), 1)
m.drawmeridians(meridians, labels=[0, 0, 0, 1])

cbar = plt.colorbar(orientation=orientation)
cbar.set_label('Depth (m)')
pf.my_savefig(fig, 'both', 'traj', sdir, fsize='double_col')

# %% Float quivers RIDGE ONLY
# ------------------

# Number of half profiles.
N = 30
N0_76 = 15
N0_77 = 10
hpids = np.empty((N, 2))
hpids[:, 0] = np.arange(N0_76, N0_76+N)
hpids[:, 1] = np.arange(N0_77, N0_77+N)

lons = np.empty_like(hpids)
lats = np.empty_like(hpids)
Us = np.empty_like(hpids)
Vs = np.empty_like(hpids)
Ws = np.empty_like(hpids)

for i, Float in enumerate([E76, E77]):
    __, idxs = Float.get_profiles(hpids[:, i], ret_idxs=True)
    lons[:, i] = Float.lon_gps[idxs]
    lats[:, i] = Float.lat_gps[idxs]

    zef = Float.zef[:, idxs]
    z = Float.z[:, idxs]
    U = Float.U_abs[:, idxs]
    V = Float.V_abs[:, idxs]
    W = Float.Ww[:, idxs]

    for j, (zefpfl, zpfl, Upfl, Vpfl, Wpfl) in enumerate(zip(zef.T, z.T, U.T,
                                                             V.T, W.T)):
        Us[j, i] = np.nanmean(Upfl[(-1400 < zefpfl) & (zefpfl < -1000)])
        Vs[j, i] = np.nanmean(Vpfl[(-1400 < zefpfl) & (zefpfl < -1000)])
        Ws[j, i] = np.nanmax(Wpfl[zpfl < -100])

llcrnrlon = np.floor(np.nanmin(lons)) - .2
llcrnrlat = np.floor(np.nanmin(lats)) - .2
urcrnrlon = np.ceil(np.nanmax(lons)) + .2
urcrnrlat = np.ceil(np.nanmax(lats)) + .2

lon_lat = np.array([llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat])

lon_grid, lat_grid, bathy_grid = sandwell.read_grid(lon_lat, bf)
bathy_grid[bathy_grid > 0] = 0

m = bm.Basemap(projection='tmerc', llcrnrlon=llcrnrlon,
               llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon,
               urcrnrlat=urcrnrlat, lon_0=0.5*(llcrnrlon+urcrnrlon),
               lat_0=0.5*(llcrnrlat+urcrnrlat), resolution='f')

fig = plt.figure(figsize=(10, 10))
x, y = m(lon_grid, lat_grid)
levels = np.arange(-4000, 0, 500)
CS = m.contour(x, y, bathy_grid, 20, cmap=plt.get_cmap('binary'),
               levels=levels)
plt.clabel(CS, inline=1, fontsize=8, fmt='%1.f')
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

marker = ['o', '*']
label = ['4976', '4977']
color = ['b', 'g']
for i, (lon, lat, U, V, W) in enumerate(zip(lons.T, lats.T, Us.T, Vs.T, Ws.T)):
    x, y = m(lon, lat)
    m.plot(x, y, marker[i], color='y', label=label[i])
    Q = m.quiver(x, y, U, V, W, scale=6, cmap=plt.get_cmap('jet'))
    plt.clim(0, 0.3)
#    for _x, _y, hpid in zip(x, y, hpids[:, i]):
#        plt.annotate("{:1.0f}".format(hpid), xy=(_x, _y),
#                     xytext=(1.03*_x, 1.03*_y), color=color[i])

plt.legend()
qk = plt.quiverkey(Q, 0.5, 0.92, 0.5, r'0.5 m s$^{-1}$', labelpos='N')
cbar = plt.colorbar(orientation='horizontal')
cbar.set_label('Maximum $W$ (m s$^{-1}$)')

pf.my_savefig(fig, 'both', 'quiver_traj', sdir, fsize='double_col')

# %% Upstream flow properties
# ------------------------

N = 10
N0_76 = 15
N0_77 = 10
E76_hpids = np.arange(N0_76, N0_76+N)
E77_hpids = np.arange(N0_77, N0_77+N)

pfls = np.hstack((E76.get_profiles(E76_hpids), E77.get_profiles(E77_hpids)))

fig, axs = plt.subplots(1, 7, sharey=True)
fig.set_size_inches(16, 5)

axs[0].set_ylabel('$z$ (m)')

for pfl in pfls:

    axs[0].plot(pfl.T, pfl.z, color='k', alpha=0.3)
    axs[1].plot(pfl.S, pfl.z, color='k', alpha=0.3)
    axs[2].plot(pfl.rho_1-1000., pfl.z, color='k', alpha=0.3)
    axs[3].plot(np.sqrt(pfl.N2_ref)*1000., pfl.z, color='k', alpha=0.3)
    axs[4].plot(pfl.U_abs*100., pfl.zef, color='k', alpha=0.3)
    axs[5].plot(pfl.V_abs*100., pfl.zef, color='k', alpha=0.3)
    axs[6].plot(pfl.Ww*100., pfl.z, color='k', alpha=0.3)

xlabels = ['$T$ ($^\circ$C)', '$S$ (-)', '$\sigma_1$ (kg m$^{-3}$)',
           '$N$ (10$^{-3}$ rad s$^{-1}$)', '$u$ (cm s$^{-1}$)',
           '$v$ (cm s$^{-1}$)', '$w$ (cm s$^{-1}$)']

#axs[2].ticklabel_format(useOffset=1000.)

for ax, xlabel in zip(axs, xlabels):
    ax.set_xticks(ax.get_xticks()[::2])
    ax.set_xlabel(xlabel)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation='vertical')
# TODO: Plot mean profiles too!

# Why is this necessary...? I don't know but it has to be done.
axs[0].ticklabel_format(useOffset=False)
plt.tight_layout()

pf.my_savefig(fig, 'both', 'upstream', sdir, fsize='double_col')

# The average flow properties below
z_max = -600.
N_mean = np.nanmean(np.hstack([np.sqrt(pfl.N2_ref)[pfl.z < z_max]
                               for pfl in pfls]))
U_mean = np.nanmean(np.hstack([pfl.U_abs[pfl.zef < z_max] for pfl in pfls]))

print("Mean buoyancy frequency below {:1.0f} m is {:.2E} rad s-1. "
      "The period is {:1.1f} min.".format(z_max, N_mean, 2*np.pi/(60*N_mean)))
print("Mean zonal speed below {:1.0f} m is {:1.2f} m s-1.".format(z_max,
      U_mean))

# %% Topography
# ----------

# Number of half profiles.
N = 15
N0_76 = 25
N0_77 = 20
hpids = np.empty((N, 2))
hpids[:, 0] = np.arange(N0_76, N0_76+N)
hpids[:, 1] = np.arange(N0_77, N0_77+N)

lons = np.empty_like(hpids)
lats = np.empty_like(hpids)

for i, Float in enumerate([E76, E77]):
    __, idxs = Float.get_profiles(hpids[:, i], ret_idxs=True)
    lons[:, i] = Float.lon_start[idxs]
    lats[:, i] = Float.lat_start[idxs]

lons = lons.flatten(order='F')
lats = lats.flatten(order='F')
p = np.polyfit(lons, lats, 1)

llcrnrlon = np.nanmin(lons) - .02
llcrnrlat = np.nanmin(lats) - .02
urcrnrlon = np.nanmax(lons) + .02
urcrnrlat = np.nanmax(lats) + .02

lon_lat = np.array([llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat])

lon_grid, lat_grid, bathy_grid = sandwell.read_grid(lon_lat, bf)
bathy_grid[bathy_grid > 0] = 0

fig = plt.figure()
plt.plot(lons, lats, 'o')
plt.plot(lons, np.polyval(p, lons))
levels = np.arange(-4000, 0, 500)
CS = plt.contour(lon_grid, lat_grid, bathy_grid, 20,
                 cmap=plt.get_cmap('binary'), levels=levels)
plt.clabel(CS, inline=1, fontsize=8, fmt='%1.f')
plt.gca().ticklabel_format(useOffset=False)
plt.xlabel('Longitude')
plt.ylabel('Latitude')

llons = np.linspace(lons.min(), lons.max(), 50)
llats = np.polyval(p, llons)
dist = np.hstack((0, np.cumsum(utils.lldist(llons, llats))))
topo = sandwell.interp_track(llons, llats, bf)

pf.my_savefig(fig, 'both', 'line_over_topo', sdir, fsize='double_col')

fig = plt.figure()

plt.plot(dist, topo, 'r', linewidth=5)
plt.xlabel('Distance (km)')
plt.ylabel('$z$ (m)')
plt.grid()

z = np.linspace(*plt.ylim(), num=10)
xg, zg = np.meshgrid(dist, z)
tg = zg - np.max(topo)
levels = np.arange(-1000, 1, 100)
CS = plt.contour(xg, zg, tg, colors='k', levels=levels)
plt.clabel(CS, inline=1, fontsize=8, fmt='%1.f')

pf.my_savefig(fig, 'topo', 'section', sdir, fsize='double_col')

# %% Characteristics of the flow
# ---------------------------

# The inertial frequency...
lat_mean = np.mean(lats)
f = np.abs(gsw.f(lat_mean))
print("Mean inertial frequency at {:1.1f} degrees is {:.2E} rad s-1. "
      "The period is {:1.1f} hours.".format(lat_mean, f, 2*np.pi/(60*60*f)))

# - $N \approx 2.3 \times 10^{-3}$ rad s$^{-1}$. (Near bottom of profiles
# upstream of the ridge.)
# - $U \approx 0.3$ m s$^{-1}$. (Near bottom of profiles upstream of the
# ridge.)
# - $f \approx 1.2 \times 10^{-4}$ rad s$^{-1}$. (57 degrees south.)
# - $h_0 \approx 2 \times 10^3$ m (Ridge height.)
# - $L \approx 2 \times 10^4$ m (Ridge width.)
#
# Periods of motion.
#
# - Buoyancy period, $T_N \approx 50$ mins.
# - Inertial period, $T_f \approx 14$ hours.

N0 = 0.8e-3
N1 = 1.0e-3
U0 = 0.2
U1 = 0.4
f = 1.2e-4
h0 = 300.
h1 = 1.5e3
L0 = 2e3
L1 = 1e4

# Assuming waves generated have a frequency near N, the horizontal scale is...
L_N0 = 2*np.pi*U0/N1
L_N1 = 2*np.pi*U1/N0
print("Assuming that generated waves have a frequency near N then their "
      "horizontal wavelength should be of order {:1.0f} -- {:1.0f} m.".format(L_N0, L_N1))
# Not assuming this the appropriate frequency is...
om_L0 = 2*np.pi*U0/L1
om_L1 = 2*np.pi*U1/L0
print("Topographic and velocity scales suggest that their frequency is in the "
      "range {:.2E} -- {:.2E} rad s-1.".format(om_L0, om_L1))
Ro0 = U1/(f*L1)
Ro1 = U1/(f*L0)
print("The Rossby number for the flow given a range of length scales is "
      "between {} -- {}.".format(Ro0, Ro1))
Fr0 = U0/(N1*h1)
Fr1 = U1/(N0*h0)
print("The Froude number given a range of height scales is  between {:1.2f} "
      "-- {:1.2f}.".format(Fr0, Fr1))
print("Consequenctly the steepness parameter (1/Fr) is in the range {:1.1f} "
      "-- {:1.1f}.".format(1./Fr1, 1./Fr0))

hc0 = 2*np.pi*0.4*U0/N1
hc1 = 2*np.pi*0.4*U1/N0
print("The height of the obstacle that would have a steepness of 0.4 is "
      "either {:1.0f} m or {:1.0f} m.\n"
      "This depends on a factor of 2 pi and possibly another factor of 2 due "
      "to uncertainty in N.".format(hc0, hc1))


# %% Mean rho_1 profile

N = 10
N0_76 = 15
N0_77 = 10
E76_hpids = np.arange(N0_76, N0_76+N)
E77_hpids = np.arange(N0_77, N0_77+N)

dz = 1.
z = np.arange(-1500, 0, dz)
rho = []

pfls = np.hstack((E76.get_profiles(E76_hpids), E77.get_profiles(E77_hpids)))

fig, axs = plt.subplots(1, 2)

axs[0].set_ylabel('$z$ (m)')

for pfl in pfls:
    axs[0].plot(pfl.rho_1, pfl.z, color='grey')
    axs[0].set_xlabel('$\sigma_1$ (kg m$^{-3}$)')
    plt.setp(axs[0].xaxis.get_majorticklabels(), rotation=45)

    axs[1].plot(pfl.srho_1, pfl.z, color='grey')
    axs[1].set_xlabel('$\sigma_1$ (kg m$^{-3}$)')
    plt.setp(axs[1].xaxis.get_majorticklabels(), rotation=45)

    rho.append(pfl.interp(z, 'z', 'rho_1'))

rho = np.transpose(np.asarray(rho))
mrho = np.mean(rho, axis=-1)

axs[0].plot(mrho, z, 'red')
axs[1].plot(mrho, z, 'red')

pfl = E76.get_profiles(32)
srhop = utils.nan_interp(pfl.z, z, mrho)
b = -gsw.grav(pfl.lat_start, pfl.P)*(pfl.rho_1 - srhop)/1031.

fig, axs = plt.subplots(1, 2)
axs[0].plot(pfl.b, pfl.z, 'grey')
axs[0].plot(b, pfl.z, 'red')
axs[0].plot(utils.nan_detrend(pfl.z, b), pfl.z, 'red', linestyle='--')
axs[1].plot(pfl.rho_1, pfl.z, 'grey')
axs[1].plot(mrho, z, 'red')
axs[1].plot(pfl.srho_1, pfl.z, 'grey', linestyle='--')
