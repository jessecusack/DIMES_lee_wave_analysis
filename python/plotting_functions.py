# -*- coding: utf-8 -*-
"""
Created on Tue Apr 08 12:17:53 2014

@author: jc3e13
"""

import matplotlib.pyplot as plt
import numpy as np
import mpl_toolkits.basemap as bm
import sandwell


def dist_section(Float, hpids, var, plot_func=plt.contourf):
    """ """
    z_vals = np.arange(-1400., -100., 10.)
    __, idxs = Float.get_profiles(hpids, ret_idxs=True)
    dists = Float.dist[idxs]
    __, __, var_grid = Float.get_interp_grid(hpids, z_vals, 'z', var)
    plt.figure()
    plot_func(dists, z_vals, var_grid)


def depth_profile(Float, hpids, var, plot_func=plt.plot):
    """ """
    profiles = Float.get_profiles(hpids)
    if np.iterable(profiles):
        for profile in profiles:
            z = getattr(profile, 'z')
            x = profile.interp(z, 'z', var)
            plt.figure()
            plot_func(x, z)
    else:
        z = getattr(profiles, 'z')
        x = profile.interp(z, 'z', var)
        plt.figure()
        plot_func(x, z)


def isosurface(Float, hpids, var_1_name, var_2_name, var_2_vals):
    """Plot the values of some property at constant surfaces of some other
    property.

    e.g. Depth of potential density surfaces.
    e.g. Temperature of potential density surfaces.
    """

    ig, __, var_1g = Float.get_interp_grid(hpids, var_2_vals,
                                           var_2_name, var_1_name)
    plt.figure()
    plt.plot(ig.T, var_1g.T)
    plt.legend(str(var_2_vals).strip('[]').split())


def timeseries(Float, hpids, var_name):
    """TODO: Docstring..."""

    t, var = Float.get_timeseries(hpids, var_name)
    plt.figure()
    plt.plot(t, var)


def track_on_bathy(Float, hpids, projection='cyl'):
    """TODO: Docstring..."""

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)
    lons = Float.lon_start[idxs]
    lats = Float.lat_start[idxs]

    llcrnrlon = np.floor(np.min(lons)) - 1.
    llcrnrlat = np.floor(np.min(lats)) - 1.
    urcrnrlon = np.ceil(np.max(lons)) + 1.
    urcrnrlat = np.ceil(np.max(lats)) + 1.

    lon_lat = np.array([llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat])

    lon_grid, lat_grid, bathy_grid = sandwell.read_grid(lon_lat)
    bathy_grid[bathy_grid > 0] = 0

    m = bm.Basemap(projection=projection, llcrnrlon=llcrnrlon,
                   llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon,
                   urcrnrlat=urcrnrlat, lon_0=0.5*(llcrnrlon+urcrnrlon),
                   lat_0=0.5*(llcrnrlat+urcrnrlat), resolution='f')

    plt.figure()
    x, y = m(lon_grid, lat_grid)
    m.pcolormesh(x, y, bathy_grid, cmap=plt.get_cmap('binary_r'))
    x, y = m(lons, lats)
    m.plot(x, y, 'r-', linewidth=2)

    m.fillcontinents()
    m.drawcoastlines()

    r = np.abs((urcrnrlon-llcrnrlon)/(urcrnrlat-llcrnrlat))

    if r > 1.:
        Nm = 8
        Np = max(3, np.round(Nm/r))
        orientation = 'horizontal'
    elif r < 1.:
        Np = 8
        Nm = max(3, np.round(Nm/r))
        orientation = 'vertical'

    parallels = np.round(np.linspace(llcrnrlat, urcrnrlat, Np), 1)
    m.drawparallels(parallels, labels=[1, 0, 0, 0])
    meridians = np.round(np.linspace(llcrnrlon, urcrnrlon, Nm), 1)
    m.drawmeridians(meridians, labels=[0, 0, 0, 1])

    cbar = plt.colorbar(orientation=orientation)
    cbar.set_label('Depth (m)')


def bathy_along_track(Float, hpids):
    """TODO: Docstring..."""

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)
    lons = Float.lon_start[idxs]
    lats = Float.lat_start[idxs]
    dist = Float.dist[idxs]
    bathy = sandwell.interp_track(lons, lats)

    plt.figure()
    plt.plot(dist, bathy)
    plt.xlabel('Distance (km)')
    plt.ylabel('Depth (m)')


    