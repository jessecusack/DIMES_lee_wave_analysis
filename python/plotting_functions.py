# -*- coding: utf-8 -*-
"""
Created on Tue Apr 08 12:17:53 2014

@author: jc3e13
"""

import matplotlib.pyplot as plt
import numpy as np
import mpl_toolkits.basemap as bm
import sandwell
import mat2py as m2p
import scipy.signal as sig


def dist_section(Float, hpids, var, plot_func=plt.contourf):
    """ """
    z_vals = np.arange(-1400., -100., 10.)
    __, idxs = Float.get_profiles(hpids, ret_idxs=True)
    dists = Float.dist[idxs]
    __, __, var_grid = Float.get_interp_grid(hpids, z_vals, 'z', var)
    plt.figure()
    plot_func(dists, z_vals, var_grid)


def scatter_section(Float, hpids, var, x_var='dist'):

    z_vals = np.arange(-1500., -40., 20.)
    __, idxs = Float.get_profiles(hpids, ret_idxs=True)
    __, z, var = Float.get_interp_grid(hpids, z_vals, 'z', var)

    if x_var == 'dist':

        __, __, d = Float.get_interp_grid(hpids, z_vals, 'z', 'dist_ctd_data')
        d = d.flatten(order='F')

    elif x_var == 'time':

        __, __, d = Float.get_interp_grid(hpids, z_vals, 'z', 'UTC')
        d = d.flatten(order='F')
        d = m2p.datenum_to_datetime(d)

    else:
        raise ValueError("Input x_var should be 'dist' or 'time'.")

    z = z.flatten(order='F')
    var = var.flatten(order='F')

    plt.figure()
    plt.scatter(d, z, c=var, edgecolor='none')
    plt.ylim(np.min(z), np.max(z))
    plt.xlim(np.min(d), np.max(d))
    plt.colorbar(orientation='horizontal', extend='both')


def depth_profile(Float, hpids, var, plot_func=plt.plot, hold='off'):
    """ """
    profiles = Float.get_profiles(hpids)
    if np.iterable(profiles):
        for i, profile in enumerate(profiles):
            z = getattr(profile, 'z')
            x = profile.interp(z, 'z', var)
            if hold == 'on':
                if i == 0:
                    plt.figure()
                plot_func(x, z)
            else:
                plt.figure()
                plot_func(x, z)

    else:
        z = getattr(profiles, 'z')
        x = profiles.interp(z, 'z', var)
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


# Source of hinton: http://wiki.scipy.org/Cookbook/Matplotlib/HintonDiagrams
def _blob(x, y, area, colour):
    """
    Draws a square-shaped blob with the given area (< 1) at
    the given coordinates.
    """
    hs = np.sqrt(area) / 2.
    xcorners = np.array([x - hs, x + hs, x + hs, x - hs])
    ycorners = np.array([y - hs, y - hs, y + hs, y + hs])
    plt.fill(xcorners, ycorners, colour, edgecolor=colour)


def hinton(W, maxWeight=None):
    """
    Draws a Hinton diagram for visualizing a weight matrix.
    Temporarily disables matplotlib interactive mode if it is on,
    otherwise this takes forever.
    """
    reenable = False
    if plt.isinteractive():
        plt.ioff()
    plt.clf()
    height, width = W.shape
    if not maxWeight:
        maxWeight = 2**np.ceil(np.log(np.max(np.abs(W)))/np.log(2))

    plt.fill(np.array([0, width, width, 0]), np.array([0, 0, height, height]),
             'gray')
    plt.axis('off')
    plt.axis('equal')
    for x in xrange(width):
        for y in xrange(height):
            _x = x+1
            _y = y+1
            w = W[y, x]
            if w > 0:
                _blob(_x - 0.5, height - _y + 0.5,
                      min(1, w/maxWeight), 'white')
            elif w < 0:
                _blob(_x - 0.5, height - _y + 0.5,
                      min(1, -w/maxWeight), 'black')
    if reenable:
        plt.ion()
    plt.show()


def welch_psd(Float, hpids, var, tz='z', hold='off'):
    """Compute power spectral density of some variable using Welch method.
    Variables are first interpolated onto a regular grid in either time or
    depth which can be specified using the tz optional argument. A time
    interval of 30 seconds is used and a depth interval of 6m."""

    dz = 6.
    dt = 30./86400.

    if tz == 'z':
        df = dz
        ivar = 'z'
        m = 1.
    elif tz == 't':
        df = dt
        ivar = 'UTC'
        m = 86400.
    else:
        raise RuntimeWarning("tz should be 't' or 'z'.")

    profiles = Float.get_profiles(hpids)
    if np.iterable(profiles):
        for i, profile in enumerate(profiles):
            f = getattr(profile, ivar)
            nans = np.isnan(f)
            f = np.unique(f[~nans])
            f = np.arange(f[0], f[-1], df)

            x = profile.interp(f, ivar, var)

            if hold == 'on':
                if i == 0:
                    plt.figure()
                freq, Pxx = sig.welch(x, 1./(m*df))
                plt.loglog(freq, Pxx)
            else:
                plt.figure()
                freq, Pxx = sig.welch(x, 1./(m*df))
                plt.loglog(freq, Pxx)

    else:
        f = getattr(profiles, ivar)
        nans = np.isnan(f)
        f = np.unique(f[~nans])

        f = np.arange(f[0], f[-1], df)

        x = profiles.interp(f, ivar, var)
        plt.figure()
        freq, Pxx = sig.welch(x, 1./(m*df))
        plt.loglog(freq, Pxx)
