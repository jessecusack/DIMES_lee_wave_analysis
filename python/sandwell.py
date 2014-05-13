# -*- coding: utf-8 -*-
"""
Created on Sun Feb 02 17:45:16 2014

Functions for reading Smith and Sandwell bathymetry.

@author: jc3e13
"""

import numpy as np
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm


def bilinear_interpolation(xg, yg, fg, x, y):
    """TODO:"""
    xa = xg[:, 0]
    ya = yg[0, :]

    i1 = np.searchsorted(xa, x)
    i2 = i1 + 1
    j1 = np.searchsorted(ya, y)
    j2 = j1 + 1

    dx = xa[i2] - xa[i1]
    dy = ya[j2] - ya[j1]

    f11, f21, f12, f22 = fg[i1, j1], fg[i2, j1], fg[i1, j2], fg[i2, j2]

    x1, y1, x2, y2 = xa[i1], ya[j1], xa[i2], ya[j2]

    return (f11*(x2 - x)*(y2 - y) + f21*(x - x1)*(y2 - y) +
            f12*(x2 - x)*(y - y1) + f22*(x - x1)*(y - y1))/(dx*dy)


def read_grid(lon_lat, file_path=None):
    """Input format [lonmin lonmax lat_min lat_max].
                     west   east   south   north

    Output:
      bathy_grid, lon_grid, lat_grid

    TODO: More input checks. Improve docstring...
    """

    # Important parameters.
    nlon = 21600         # Number of longitude points.
    nlat = 17280         # Number of latitude points.
    lat_min = -80.738    # Most southern extent of grid.
    lat_max = 80.738     # Most northern extent of grid.
    arcmin = 1./60.      # A single arcminute. (1/60 of a degree)
    rad = np.pi/180.     # A single radian.
    bytes_per_val = 2    # Number of bytes for each datum in file.
    dtype = '>i2'        # Data are big-endian (>) 2 byte signed integers.
    cross_0 = False      # Flag if grid crosses Greenwich meridian.

    west, east, south, north = lon_lat

    if (west < -180.) | (west > 180.) | (east < -180.) | (east > 180.):
        raise ValueError('Longitude out of bounds (-180 to 180).')

    if ((south < lat_min) | (south > lat_max) |
        (north < lat_min) | (north > lat_max)):

        raise ValueError(
            'Latitude out of bounds ({} to {}).'.format(lat_min, lat_max)
            )

    if north < south:
        raise ValueError('The north latitude is less than the south.')

    south = 90. - south
    north = 90. - north

    if west < 0:
        west = west + 360.
    if east < 0:
        east = east + 360.

    # Mercator projection transformations. (Ref: Wikipedia)
    y = lambda phi: np.log(np.tan(np.pi/4. + phi/2.))
    phi = lambda y: 2*np.arctan(np.exp(y)) - np.pi/2.

    all_lons = np.arange(0., 360., arcmin)
    all_lats = 90. - phi(np.linspace(y(lat_max*rad), y(lat_min*rad), nlat))/rad

    loni1, loni2 = all_lons.searchsorted([west, east]) + 1
    lati1, lati2 = all_lats.searchsorted([north, south]) + 1
    Nlats = lati2 - lati1
    lats = all_lats[lati1:lati2]

    if east < west:
        cross_0 = True
        lons = np.concatenate((all_lons[(loni1 - nlon):], all_lons[:loni2]))
    else:
        Nlons = loni2 - loni1
        lons = all_lons[loni1:loni2]

    lat_grid, lon_grid = np.meshgrid(lats, lons)
    bathy_grid = np.ndarray(lat_grid.shape, dtype='i2')

    if file_path is None:
        file_path = '../../data/sandwell_bathymetry/topo_16.1.img'

    with open(file_path, 'rb') as f:
        for i in xrange(Nlats):
            if cross_0:
                f.seek(bytes_per_val*((lati1 + i)*nlon + loni1))
                N = nlon - loni1
                bathy_grid[:N, i] = np.fromfile(f, dtype=dtype, count=N)

                f.seek(bytes_per_val*(lati1 + i)*nlon)
                bathy_grid[N:, i] = np.fromfile(f, dtype=dtype, count=loni2)
            else:
                f.seek(bytes_per_val*((lati1 + i)*nlon + loni1))
                bathy_grid[:, i] = np.fromfile(f, dtype=dtype, count=Nlons)

    lat_grid = 90 - lat_grid
    lon_grid[lon_grid > 180.] = lon_grid[lon_grid > 180.] - 360.

    return lon_grid, lat_grid, bathy_grid


def interp_track(lons, lats, file_path=None, kind='linear'):
    """Interpolates bathymetry data to given longitude and latitude
    coordinates. The inputs should be one dimensional arrays."""

    margin = 1
    lon_lat = np.array([np.min(lons)-margin, np.max(lons)+margin,
                        np.min(lats)-margin, np.max(lats)+margin])
    lon_grid, lat_grid, bathy_grid = read_grid(lon_lat, file_path)
    return bilinear_interpolation(lon_grid, lat_grid, bathy_grid, lons, lats)


if __name__ == '__main__':

    # Plotting bathymetry example.
    lon_lat = np.array([-72, -27, -68, -47])
    west, east, south, north = lon_lat

    lon_grid, lat_grid, bathy_grid = read_grid(lon_lat)
    bathy_grid[bathy_grid > 0] = 0

    m = bm.Basemap(projection='cyl', llcrnrlon=west,
                   llcrnrlat=south, urcrnrlon=east,
                   urcrnrlat=north, lon_0=0.5*(west+east),
                   lat_0=0.5*(south+north), resolution='f')

    x_grid, y_grid = m(lon_grid, lat_grid)
    m.pcolormesh(x_grid, y_grid, bathy_grid, cmap=plt.get_cmap('binary_r'))

    m.drawcoastlines()
    m.fillcontinents()

    meridians = np.round(np.linspace(west, east, 7))
    m.drawmeridians(meridians, labels=[0, 0, 0, 1])
    parallels = np.round(np.linspace(south, north, 3))
    m.drawparallels(parallels, labels=[1, 0, 0, 0])
    plt.title('Drake Passage and Scotia Sea Bathymetry')

    #Plotting bathymetry along a track example.
    lons = np.linspace(-60., -50., 1000)
    lats = np.linspace(-60., -55., 1000)
    bathy_track = interp_track(lons, lats)

    q, p = m(lons, lats)
    m.plot(p, q, 'r--', linewidth=2)

    plt.figure()
    plt.plot(lons, bathy_track)
    plt.xlabel('Longitude')
    plt.ylabel('Depth (m)')
