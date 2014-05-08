# -*- coding: utf-8 -*-
"""
Created on Sun Feb 02 17:45:16 2014

Functions for reading Smith and Sandwell bathymetry.

@author: jc3e13
"""

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm


def read_grid(lon_lat, file_path=None):
    """Input format [lonmin lonmax lat_min lat_max].
                     west   east   south   north

    Output:
      bathy_grid, lon_grid, lat_grid

    TODO: More input checks.
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


if __name__ == '__main__':

    lon_lat = np.array([-72, -27, -68, -47])
    west, east, south, north = lon_lat

    lon_grid, lat_grid, bathy_grid = read_grid(lon_lat)
    bathy_grid[bathy_grid > 0] = 0

    m = bm.Basemap(projection='cyl', llcrnrlon=west,
                   llcrnrlat=south, urcrnrlon=east,
                   urcrnrlat=north, lon_0=0.5*(west+east),
                   lat_0=0.5*(south+north), resolution='f')

    x, y = m(lon_grid, lat_grid)
    m.pcolormesh(x, y, bathy_grid, cmap=plt.get_cmap('binary_r'))

    m.drawcoastlines()
    m.fillcontinents()

    meridians = np.round(np.linspace(west, east, 7))
    m.drawmeridians(meridians, labels=[0, 0, 0, 1])
    parallels = np.round(np.linspace(south, north, 3))
    m.drawparallels(parallels, labels=[1, 0, 0, 0])
    plt.title('Drake Passage and Scotia Sea Bathymetry')
