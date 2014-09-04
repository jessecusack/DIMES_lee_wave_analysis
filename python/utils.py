# -*- coding: utf-8 -*-
"""
Created on Wed Aug  6 12:09:43 2014

@author: jc3e13

Too many little modules were cluttering the directory so I have shoved them
all into one miscellaneous 'utilities' module.

"""

import numpy as np
import datetime as dt


class Bunch(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def datenum_to_datetime(datenum):
    """
    Convert a MATLAB datenums into python datetimes.

    Generate smooth buoyancy frequency profile by applying the adiabatic
    levelling method of Bray and Fofonoff (1981).

    Parameters
    ----------
    datenum : array_like
        MATLAB datenumber which is the number of days since 0001-01-01.

    Returns
    -------
    pydt : ndarray
        Python datetime. See datetime module.


    """

    def convert(datenum):
        try:
            return dt.datetime.fromordinal(int(datenum)) + \
                dt.timedelta(days=datenum % 1) - dt.timedelta(days=366)
        except ValueError:
            return np.nan

    try:
        # If datenum is not iterable it will raise a TypeError. I could just
        # check whether it is iterable first... !
        pydt = np.array([convert(el) for el in datenum])

    except TypeError:

        pydt = convert(datenum)

    return pydt


def lldist(lon, lat):
    """Calculates the distance between longitude and latitude coordinates on a
    spherical earth with radius using the Haversine formula. Code modified from
    the MATLAB m_map toolbox function m_lldist.m.

    Parameters
    ----------
    lon : 1-D numpy.ndarray of floats.
        Longitude values. [degrees]
    lat : 1-D numpy.ndarray of floats.
        Latitude values. [degrees]

    Returns
    -------
    dist : 1-D numpy.ndarray of floats.
        Distance between lon and lat positions. [km]
    """

    pi180 = np.pi/180.
    earth_radius = 6378.137  # [km]

    lat1 = lat[:-1]*pi180
    lat2 = lat[1:]*pi180

    dlon = np.diff(lon)*pi180
    dlat = lat2 - lat1

    a = (np.sin(dlat/2.))**2 + np.cos(lat1)*np.cos(lat2)*(np.sin(dlon/2.))**2
    angles = 2.*np.arctan2(np.sqrt(a), np.sqrt(1-a))
    dist = earth_radius*angles
    return dist
