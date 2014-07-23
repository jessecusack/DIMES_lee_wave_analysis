# -*- coding: utf-8 -*-
"""
Created on Sun Feb 02 15:48:02 2014

@author: jc3e13
"""

import numpy as np
import datetime as dt


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
