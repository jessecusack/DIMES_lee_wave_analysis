# -*- coding: utf-8 -*-
"""
Created on Sun Feb 02 15:48:02 2014

@author: jc3e13
"""

import numpy as np
import datetime as dt


def datenum_to_datetime(datenum):
    """
    Convert a MATLAB datenum, or an array of datenums into a python datetime
    object or array of objects.

    If conversion fails the output is nan.

    TODO: rewrite function more explicitly, not clear what second exceptions
    are dealing with.
    """

    def convert(datenum):
        try:
            return dt.datetime.fromordinal(int(datenum)) + \
                dt.timedelta(days=datenum % 1) - dt.timedelta(days=366)
        except ValueError:
            return np.nan

    try:

        pydt = np.array([convert(el) for el in datenum])

    except TypeError:

        pydt = convert(datenum)

    return pydt
