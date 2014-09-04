# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 14:16:02 2014

@author: jc3e13
"""

import numpy as np


def chunk(x, x_range, y):
    """x should be monotonically increasing.
    y should be of the same size as x. (it will still work if this is not the
    case but beware indexing problems)
    """

    if len(x_range) != 2:
        raise ValueError('x_range must be a sequence of two numbers only.')

    s = slice(*np.searchsorted(x, x_range))

    return x[s], y[s]


def window(x, y, width=25, overlap=None, expansion=0.):
    """Can't currently deal with expanding windows..."""

    if x.size != y.size:
        raise ValueError('x and y must be of equal size.')

    x_nans = np.isnan(x)
    x = x[~x_nans]
    y = y[~x_nans]

    if overlap is None:
        overlap = 0.

    step = width - overlap

    left = np.arange(x[0], x[-1], step)
    right = left + width

    bins = np.transpose(np.vstack((left, right)))

    vals = np.asarray([chunk(x, b, y) for b in bins])

    return vals


def window_func(x, y, func, width=25, overlap=None):

    vals = window(x, y, width=width, overlap=overlap)

    X = np.empty_like(vals[:, 0])
    Y = np.empty_like(vals[:, 0])

    for i, row in enumerate(vals):
        X[i] = func(row[0])
        Y[i] = func(row[1])

    return X, Y


def window_func_2(x, y, func, func_args=None, width=25, overlap=None):

    vals = window(x, y, width=width, overlap=overlap)

    X = np.empty_like(vals[:, 0])
    Y = np.empty_like(vals[:, 0])

    for i, row in enumerate(vals):
        X[i] = np.mean(row[0])
        Y[i] = func(row[0], row[1], *func_args)

    return X, Y