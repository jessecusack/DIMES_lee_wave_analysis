# -*- coding: utf-8 -*-
"""
Created on Tue Feb  3 10:40:31 2015

@author: jc3e13
"""

import numpy as np
import scipy as sp
import sys
import os
import matplotlib.pyplot as plt

lib_path = os.path.abspath('../modules')
if lib_path not in sys.path:
    sys.path.append(lib_path)

import emapex
import utils


try:
    print("Floats {} and {} exist!.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.load(4976)
    E77 = emapex.load(4977)


# %%

# The portions of the profiles that contain the wave. Found by eye.
zlims = {4976: {29: (-1600, -200),
                30: (-1000, -200),
                31: (-1600, -600),
                32: (-1600, -400)},
         4977: {24: (-1600, -200),
                25: (-1400, -600),
                26: (-1600, -600),
                27: (-1200, -200)}}

hpids_76 = np.array([29, 30, 31, 32])
hpids_77 = np.array([24, 25, 26, 27])


# Detrend degre
deg = 2

plt.figure()

for Float, hpids in zip([E76, E77], [hpids_76, hpids_77]):
    for hpid in hpids:

        zmin, zmax = zlims[Float.floatID][hpid]
        pfl = Float.get_profiles(hpid)

        use = (pfl.z > zmin) & (pfl.z < zmax)
        useef = (pfl.zef > zmin) & (pfl.zef < zmax)

        t = pfl.UTC[use]
        tef = pfl.UTCef[useef]
        z = pfl.zef[useef]
        u = pfl.U[useef]
        v = pfl.V[useef]
        w = pfl.Ww[use]

        w = np.interp(tef, t, w)

        u = utils.nan_detrend(z, u, deg)
        v = utils.nan_detrend(z, v, deg)

        plt.plot(u*w, z, v*w, z)
