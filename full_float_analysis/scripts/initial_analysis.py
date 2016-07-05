# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 12:21:44 2016

@author: jc3e13
"""

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from glob import glob
from scipy import io

import emapex
import utils
from my_savefig import my_savefig

# Figure save path.
sdir = '../figures/initial_analysis'
if not os.path.exists(sdir):
    os.makedirs(sdir)
# Universal figure font size.
matplotlib.rc('font', **{'size': 8})

files = glob('/noc/users/jc3e13/storage/DIMES/EM-APEX/allprof*.mat')
fids = np.array([])
for f in files:
    fid = io.loadmat(f, squeeze_me=True, variable_names='flid')['flid']
    fids = np.hstack((fids, fid[~np.isnan(fid)]))

fids = fids.astype(np.uint16)
ufids = np.unique(fids)

print("Floats found:")
for ufid in ufids:
    print("{}".format(ufid))

for ufid in ufids:
    Float = emapex.load(ufid, apply_w=False, apply_strain=False,
                        apply_iso=False)

    UTC = Float.UTC.flatten(order='F')
    UTC[UTC < 733774.] = np.NaN
    t = utils.datenum_to_datetime(UTC)
    z = Float.z.flatten(order='F')
    z[z > -1.] = np.NaN

    fig, ax = plt.subplots(1, 1, figsize=(6.5, 3))
    ax.plot(t, z)
    ax.set_xlabel('Time')
    ax.set_ylabel('Height')
    ax.set_title("{}".format(ufid))
    my_savefig(fig, "z_{}".format(ufid), sdir)
    plt.close(fig)
