# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 14:36:44 2016

@author: jc3e13
"""

import os
import numpy as np
import matplotlib
import fnmatch
from scipy.io import loadmat

import emapex
import vertical_velocity_fitter as vvf

# Figure save path.
fsdir = '../figures/all_fit_specs'
if not os.path.exists(fsdir):
    os.makedirs(fsdir)
psdir = '../processed_data/'
if not os.path.exists(psdir):
    os.makedirs(psdir)
# Universal figure font size.
matplotlib.rc('font', **{'size': 8})


Floats = emapex.load_DIMES()

# %% Test
basepath = '/noc/users/jc3e13/storage/DIMES/EM-APEX/'

for Float in Floats:
    floatfolder = str(Float.floatID) + 'a'
    dirpath = os.path.join(basepath, floatfolder)

    print('Searching for mission file.')
    filesdict = {}
    mis_file = None
    gps_file = None
    single_mis_file = False
    single_gps_file = False
    searchstr = '*{}*mis.mat'.format(Float.floatID)
    for root, dirnames, filenames in os.walk(dirpath):
        for filename in fnmatch.filter(filenames, searchstr):
            nameparts = filename.split('-')

            try:
                hpid = int(nameparts[2])
            except ValueError:
                if nameparts[2] == 'mis.mat':
                    single_mis_file = True
                    mis_file = os.path.join(root, filename)
                    continue

            filetype = nameparts[3].split('.')[0]
            fullname = os.path.join(root, filename)
            if hpid in filesdict.keys():
                filesdict[hpid][filetype] = fullname
            else:
                filesdict[hpid] = {filetype: fullname}

    # Initial parameters
    V_0 = 2.62e-2
    CA = 5e-2

    if single_mis_file:
        mis = loadmat(mis_file, squeeze_me=True)

        alpha_p = -mis['FloatAlpha'][0]
        p_0 = mis['PressureBallastPoint'][0]
        alpha_ppos = mis['PistonCountsPerCC'][0]/1e6
        ppos_0 = mis['PistonParkPosition'][0]
        M = mis['FloatMass'][0]/1000.
    else:
        key0 = filesdict.keys()[0]
        mis = loadmat(filesdict[key0]['mis'], squeeze_me=True)

        alpha_p = -mis['FloatAlpha']
        p_0 = mis['PressureBallastPoint']
        alpha_ppos = mis['PistonCountsPerCC']/1e6
        ppos_0 = mis['PistonParkPosition']
        M = mis['FloatMass']/1000.

    print("Fitting float {}".format(Float.floatID))
    print("Initial parameters:\n"
          "  V_0 = {:1.2e}\n"
          "  CA = {:1.2e}\n"
          "  alpha_p = {:1.2e}\n"
          "  p_0 = {:1.2e}\n"
          "  alpha_ppos = {:1.2e}\n"
          "  ppos_0 = {:1.0f}\n"
          "  M = {:1.2e}".format(V_0, CA, alpha_p, p_0, alpha_ppos, ppos_0, M))

    name = "{}_wfi.p".format(Float.floatID)
    save_name = os.path.join(psdir, name)

    p0 = np.array([V_0, CA, alpha_p, p_0, alpha_ppos, ppos_0, M])
    pfixed = [V_0, None, None, p_0, None, None, M]
    wfi = vvf.fitter(Float, p0, pfixed, save_name=save_name)
    print("Fitting completed, starting assessment.")

    Float.apply_w_model(wfi)
    vvf.assess_w_fit(Float, save_figures=True, save_id=str(Float.floatID),
                     save_dir=fsdir)
