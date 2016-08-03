# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 15:18:09 2016

@author: jc3e13
"""

import os
import argparse
import numpy as np
import matplotlib
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


parser = argparse.ArgumentParser(description='Run w fitting on a float.')
parser.add_argument('--floatID', type=int, help='EM-APEX float ID number')
args = parser.parse_args()

Float = emapex.load(args.floatID, apply_w=False, apply_strain=False,
                    apply_iso=False, verbose=False)


# %% Test

# Ballast information provided to me by James Girton was contain in a matlab
# file called float_ballast.m. I have converted it into a dictionary and
# converted the unit of mass (last two columns) to kg from g.
ballast_info = {
    1632: [1, 1000., 5.102, 34.878, 31., 2.72e-6, 1.156, 27.885, 0.],
    1633: [1, 1000., 5.102, 34.878, 31., 2.72e-6, 1.156, 28.008, 0.],
    1634: [1, 1000., 5.102, 34.878, 31., 2.72e-6, 1.156, 27.948, 0.],
    1635: [1, 1000., 5.102, 34.878, 31., 2.72e-6, 1.156, 27.964, 0.],
    1636: [1, 1000., 5.102, 34.878, 31., 2.72e-6, 1.156, 27.948, 0.],
    3304: [1, 1500., 3.0, 34.6, 16., 3.67e-6, 1.156, 26.988, 0.],
    3305: [1, 1500., 3.0, 34.6, 16., 3.67e-6, 1.156, 27.016, 0.],
    3760: [1, 2000., 1.773, 34.766, 16., 3.6394E-6, 1.156, 27.146, 0.0455],
    3761: [1, 2000., 1.915, 34.770, 16., 3.7314E-6, 1.156, 27.201, 0.0455],
    3762: [1, 2000., 2.017, 34.766, 16., 3.7600E-6, 1.156, 27.193, 0.0455],
    3764: [1, 2000., 2.147, 34.756, 16., 3.8227E-6, 1.156, 27.201, 0.0455],
    3950: [1, 2000., 1.773, 34.766, 16., 3.7249E-6, 1.156, 27.226, 0.0455],
    3951: [1, 2000., 1.915, 34.770, 16., 3.6002E-6, 1.156, 27.239, 0.0455],
    3952: [1, 2000., 2.017, 34.766, 16., 3.6316E-6, 1.156, 27.210, 0.0455],
    4051: [1, 2000., 2.147, 34.756, 16., 3.6737E-6, 1.156, 27.171, 0.0455],
    3763: [1, 300., 13.5, 36.0, 16., 3.80e-6, 1.156, 27.900, 0.030],
    3765: [1, 300., 13.5, 36.0, 16., 3.69e-6, 1.156, 27.900, 0.030],
    3766: [1, 300., 13.5, 36.0, 16., 3.73e-6, 1.156, 27.900, 0.030],
    # DIMES FLOATS
    3767: [1, 2000., 0.400, 34.71, 16., 3.6736e-6, 1.156, 27.178, 0.049],
    4086: [1, 2000., 0.400, 34.71, 16., 3.7100e-6, 1.156, 27.177, 46.],
    4087: [1, 2000., 0.400, 34.71, 16., 3.6700e-6, 1.156, 27.175, 33.],
    4088: [1, 2000., 0.400, 34.71, 16., 3.6700e-6, 1.156, 27.154, 0.],
    4089: [1, 2000., 0.400, 34.71, 16., 3.6700e-6, 1.156, 27.136, 0.],
    4090: [1, 2000., 0.400, 34.71, 16., 3.6700e-6, 1.156, 27.197, 0.],
    4812: [1, 2000., 0.400, 34.71, 16., 3.6700e-6, 1.156, 27.197, 0.],
    4813: [1, 2000., 0.400, 34.71, 16., 3.6700e-6, 1.156, 27.197, 0.],
    4814: [1, 2000., 0.400, 34.71, 16., 3.6700e-6, 1.156, 27.197, 0.],
    4815: [1, 2000., 0.400, 34.71, 16., 3.6700e-6, 1.156, 27.197, 0.],
    4976: [1, 2000., 0.400, 34.71, 16., 3.6700e-6, 1.156, 27.197, 0.],
    4977: [1, 2000., 0.400, 34.71, 16., 3.6700e-6, 1.156, 27.197, 0.],
    4594: [1, 2000., 0.400, 34.71, 16., 3.6700e-6, 1.156, 27.197, 0.],
    4595: [1, 2000., 0.400, 34.71, 16., 3.6700e-6, 1.156, 27.197, -0.050],
    4596: [1, 2000., 0.400, 34.71, 16., 3.6700e-6, 1.156, 27.197, 0.],
    4597: [1, 2000., 0.400, 34.71, 16., 3.6700e-6, 1.156, 27.197, -0.070],
    6478: [1, 2000., 0.400, 34.71, 39., 3.8100e-6, 0.8, 27.2626, 0.],
    6480: [1, 2000., 0.400, 34.71, 39., 3.6700e-6, 0.8, 27.197, 0.],
    6481: [1, 2000., 0.400, 34.71, 39., 3.6700e-6, 0.8, 27.197, 0.],
    6625: [1, 2000., 0.400, 34.71, 39., 3.6700e-6, 0.8, 27.197, 0.],
    6626: [1, 2000., 0.400, 34.71, 39., 3.6700e-6, 0.8, 27.197, 0.]
    }

# Limits found by looking at profiling strategy.
hpid_lims = {
    3767: [(1, 46)],
    4086: [(1, 20)],
    4087: [(1, 15)],
    4089: [(1, 600)],  # Changed from (1, 1)
    4090: [(1, 15)],
    4594: [(1, 70), (98, 370)],
    4595: [(1, 60), (140, 170)],
    4596: [(1, 600)],  # Changed from (1, 1)
    4597: [(1, 204)],
    4812: [(1, 26)],
    4813: [(1, 76)],
    4814: [(1, 24)],
    4815: [(1, 38)],
    4976: [(1, 600)],
    4977: [(1, 600)],
    6478: [(1, 50), (76, 114)],
    6480: [(1, 50)],
    6481: [(1, 36)],
    6625: [(1, 40), (160, 240)],
    6626: [(1, 50)]
    }

hpids = np.hstack([np.arange(*hpid_lim) for hpid_lim in hpid_lims[Float.floatID]])

profiles = 'updown'

# Initial parameters
V_0 = 2.615e-2
CA = 2e-2
alpha_p = ballast_info[Float.floatID][5]
p_0 = ballast_info[Float.floatID][1]
alpha_ppos = ballast_info[Float.floatID][6]/1e6
ppos_0 = ballast_info[Float.floatID][4]
M = ballast_info[Float.floatID][7] + ballast_info[Float.floatID][8]

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

if profiles == 'all':
    print('Fitting drag to up and down profiles together.')
    p0 = np.array([V_0, CA, alpha_p, p_0, alpha_ppos, ppos_0, M])
    pfixed = [None, None, None, p_0, None, ppos_0, M]
if profiles == 'updown':
    print('Fitting drag to up and down profiles separately.')
    p0 = np.array([V_0, CA, CA, alpha_p, p_0, alpha_ppos, ppos_0, M])
    pfixed = [None, None, None, None, p_0, None, ppos_0, M]

wfi = vvf.fitter(Float, p0, pfixed, profiles=profiles, save_name=save_name,
                 N_bootstrap=200, method='Nelder-Mead')
print("Fitting completed, starting assessment.")

Float.apply_w_model(wfi)
vvf.assess_w_fit(Float, save_figures=True, save_id=str(Float.floatID),
                 save_dir=fsdir)
