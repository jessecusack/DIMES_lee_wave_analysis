# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 12:15:02 2016

@author: Jesse
"""

import numpy as np
import scipy.signal as sig
from scipy.io import loadmat
import matplotlib
import matplotlib.pyplot as plt
import os
import fnmatch
from my_savefig import my_savefig

import argparse

# Figure save path.
fsdir = '../figures/hf_noise'
if not os.path.exists(fsdir):
    os.makedirs(fsdir)

# Universal figure font size.
matplotlib.rc('font', **{'size': 8})

parser = argparse.ArgumentParser(
    description='Analyse high frequency noise content of pressure sensor.')

parser.add_argument('--floatID', type=int, help='EM-APEX float ID number')
parser.add_argument('--dirpath', type=str, help='float data directory')
args = parser.parse_args()

floatID = args.floatID
dirpath = args.dirpath

# Find data files.
ctd_files = []
searchstr = '*{}*ctd.mat'.format(floatID)
for root, dirnames, filenames in os.walk(dirpath):
    for filename in fnmatch.filter(filenames, searchstr):
        nameparts = filename.split('-')
        hpid = int(nameparts[2])
        fullname = os.path.join(root, filename)
        ctd_files.append(fullname)

np.sort(ctd_files)

Npfls = len(ctd_files)

# Spectra things.
nperseg = 2**10
noverlap = nperseg/2
f = np.abs(np.fft.fftfreq(nperseg)[:nperseg/2+1])
dtm = 20.  # Average time step.
use = (f < 1./dtm) & (f > 0)
f = f[use]
Nfreqs = f.size

# Initialise arrays.
PP = np.NaN*np.zeros((Nfreqs, Npfls))
Pwf = np.NaN*np.zeros_like(PP)
hpid = np.NaN*np.zeros((Npfls,))

for i, ctd_file in enumerate(ctd_files):
    data = loadmat(ctd_file, squeeze_me=True)

    nout_ctd = data['nout_ctd']

    if nout_ctd < 3:  # Not enough data.
        continue

    UXT = data['UXT']
    P = data['P']
#    S = data['S']
    hpid[i] = data['hpid']
#    T = data['T']
#    ctdtime = data['ctdtime']
#    pcb = data['pcb']
#    pca = data['pca']

    wf = -np.gradient(P)/np.gradient(UXT)

    dt = 1.
    UXTi = np.arange(UXT[0], UXT[-1], dt)
    Pi = np.interp(UXTi, UXT, P)
    wfi = np.interp(UXTi, UXT, wf)

    __, PPi = sig.welch(Pi, nperseg=nperseg, noverlap=noverlap)
    __, Pwfi = sig.welch(wfi, nperseg=nperseg, noverlap=noverlap)

    PP[:, i] = PPi[use]
    Pwf[:, i] = Pwfi[use]

# %%

PPm = np.ma.masked_invalid(PP)
Pwfm = np.ma.masked_invalid(Pwf)
hpidm = np.ma.masked_invalid(hpid)
logf = np.log10(f)

fig, axs = plt.subplots(2, 1, sharex='col')
CPP = axs[0].pcolormesh(hpidm, logf, np.log10(PPm), cmap=plt.get_cmap('viridis'))
Cwfm = axs[1].pcolormesh(hpidm, logf, np.log10(Pwfm), cmap=plt.get_cmap('viridis'))

axs[0].set_xlim(np.min(hpidm), np.max(hpidm))
axs[0].set_ylim(np.min(logf), np.max(logf))
axs[1].set_ylim(np.min(logf), np.max(logf))

cb0 = plt.colorbar(CPP, ax=axs[0])
cb1 = plt.colorbar(Cwfm, ax=axs[1])

fig.suptitle("{}".format(floatID))
axs[1].set_xlabel('Half profile number')
axs[0].set_ylabel('log$_{10}$ Frequency (Hz)')
axs[1].set_ylabel('log$_{10}$ Frequency (Hz)')
cb0.set_label('log$_{10}$ PSD($P$) (dbar$^2$ Hz$^{-1}$)')
cb1.set_label('log$_{10}$ PSD($w_f$) (m$^2$ s$^{-2}$ Hz$^{-1}$)')

my_savefig(fig, '{}_PSD'.format(floatID), fsdir)
plt.close(fig)
