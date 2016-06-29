# -*- coding: utf-8 -*-
"""
Created on Fri May 15 13:39:41 2015

@author: jc3e13
"""

import numpy as np


# Typically measured values.

w0 = 0.15
N = 2e-3
b0 = 5e-4
om0 = 1.5e-3
V = 0.
U = 0.4
u0 = 0.1
v0 = 0.1

# Squared values.

w2 = w0**2
N2 = N**2
b2 = b0**2
om2 = om0**2
V2 = V**2
U2 = U**2

# Simple case first.

kh1 = np.sqrt((w2*N2**2 + om2*b2)/(U2*b2))

# Complicated case.

#theta = np.arctan2(v0, u0)
#
#den = b2*np.cos(theta)**2*(U + u0)**2 + b2*np.sin(theta)**2*(V + v0)**2 \
#    + 2*b2*np.cos(theta)*np.sin(theta)*(U + u0)*(V + v0)
#num = w2*N2**2 + om2*b2
#
#kh2 = np.sqrt(num/den)
