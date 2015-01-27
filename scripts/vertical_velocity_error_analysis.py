# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 16:27:36 2014

@author: jc3e13
"""

import numpy as np
import matplotlib.pyplot as plt
import emapex

try:
    print("Floats {} and {}.".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.load(4976)
    E77 = emapex.load(4977)