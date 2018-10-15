# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 14:36:44 2016

@author: jc3e13
"""

import os
import emapex


for floatID in emapex.FIDS_DIMES:
    command = 'source activate dimes && python fit_w_model_alt.py --floatID ' + str(floatID) + ' &'
    print(command)
    os.system(command)
