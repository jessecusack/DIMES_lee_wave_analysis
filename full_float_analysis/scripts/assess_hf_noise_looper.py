# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 15:48:44 2016

@author: jc3e13
"""

import os


d = '/noc/soes/physics/jc3e13/DIMES/EM-APEX'
dirpaths = [os.path.join(d,o) for o in os.listdir(d) if os.path.isdir(os.path.join(d,o))]

for dirpath in dirpaths:
    floatID = os.path.basename(dirpath)[:4]
    command = "python assess_hf_noise.py --floatID " + floatID + " --dirpath " + dirpath + " &"
    print(command)
    os.system(command)
