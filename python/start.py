# -*- coding: utf-8 -*-
"""
This script loads floats and imports useful modules into the namespace.

Created on Thu May 29 14:52:32 2014

@author: jc3e13
"""

import scipy as sp
import gsw
import emapex


E76 = emapex.EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4976)
E77 = emapex.EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4977)

E76.apply_w_model('../../data/EM-APEX/4976_fix_alphakM_fit_info.p')
E77.apply_w_model('../../data/EM-APEX/4977_fix_alphakM_fit_info.p')
