# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 12:38:27 2015

@author: jc3e13
"""

import itertools
import os

info = [('4976', ['31', '32'], [['-1600', '-600'], ['-1600', '-400']]),
        ('4977', ['26', '27'], [['-1600', '-600'], ['-1100', '-150']])]

opts = ['-5000', '5000']

for floatID, hpids, zranges in info:
    for hpid, zrange in zip(hpids, zranges):
        for XYZ in itertools.product(opts, opts, opts):
            zmin, zmax = zrange
            X0, Y0, Z0 = XYZ
            command = 'python run_pymc.py --floatID ' + floatID + ' --hpid ' \
                + hpid + ' --zrange ' + zmin + ' ' + zmax + ' --xyz ' + X0 + \
                ' ' + Y0 + ' ' + Z0 + ' &'
            print(command)
            os.system(command)
