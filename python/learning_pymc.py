# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 10:56:19 2014

@author: jc3e13

Can I fit this model using MCMC...
"""

import numpy
import pymc


def Wrel_model(rho, k, P, drag, alpha, k0, kappa, rho0=1030):

    g = 9.81
    drho = rho/rho0
    aWrel2 = numpy.abs(g/drag*(1 - drho*(1 + alpha*(k - k0) - kappa*P)))
    return numpy.sqrt(aWrel2)

# Import real data...
k = numpy.genfromtxt('pc.txt', dtype=None, delimiter=',')
P = numpy.genfromtxt('P.txt', dtype=None, delimiter=',')
rho = numpy.genfromtxt('dens.txt', dtype=None, delimiter=',')
Wf = numpy.genfromtxt('Wf.txt', dtype=None, delimiter=',')
aWf = numpy.abs(Wf)

# Estimate priors...
drag = pymc.Uniform('drag', 0., 10000., value=0.2)
alpha = pymc.Uniform('alpha', -100., 100., value=0.1)
k0 = pymc.Uniform('k0', -1000., 1000., value=50.)
kappa = pymc.Uniform('kappa', 0., 100., value=0.1)
rho0 = pymc.Uniform('rho0', 0., 4000., value=1030.)
sig = numpy.std(aWf)

# The model is a cost function, I think.


@pymc.deterministic
def costFunc(k=k, P=P, rho=rho, aWf=aWf, drag=drag, alpha=alpha, k0=k0,
             kappa=kappa, rho0=rho0):

    Wr = Wrel_model(rho, k, P, drag, alpha, k0, kappa, rho0)
    resid2 = (aWf - Wr)**2
    return resid2.sum()

# The observation is that the cost function is 0....
obs = pymc.Normal('obs', mu=costFunc, tau=1/sig**2, value=0., observed=True)
