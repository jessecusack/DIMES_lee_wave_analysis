# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 12:58:56 2014

@author: jc3e13
"""

import numpy as np
import sys
import gsw
from scipy.integrate import odeint
import matplotlib.pyplot as plt

lib_path = '/noc/users/jc3e13/emapex/python'
pymc3_path = '/noc/users/jc3e13/envs/my_root/lib/python2.7/site-packages/pymc-3.0-py2.7.egg'
if lib_path not in sys.path:
    sys.path.append(lib_path)

# We want pymc 2.3
if pymc3_path in sys.path:
    sys.path.remove(pymc3_path)

import pymc
import emapex
import utils
import gravity_waves as gw

try:
    print("Floats {} and {} exist!".format(E76.floatID, E77.floatID))
except NameError:
    E76 = emapex.EMApexFloat('/noc/users/jc3e13/data/EM-APEX/allprofs11.mat', 4976)
    E76.apply_w_model('/noc/users/jc3e13/data/EM-APEX/4976_fix_p0k0M_fit_info.p')
    E76.apply_strain('/noc/users/jc3e13/data/EM-APEX/4976_N2_ref_300dbar.p')
    E76.apply_isopycnal_displacement('/noc/users/jc3e13/data/EM-APEX/srho_4976_100mbin.p')
    E77 = emapex.EMApexFloat('/noc/users/jc3e13/data/EM-APEX/allprofs11.mat', 4977)
    E77.apply_w_model('/noc/users/jc3e13/data/EM-APEX/4977_fix_p0k0M_fit_info.p')
    E77.apply_strain('/noc/users/jc3e13/data/EM-APEX/4977_N2_ref_300dbar.p')
    E77.apply_isopycnal_displacement('/noc/users/jc3e13/data/EM-APEX/srho_4977_100mbin.p')


def drdt(r, t, phi_0, U, Wf_pvals, k, l, m, om, N, f, phase):
    x = r[0]
    y = r[1]
    z = r[2]

    Wf_g = Wf_pvals[0]
    Wf_0 = Wf_pvals[1]

#    U = np.polyval(U_pvals, z)

    om2 = om**2
    f2 = f**2
    K2 = k**2 + l**2 + m**2

    dxdt = U + np.real(((k*om + 1j*l*f)/(om2 - f2))*phi_0*np.exp(1j*(k*x + l*y + m*z - om*t + phase)))
    dydt = np.real(((l*om - 1j*k*f)/(om2 - f2))*phi_0*np.exp(1j*(k*x + l*y + m*z - om*t + phase)))
    dzdt = (Wf_0 + np.real(((-om*K2)/((N**2 - f2)*m))*phi_0*np.exp(1j*(k*x + l*y + m*z - om*t + phase))))/(1 - Wf_g)

    return np.array([dxdt, dydt, dzdt])


def wave_vel(r, t, phi_0, k, l, m, om, N, f, phase):
    x = r[..., 0]
    y = r[..., 1]
    z = r[..., 2]

    om2 = om**2
    f2 = f**2
    K2 = k**2 + l**2 + m**2

    u_x = np.real(((k*om + 1j*l*f)/(om2 - f2))*phi_0*np.exp(1j*(k*x + l*y + m*z - om*t + phase)))
    u_y = np.real(((l*om - 1j*k*f)/(om2 - f2))*phi_0*np.exp(1j*(k*x + l*y + m*z - om*t + phase)))
    u_z = np.real(((-om*K2)/((N**2 - f2)*m))*phi_0*np.exp(1j*(k*x + l*y + m*z - om*t + phase)))

    return (np.vstack((u_x, u_y, u_z))).T


def buoy(r, t, phi_0, k, l, m, om, N, f, phase):
    x = r[:, 0]
    y = r[:, 1]
    z = r[:, 2]

    om2 = om**2
    N2 = N**2

    b = np.real((1j*m*N2/(N2 - om2))*phi_0*np.exp(1j*(k*x + l*y + m*z - om*t + phase)))

    return b


def wave_profile(zf, X, Y, Z, phase):
    U = 0.5
    V = -0.0
    f = gsw.f(-57.5)
    N = 1.8e-3

    # Float change in buoyancy with velocity.
    Wf_pvals = np.polyfit([0., 0.06], [0.14, 0.12], 1)

    # Wave parameters
    W_0 = 0.17
    k = 2*np.pi/X
    l = 2*np.pi/Y
    m = 2*np.pi/Z

    om = gw.omega(N, k, m, l, f) + k*U + l*V
    phi_0 = W_0*(N**2 - f**2)*m/(om*(k**2 + l**2 + m**2))

    args = (phi_0, U, Wf_pvals, k, l, m, om, N, f, phase)
    uargs = (phi_0, k, l, m, om, N, f, phase)

    # Integration parameters.
    dt = 20.
    t_0 = 0.
    t_1 = 13000.
    t = np.arange(t_0, t_1, dt)

    # Initial conditions.
    x_0 = 0.
    y_0 = 0.
    z_0 = zmin
    r_0 = np.array([x_0, y_0, z_0])

    # This integrator calls FORTRAN odepack to solve the problem.
    r = odeint(drdt, r_0, t, args)
    u = wave_vel(r, t, *uargs)
    u[:, 0] += U
    u[:, 1] += V
    b = buoy(r, t, *uargs)

    um = np.interp(zf, r[:, 2], u[:, 0])
    vm = np.interp(zf, r[:, 2], u[:, 1])
    wm = np.interp(zf, r[:, 2], u[:, 2])
    bm = 250.*np.interp(zf, r[:, 2], b)

#    # Variable to return.
#    var_dict = {'w': u[:, 2], 'u': u[:, 0], 'v': u[:, 1], 'b': b}
#    var = var_dict['w']
#    ivar = np.interp(zf, r[:, 2], var)

    return np.hstack((um, vm, wm, bm))

N = 10
N0_76 = 15
N0_77 = 10
E76_hpids = np.arange(N0_76, N0_76+N)
E77_hpids = np.arange(N0_77, N0_77+N)

dz = 1.
z = np.arange(-1500, 0, dz)
rho = []

pfls = np.hstack((E76.get_profiles(E76_hpids), E77.get_profiles(E77_hpids)))

for pfl in pfls:
    rho.append(pfl.interp(z, 'z', 'rho_1'))

rho = np.transpose(np.asarray(rho))
mrho = np.mean(rho, axis=-1)

pfl26 = E77.get_profiles(26)
srhop = utils.nan_interp(pfl26.z, z, mrho)

zmax = -650.
use = ~np.isnan(pfl26.z) & (pfl26.z < zmax)
zf = pfl26.z[use]
wf = pfl26.Ww[use]
uf = pfl26.U_abs[use]
vf = pfl26.V_abs[use]
bf = utils.nan_detrend(zf, (-gsw.grav(pfl26.lat_start, pfl26.P)*(pfl26.rho_1 - srhop)/1031.)[use])
zmin = np.min(zf)

data_stack = np.hstack((uf, vf, wf, 250*bf))


def model():

    # Priors.
    sig = pymc.Uniform('sig', 0.0, 5., value=0.01)
    X = pymc.Uniform('X', -2e4, 2e4, value=5e3)
    Y = pymc.Uniform('Y', -1e5, 1e5, value=1e4)
    Z = pymc.Uniform('Z', -2e4, 2e4, value=6e3)
    phase = pymc.Uniform('phase', 0., 2.*np.pi, value=0.)

    @pymc.deterministic()
    def wave_model(zf=zf, X=X, Y=Y, Z=Z, phase=phase):
        return wave_profile(zf, X, Y, Z, phase)

    # Likelihood
    y = pymc.Normal('y', mu=wave_model, tau=1./sig**2, value=data_stack, observed=True)

    return locals()

M = pymc.MCMC(model(), db='pickle', dbname='trace.p')
samples = 50000
burn = 20000
thin = 10
M.sample(samples, burn, thin)
pymc.Matplot.plot(M, common_scale=False)

plt.figure()
plt.plot(data_stack)
plt.plot(wave_profile(zf, np.median(M.trace('X')[:]), np.median(M.trace('Y')[:]),
                      np.median(M.trace('Z')[:]), np.median(M.trace('phase')[:])))
plt.savefig('output.png')

plt.figure()
plt.plot(data_stack)
plt.plot(wave_profile(zf, 5000., 15000, 8000, 0.))
