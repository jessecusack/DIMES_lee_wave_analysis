# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 11:00:28 2015

@author: jc3e13
"""

import numpy as np
#from scipy import optimize as op
import matplotlib.pyplot as plt
import pymc
import triangle

np.random.seed(12345)

# True parameters.
phi_0 = 0.5
lx = -3000.
ly = -3000.
lz = -3000.
phase_0 = 1.
N = 2e-3
U = 0.5

k = np.pi*2./lx
l = np.pi*2./ly
m = np.pi*2./lz
om = np.sqrt(N**2*(k**2 + l**2)/(k**2 + l**2 + m**2))


#def model(data, amp, k, m, phase):
#    x, z = data
#    return amp*np.sin(k*x + m*z + phase)


def u_model(data, phi_0, k, l, m, phase_0, N, U):
    x, y, z, t = data
    om = np.sqrt(N**2*(k**2 + l**2)/(k**2 + l**2 + m**2))
    u_0 = phi_0*k/om
    phase = 1j*(k*x + l*y + m*z - (om + k*U)*t + phase_0)
    return np.real(u_0*np.exp(phase))


def v_model(data, phi_0, k, l, m, phase_0, N, U):
    x, y, z, t = data
    om = np.sqrt(N**2*(k**2 + l**2)/(k**2 + l**2 + m**2))
    v_0 = phi_0*l/om
    phase = 1j*(k*x + l*y + m*z - (om + k*U)*t + phase_0)
    return np.real(v_0*np.exp(phase))


def w_model(data, phi_0, k, l, m, phase_0, N, U):
    x, y, z, t = data
    om = np.sqrt(N**2*(k**2 + l**2)/(k**2 + l**2 + m**2))
    w_0 = -phi_0*m*om/(N**2 - om**2)
    phase = 1j*(k*x + l*y + m*z - (om + k*U)*t + phase_0)
    return np.real(w_0*np.exp(phase))


def b_model(data, phi_0, k, l, m, phase_0, N, U):
    x, y, z, t = data
    om = np.sqrt(N**2*(k**2 + l**2)/(k**2 + l**2 + m**2))
    b_0 = phi_0*1j*m*N**2/(N**2 - om**2)
    phase = 1j*(k*x + l*y + m*z - (om + k*U)*t + phase_0)
    return 250.*np.real(b_0*np.exp(phase))


def model(data, phi_0, k, l, m, phase_0, N, U):

    u = u_model(data, phi_0, k, l, m, phase_0, N, U)
    v = v_model(data, phi_0, k, l, m, phase_0, N, U)
    w = w_model(data, phi_0, k, l, m, phase_0, N, U)
    b = b_model(data, phi_0, k, l, m, phase_0, N, U)

    return np.hstack((u, v, w, b))


# Fake data.
Np = 100
noise_amp = 0.5
x = np.linspace(0., np.abs(lx), Np)
y = np.linspace(0., np.abs(ly), Np)
z = np.linspace(0., np.abs(lz), Np)
t = np.linspace(0., 5000., Np)

u = u_model([x, y, z, t], phi_0, m, l, k, phase_0, N, U)
v = v_model([x, y, z, t], phi_0, m, l, k, phase_0, N, U)
w = w_model([x, y, z, t], phi_0, m, l, k, phase_0, N, U)
b = b_model([x, y, z, t], phi_0, m, l, k, phase_0, N, U)
uvwb = model([x, y, z, t], phi_0, m, l, k, phase_0, N, U)

u_noise = noise_amp*np.random.randn(Np)
v_noise = noise_amp*np.random.randn(Np)
w_noise = noise_amp*np.random.randn(Np)
b_noise = noise_amp*np.random.randn(Np)

u_noisy = u + u_noise
v_noisy = v + v_noise
w_noisy = w + w_noise
b_noisy = b + b_noise

uvwb_noise = np.hstack((u_noise, v_noise, w_noise, b_noise))
uvwb_noisy = uvwb + uvwb_noise

# Solve using MCMC sampling.
def pymc_model(data, uvwb_noisy, noise_amp, N, U):

    sig = noise_amp
    phi_0 = pymc.Uniform('phi_0', 0, 50., value=4.)
    X = pymc.Uniform('X', -100000., 100000., value=-5000.)
    Y = pymc.Uniform('Y', -100000., 100000., value=-1000.)
    Z = pymc.Uniform('Z', -100000., 100000., value=-10000.)
    phase_0 = pymc.Uniform('phase_0', -100., 100., value=0.)

    @pymc.deterministic()
    def mmodel(phi_0=phi_0, X=X, Y=Y, Z=Z, phase_0=phase_0):
        k = np.pi*2./X
        l = np.pi*2./Y
        m = np.pi*2./Z
        return model(data, phi_0, k, l, m, phase_0, N, U)

    # Likelihood
    uvwb_fit = pymc.Normal('uvwb_fit', mu=mmodel, tau=1./sig**2,
                           value=uvwb_noisy, observed=True)
    return locals()

# Run the sampler.
M = pymc.MCMC(pymc_model([x, y, z, t], uvwb_noisy, noise_amp, N, U))
samples = 300000
burn = 250000
thin = 10
M.sample(samples, burn, thin)

# Plot fit comparison.
fig, axs = plt.subplots(1, 4, sharey='col')
axs[0].plot(u, z, color='black')
axs[0].errorbar(u_noisy, z, xerr=noise_amp, marker='o', color='black',
                linestyle='none')
axs[0].set_xlabel('$u$')
axs[1].plot(v, z, color='black')
axs[1].errorbar(v_noisy, z, xerr=noise_amp, marker='o', color='black',
                linestyle='none')
axs[1].set_xlabel('$u$')
axs[2].plot(w, z, color='black')
axs[2].errorbar(w_noisy, z, xerr=noise_amp, marker='o', color='black',
                linestyle='none')
axs[2].set_xlabel('$w$')
axs[3].plot(b, z, color='black')
axs[3].errorbar(b_noisy, z, xerr=noise_amp, marker='o', color='black',
                linestyle='none')
axs[3].set_xlabel('$b$')

Ns = (samples - burn)/thin

for i in xrange(0, Ns, 20):
    params = [M.trace('phi_0')[i], np.pi*2./M.trace('X')[i],
              np.pi*2./M.trace('Y')[i], np.pi*2./M.trace('Z')[i],
              M.trace('phase_0')[i], N, U]
    axs[0].plot(u_model([x, y, z, t], *params), z, color='red', alpha=0.03)
    axs[1].plot(v_model([x, y, z, t], *params), z, color='red', alpha=0.03)
    axs[2].plot(w_model([x, y, z, t], *params), z, color='red', alpha=0.03)
    axs[3].plot(b_model([x, y, z, t], *params), z, color='red', alpha=0.03)

triangle.corner(np.transpose(np.asarray([M.trace('phi_0')[:], M.trace('X')[:],
                                         M.trace('Y')[:], M.trace('Z')[:],
                                         M.trace('phase_0')[:]])),
                labels=['$\phi_0$', '$\lambda_x$', '$\lambda_y$',
                        '$\lambda_z$', 'phase'],
                truths=[phi_0, lx, ly, lz, phase_0],
                quantiles=[.05, .95])
