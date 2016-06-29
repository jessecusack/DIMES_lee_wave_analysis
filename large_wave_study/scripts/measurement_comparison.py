# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 15:20:19 2016

@author: jc3e13
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import gravity_waves as gw
import detect_peaks as dp


def drdt(r, t, phi_0, k, l, m, N=2e-3, U=0.3, V=0., Wf=0.1, f=0., phase_0=0.):
    x = r[0]
    y = r[1]
    z = r[2]

    om = gw.omega(N, k, m, l, f)
    dxdt = U + gw.u(x, y, z, t, phi_0, k, l, m, om, f=f, U=U, phase_0=phase_0)
    dydt = V + gw.v(x, y, z, t, phi_0, k, l, m, om, f=f, U=U, phase_0=phase_0)
    dzdt = Wf + gw.w(x, y, z, t, phi_0, k, l, m, om, N, U=U, phase_0=phase_0)

    return np.array([dxdt, dydt, dzdt])


# Define the wave.
X = -1500.
Y = -1500.
Z = -600.
phi_0 = 0.004
N = 2e-3
f = 1.2e-4
U = 0.
V = 0.
Wf = -0.12
phase_0 = 0.
rho0 = 1025.

k = 2.*np.pi/X
l = 2.*np.pi/Y
m = 2.*np.pi/Z
om = gw.omega(N, k, m, l, f)
alpha = gw.alpha(k, m, l)
w_0 = gw.W_0(phi_0, m, om, N)
Edens = gw.Edens(w_0, k, m, l, rho0)
Efluxz = gw.Efluxz(w_0, k, m, N, l, f, rho0)
Mfluxz = gw.Mfluxz(phi_0, k, l, m, om, N, f, rho0)

t = np.arange(0., 50000.)
r0 = np.array([0., 0., 0.])
args = (phi_0, k, l, m, N, U, V, Wf, f, phase_0)
r = sp.integrate.odeint(drdt, r0, t, args=args)
x, y, z = r[:, 0], r[:, 1], r[:, 2]

vel = gw.wave_vel(r, t, phi_0, N, f, k, l, m, om, U=U, V=V, phase_0=phase_0)
u, v, w = vel[:, 0], vel[:, 1], vel[:, 2]
b = gw.buoy(r, t, phi_0, N, k, l, m, om, U=U, V=V, phase_0=phase_0)
p = gw.phi(r[:, 0], r[:, 1], r[:, 2], t, phi_0, k, l, m, om, U=U, V=V)

# Nash method of estimating pressure perturbation
# z should be increasing.
if z[0] > z[-1]:
    print("Downward Profile!\n")
    zud = np.flipud(z)
    bud = np.flipud(b)
    bi = sp.integrate.cumtrapz(bud, zud, initial=0.)
    bii = sp.integrate.cumtrapz(bi, zud, initial=0.)

    H = zud.max() - zud.min()

    pi = bi + (bii[0] - bii[-1])/H

    pi = np.flipud(pi)

else:
    bi = sp.integrate.cumtrapz(b, z, initial=0.)
    bii = sp.integrate.cumtrapz(bi, z, initial=0.)

    H = zud.max() - zud.min()

    pi = bi + (bii[0] - bii[-1])/H

nhs = 1.  # Non hydrostatic factor...
pi *= nhs

# My method of estimating w
wi = np.gradient(r[:, 2]) - Wf

# Measurements of momentum flux
idxs = dp.detect_peaks(w)
idx1, idx2 = idxs[0], idxs[1]
use = slice(idx1, idx2)
dT = t[idx2] - t[idx1]

uwbar = rho0*sp.integrate.trapz(wi[use]*u[use], t[use])/dT
vwbar = rho0*sp.integrate.trapz(wi[use]*v[use], t[use])/dT
pwbar = rho0*sp.integrate.trapz(wi[use]*pi[use], t[use])/dT

tau = np.sqrt(uwbar**2 + vwbar**2)

kinetic = 0.5*rho0*sp.integrate.trapz(u[use]**2 + v[use]**2 + w[use]**2, t[use])/dT
potential = 0.5*rho0*sp.integrate.trapz(b[use]**2/N**2, t[use])/dT
E = kinetic + potential

ER = pwbar/Efluxz

print("Reynolds stress:\n"
      " Measured components: ({:1.2f}, {:1.2f}) N m-2\n"
      " Measured magnitude: {:1.2f} N m-2\n"
      "".format(uwbar, vwbar, tau, Mfluxz))
print("Energy density:\n"
      " Measured {:1.2f} J m-3\n"
      " Actual {:1.2f} J m-3\n".format(E, Edens))
print("Vertical energy flux:\n"
      " Measured {:1.2f} W m-2\n"
      " Actual {:1.2f} W m-2\n"
      " Ratio {:1.2f}".format(pwbar, Efluxz, ER))


fig, axs = plt.subplots(1, 6, sharey=True, figsize=(16, 6))
axs[0].set_ylabel('$z$')
axs[0].plot(1e2*u, z, 'b')
axs[0].plot(1e2*u[use], z[use], 'g')
axs[0].set_xlabel('$u$ (cm s$^{-1}$)')
plt.setp(axs[0].xaxis.get_majorticklabels(), rotation=60)
axs[1].plot(1e2*v, z, 'b')
axs[1].plot(1e2*v[use], z[use], 'g')
axs[1].set_xlabel('$v$ (cm s$^{-1}$)')
plt.setp(axs[1].xaxis.get_majorticklabels(), rotation=60)
axs[2].plot(1e2*w, z, 'b')
axs[2].plot(1e2*wi, z, 'r')
axs[2].plot(1e2*wi[use], z[use], 'g')
axs[2].set_xlabel('$w$ (cm s$^{-1}$)')
plt.setp(axs[2].xaxis.get_majorticklabels(), rotation=60)
axs[3].plot(1e4*b, z, 'b')
axs[3].plot(1e4*b[use], z[use], 'g')
axs[3].set_xlabel('$b$ ($10^{-4}$ m s$^{-2}$)')
plt.setp(axs[3].xaxis.get_majorticklabels(), rotation=60)
axs[4].plot(1e2*p, z, 'b')
axs[4].plot(1e2*pi, z, 'r')
axs[4].plot(1e2*pi[use], z[use], 'g')
axs[4].set_xlabel('$\phi$ ($10^{-2}$ m$^2$ s$^{-2}$)')
plt.setp(axs[4].xaxis.get_majorticklabels(), rotation=60)
axs[5].plot(x, z, 'b')
axs[5].set_xlabel('$x$ (m)')
plt.setp(axs[5].xaxis.get_majorticklabels(), rotation=60)
plt.ylim(z.min(), z.max())

for ax in axs:
    ax.hlines([z[idx1], z[idx2]], *ax.get_xlim(), color='k')


# %% Error in hydrostatic approximation.

# Define the wave.
Xl = np.arange(-4000., 0., 200.)
Y = 1e7
Z = -2000.
phi_0l = np.arange(0.001, 0.101, 0.001)
N = 2e-3
f = 1.2e-4
U = 0.3
V = 0.
Wf = 0.12
phase_0 = 0.
rho0 = 1025.

Xg, phi_0g = np.meshgrid(Xl, phi_0l)
shp = Xg.shape
Ni = Xg.size

alphag = np.zeros_like(Xg)
Erat = np.zeros_like(Xg)

print("Progress:")
for il, (X, phi_0) in enumerate(zip(Xg.flatten(), phi_0g.flatten())):

    print("{:1.1f}%".format(100.*float(il)/float(Ni)))

    i, j = np.unravel_index(il, shp)

    k = 2.*np.pi/X
    l = 2.*np.pi/Y
    m = 2.*np.pi/Z
    om = gw.omega(N, k, m, l, f)
    alpha = gw.alpha(k, m, l)
    alphag[i, j] = alpha
    w_0 = gw.W_0(phi_0, m, om, N)
    Edens = gw.Edens(w_0, k, m, l, rho0)
    Efluxz = gw.Efluxz(w_0, k, m, N, l, f, rho0)

    t = np.arange(0., 50000.)
    r0 = np.array([0., 0., 0.])
    args = (phi_0, k, l, m, N, U, V, Wf, f, phase_0)
    r = sp.integrate.odeint(drdt, r0, t, args=args)
    x, y, z = r[:, 0], r[:, 1], r[:, 2]

    vel = gw.wave_vel(r, t, phi_0, N, f, k, l, m, om, U=U, V=V, phase_0=phase_0)
    u, v, w = vel[:, 0], vel[:, 1], vel[:, 2]
    b = gw.buoy(r, t, phi_0, N, k, l, m, om, U=U, V=V, phase_0=phase_0)
    p = gw.phi(r[:, 0], r[:, 1], r[:, 2], t, phi_0, k, l, m, om, U=U, V=V)

    # Nash method of estimating pressure perturbation
    if z[0] > z[-1]:
        print("Downward Profile!\n")
        zud = np.flipud(z)
        bud = np.flipud(b)
        bi = sp.integrate.cumtrapz(bud, zud, initial=0.)
        bii = sp.integrate.cumtrapz(bi, zud, initial=0.)

        H = zud.max() - zud.min()

        pi = bi + (bii[0] - bii[-1])/H

        pi = np.flipud(pi)

    else:
        bi = sp.integrate.cumtrapz(b, z, initial=0.)
        bii = sp.integrate.cumtrapz(bi, z, initial=0.)

        H = zud.max() - zud.min()

        pi = bi + (bii[0] - bii[-1])/H

    # Measurements of momentum flux
    idxs = dp.detect_peaks(w)
    idx1, idx2 = idxs[0], idxs[1]
    use = slice(idx1, idx2)
    dT = t[idx2] - t[idx1]

    pwbar = rho0*sp.integrate.trapz(w[use]*pi[use], t[use])/dT
    Erat[i, j] = pwbar/Efluxz

plt.figure()
plt.plot(alphag[0, :-4], np.mean(Erat[:, :-4], axis=0))
plt.xlabel('Aspect Ratio')
plt.ylabel('Energy flux error factor')