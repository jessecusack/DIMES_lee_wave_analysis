# -*- coding: utf-8 -*-
"""
Created on Mon Jul 14 14:42:44 2014

@author: jc3e13


"""

import numpy as np


N_0 = 5.2  # Buoyancy frequency (angular).
b = 1300.  # e-folding scale of N with depth [m].
E_0 = 6.3e-5  # Internal wave energy parameter.
f_30 = 7.3e-5  # Coriolis frequency (angular) at 30N.


def H(j, j_star=3., N_sum=100000):

    # The number over which to sum if j_star is not 3.
    if j_star == 3.:

        # The factor 0.468043 comes from summing denominator from j = 1 to
        # j = 1e+8 using j_star = 3.
        return (j**2 + j_star**2)**(-1)/0.468043

    else:

        j_sum = np.arrange(1, N_sum)
        return (j**2 + j_star**2)**(-1)/np.sum((j_sum**2 + j_star**2)**(-1))


def B(om, f=f_30):
    """The frequency part of the GM spectrum."""
    return 2.*f/(np.pi*om*np.sqrt(om**2 + f**2))


def E(om, j):
    return B(om)*H(j)*E_0


def F_disp(om, N, j, f=f_30):
    """Displacement spectra."""
    return b**2*N_0*(om**2 - f**2)*E(om, j)/(N*om**2)


def F_vel(om, N, j, f=f_30):
    """Horizontal velocity spectra."""
    return b**2*N_0*N*(om**2 + f**2)*E(om, j)/om**2


def F_eng(om, N, j):
    """Energy per unit mass spectra."""
    return b**2*N_0*N*E(om, j)

#def F_str(om, N, j, f=f_30):
#
#def F_she(om, N, j, f=f_30):
#
#
# case upper('Str')
#  R = (2*pi*kz).^2*(b.^2*N0/N.*(om.^2-f.^2)./om.^2);
# case upper('She')
#  R = (2*pi*kz).^2*(b.^2*N0*N*(om.^2+f.^2)./om.^2);


def m(om, N, j):
    """Convert from frequency space to vertical wavenumber space."""
    return (np.pi/b)*np.sqrt((N**2 - om**2)/(N_0**2 - om**2))*j


def k(om, N, j, f=f_30):
    """Convert from frequency space to horizontal wavenumber space."""
    return (np.pi/b)*np.sqrt((om**2 - f**2)/(N_0**2 - om**2))*j


def Emk(k, m, E_star=E_0, N=N_0, f=f_30, m_star=3*np.pi/b):
    """The GM spectra in k and m space as defined in Cushman-Roisin."""
    num = 3*f*N*E_star*m/m_star
    den = np.pi*(1 + m/m_star)**(2.5) * (N**2 * k**2 + f**2 * m**2)
    return num/den


#% for normalization purposes compute power spectral density of
#% normalized vertical shear for the GM76 model
#%epsilon0=7.8e-10;
#epsilon0=6.73*10^(-10);
#N0=5.2e-3;
#f0=sw_f(30);
#
#N2mean = n_mean^2;
#
#% power spectral density of vertical shear normalised by N for the GM76 model
#E=6.3e-5; % dimensionless energy level
#b=1300; % scale depth of thermocline
#jstar=3; % mode number
#
#%betastar=pi*jstar/b*sqrt(nanmean(N2mean^2))/N0;
#betastar=pi*jstar/b*sqrt(N2mean)/N0;
#
#
#clear phi_u phi_sn
#phi_u = 3*E*b^3*N0^2/(2*jstar*pi) ./(1+kzax/betastar).^2; % power spectral density of horizontal velocity
#phi_sn = kzax.^2.*phi_u/N2mean; % power spectral density of vertical shear normalised by N