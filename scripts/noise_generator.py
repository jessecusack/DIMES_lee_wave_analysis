# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 11:53:40 2016

@author: jc3e13
"""

import numpy as np
import numpy.random as random
import scipy.signal as sig
import matplotlib.pyplot as plt


def noise(N, dx, beta):
    """Generate a noise with a given spectral slope, beta.

    Parameters
    ----------
    N : scalar
        Number of noise points to generate.
    dx : scalar
        Sample spacing, i.e. if time series in seconds, the number of seconds
        between consecutive measurements.
    beta :scalar
        Power spectral slope of the noise. Negative for red noise, positive for
        blue noise.

    Returns
    -------
    y : 1D array
        Noise.

    Notes
    -----
    The function first generates a spectrum with given slope beta, with random
    phase, then performs an inverse FFT to generate the series.
    """

    f = np.fft.fftfreq(N, dx)[1:N/2]
    fNy = (1./(2.*dx))
    b = beta/2.

    mag = (0.1*random.randn(N/2-1)+ 0.5)*f**b
    magNy = np.sign(random.randn())*(0.1*random.randn() + 0.5)*fNy**b

    # Normalise the spectra so the intagrel is 1. There might be a bug here in
    # that the Nyquist freq doesn't want to be multiplied by 2.
    I = 2*np.trapz(np.hstack((mag, magNy)), np.hstack((f, fNy)))
    mag /= I
    magNy /= I

    phase = random.rand(N/2-1)*2.*np.pi

    real = np.zeros(N)
    imag = np.zeros(N)

    real[1:N/2] = mag*np.cos(phase)
    imag[1:N/2] = mag*np.sin(phase)
    real[:N/2:-1] = real[1:N/2]
    imag[:N/2:-1] = -imag[1:N/2]

    # Nyquist frequency has magnitude too!
    real[N/2] = magNy

    Fy = real + 1j*imag

    return np.fft.ifft(Fy).real


dx = 2.5
N = 660
beta = -2.

x = np.arange(0., N*dx, dx)
y = noise(N, dx, beta)

f, Py = sig.periodogram(y, 1./dx)
Py[0] = 0.

fig, axs = plt.subplots(1, 2)
axs[0].loglog(f, Py)
axs[1].plot(x, y)