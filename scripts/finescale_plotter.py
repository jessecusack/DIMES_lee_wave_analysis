# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 12:21:05 2014

@author: jc3e13
"""

hpids = np.arange(1, 60)
__, idxs = E76.get_profiles(hpids, ret_idxs=True)
dists = E76.dist[idxs]

pf.scatter_section(E76, hpids, 'Ww', cmap=plt.get_cmap('bwr'))
plt.clim(-0.1, 0.1)
plt.savefig('../figures/finescale/Ww.png', bbox_inches='tight')

plt.figure()
plt.title('log10 R_pol')
for result, dist in zip(results, dists):
    z_bins, z_mean, EK, R_pol, R_om, epsilon, kappa = result
    d = dist*np.ones_like(z_mean)
    plt.scatter(d, z_mean, c=np.log10(R_pol), edgecolor='none',
                cmap=plt.get_cmap('bwr'))

plt.clim(-1, 1)
plt.colorbar()
plt.ylim(-1400., -200.)
plt.xlim(np.min(dists), np.max(dists))
plt.savefig('../figures/finescale/R_pol.png', bbox_inches='tight')

plt.figure()
plt.title('log10 epsilon')
for result, dist in zip(results, dists):
    z_bins, z_mean, EK, R_pol, R_om, epsilon, kappa = result
    d = dist*np.ones_like(z_mean)
    plt.scatter(d, z_mean, c=np.log10(epsilon), edgecolor='none',
                cmap=plt.get_cmap('jet'))

plt.colorbar()
plt.ylim(-1400., -200.)
plt.xlim(np.min(dists), np.max(dists))
plt.savefig('../figures/finescale/epsilon.png', bbox_inches='tight')

plt.figure()
plt.title('log10 kappa')
for result, dist in zip(results, dists):
    z_bins, z_mean, EK, R_pol, R_om, epsilon, kappa = result
    d = dist*np.ones_like(z_mean)
    plt.scatter(d, z_mean, c=np.log10(kappa), edgecolor='none',
                cmap=plt.get_cmap('jet'))

plt.colorbar()
plt.ylim(-1400., -200.)
plt.xlim(np.min(dists), np.max(dists))
plt.savefig('../figures/finescale/kappa.png', bbox_inches='tight')

plt.figure()
plt.title('R_om')
for result, dist in zip(results, dists):
    z_bins, z_mean, EK, R_pol, R_om, epsilon, kappa = result
    d = dist*np.ones_like(z_mean)
    plt.scatter(d, z_mean, c=R_om, edgecolor='none',
                cmap=plt.get_cmap('jet'))

plt.colorbar()
plt.ylim(-1400., -200.)
plt.xlim(np.min(dists), np.max(dists))
plt.savefig('../figures/finescale/R_om.png', bbox_inches='tight')

plt.figure()
plt.title('log10 EK')
for result, dist in zip(results, dists):
    z_bins, z_mean, EK, R_pol, R_om, epsilon, kappa = result
    d = dist*np.ones_like(z_mean)
    plt.scatter(d, z_mean, c=np.log10(EK), edgecolor='none',
                cmap=plt.get_cmap('jet'))

plt.colorbar()
plt.ylim(-1400., -200.)
plt.xlim(np.min(dists), np.max(dists))
plt.savefig('../figures/finescale/EK.png', bbox_inches='tight')
