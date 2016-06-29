# -*- coding: utf-8 -*-
"""
Created on Sun Nov 22 19:56:26 2015

@author: jesse
"""


def smooth_buoyancy(Float, P_bin_width=100., save_dir='../../data/EM-APEX'):
    """Smooth buoyancy frequency and save to file."""

    Pg = getattr(Float, 'P')
    SAg = getattr(Float, 'SA')
    Tg = getattr(Float, 'T')
    lats = getattr(Float, 'lat_start')

    N2_ref = np.NaN*Pg.copy()

    for i, (P, SA, T, lat) in enumerate(zip(Pg.T, SAg.T, Tg.T, lats)):
        print("hpid: {}".format(Float.hpid[i]))
        N2_ref[:, i] = adiabatic_level(P, SA, T, lat, P_bin_width)

    save_name = "{:g}_N2_ref_{:g}dbar.p".format(Float.floatID, P_bin_width)
    file_path = os.path.join(save_dir, save_name)

    pickle.dump(N2_ref, open(file_path, 'wb'))


def smooth_density(Float, z_bin_width=100., save_dir='../../data/EM-APEX'):
    """Smooth potential density and save to a file."""

    srho_1 = np.nan*Float.rho_1.copy()

    for i in xrange(len(Float.hpid)):
        print("hpid: {}".format(Float.hpid[i]))
        srho_1[:, i] = wdw.moving_polynomial_smooth(
            Float.z[:, i], Float.rho_1[:, i], width=100., deg=1.)

    save_name = "srho_{:g}_{:g}mbin.p".format(Float.floatID, z_bin_width)
    file_path = os.path.join(save_dir, save_name)

    pickle.dump(srho_1, open(file_path, 'wb'))


def analyse_profile(Pfl, params=default_params):
    """ """

    if params['zmin'] is None:
        params['zmin'] = np.nanmin(Pfl.z)

    # First remove NaN values and interpolate variables onto a regular grid.
    dz = params['dz']
    z = np.arange(params['zmin'], params['zmax']+dz, dz)
    U = Pfl.interp(z, 'zef', 'U_abs')
    V = Pfl.interp(z, 'zef', 'V_abs')
    dUdz = Pfl.interp(z, 'zef', 'dUdz')
    dVdz = Pfl.interp(z, 'zef', 'dVdz')
    strain = Pfl.interp(z, 'z', 'strain_z')
    N2_ref = Pfl.interp(z, 'z', 'N2_ref')
    lat = (Pfl.lat_start + Pfl.lat_end)/2.

    return analyse(z, U, V, dUdz, dVdz, strain, N2_ref, lat, params)


def analyse_float(Float, hpids, params=default_params):
    """ """
    # Nothing special for now. It doesn't even work.
    __, idxs = Float.get_profiles(hpids, ret_idxs=True)
    return [analyse_profile(Pfl, params) for Pfl in Float.Profiles[idxs]]
    

def w_scales_float(Float, hpids, c=0.1, eff=0.2, lc=30.):

    __, idxs = Float.get_profiles(hpids, ret_idxs=True)

    w = Float.r_Ww[:, idxs]
    z = Float.r_z[:, 0]
    N2 = Float.r_N2_ref[:, idxs]

    dz = z[0] - z[1]

    epsilon = np.zeros_like(w)
    kappa = np.zeros_like(w)

    for i, (w_row, N2_row) in enumerate(zip(w.T, N2.T)):
        epsilon[:, i], kappa[:, i] = w_scales(w_row, z, N2_row, dz, c, eff, lc)

    return epsilon, kappa
