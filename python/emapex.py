# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 22:49:23 2013

@author: Jesse

Contains functions and classes for investigating and manipulating EM-APEX float
data.
"""

import numpy as np
import scipy.io as io
import gsw
import mapping_tools as mptls
import mat2py as m2p
#import vertical_velocity_model as vvm
#import matplotlib.pyplot as pl


class Profile:
    """Provide to this class its parent float and half profile number and it
    will extract its data.
    """
    def __init__(self, parent_float, hpid):

        self.hpid = hpid

        # Convert the parent float into a dictionary.
        data = vars(parent_float)
        self.floatID = data['floatID']

        self.update(parent_float)

#        print("Profile {} has been created.".format(hpid))

    def interp(self, var_2_vals, var_2_name, var_1_name):
        """Linear 1-D interpolation of variables stored by a Profile.

          Parameters
          ----------
          var_2_vals : 1-D numpy.ndarray of floats
              The values at which to return interpolant.
          var_2_name : string
              The variable to which var_2_vals correspond.
          var_1_name : string
              The variable to be interpolated.

          Returns
          -------
          var_1_vals : 1-D numpy.ndarray of floats
              The interpolated values of variable 1 at given positions in
              variable 2.

          Raises
          ------
          ValueError
              If array sizes cannot be matched with eachother, the ctd time or
              the ef time.
          RuntimeWarning
              If numpy.interp raises a value error which could be because an
              empty array was passed.
              TODO: Make this work property without exiting program.

          Notes
          -----
          There may be issues with the orientation of the returned array
          because numpy.interp insists on increasing points a certain amount of
          sorting has to be done so that interpolation points are monotonically
          increasing.

          Examples
          --------
          The following example returns interpolated temperatures every 10
          meters from -1000 m to -100 m depth.
          >>> import emapex
          >>> Float = emapex.EMApexFloat('allprofs11.mat',4977)
          >>> P = Float.get_profiles(100)
          >>> z_vals = np.arange(-1000., -100., 10)
          >>> T10 = P.interp(z_vals, 'z', 'T')
        """

        # Shorten some names.
        var_1 = getattr(self, var_1_name)
        var_2 = getattr(self, var_2_name)
        t = self.UTC
        tef = self.UTCef

        # Check array sizes *before removing NaNs*. Gives boolean result.
        equal_size = var_2.size == var_1.size
        var_1_ef_var_2_ctd = ((var_1.size < var_2.size) &
                              (var_1.size == tef.size) &
                              (var_2.size == t.size))
        var_1_ctd_var_2_ef = ((var_1.size > var_2.size) &
                              (var_1.size == t.size) &
                              (var_2.size == tef.size))

        # Find NaN values.
        nans_var_1 = np.isnan(var_1)
        nans_var_2 = np.isnan(var_2)
        nans_t = np.isnan(t)
        nans_tef = np.isnan(tef)

        try:

            if equal_size:
                # Both arrays are same length.
                nans = nans_var_1 | nans_var_2
                var_1, var_2 = var_1[~nans], var_2[~nans]
                var_2_sorted, idxs = np.unique(var_2, return_index=True)
                var_1_vals = np.interp(var_2_vals, var_2_sorted, var_1[idxs])

            elif var_1_ef_var_2_ctd:
                # Variable 1 is of length ef and variable 2 is of length ctd.
                # Remove all possible NaNs.
                nans_ef = nans_var_1 | nans_tef
                nans_ctd = nans_var_2 | nans_t
                var_1, tef = var_1[~nans_ef], tef[~nans_ef]
                var_2, t = var_2[~nans_ctd], t[~nans_ctd]
                # np.unique is necessary to make sure inputs to interp are
                # monotonically increasing!
                var_2_sorted, idxs = np.unique(var_2, return_index=True)
                t_interp = np.interp(var_2_vals, var_2_sorted, t[idxs])
                tef_sorted, idxs = np.unique(tef, return_index=True)
                var_1_vals = np.interp(t_interp, tef_sorted, var_1[idxs])

            elif var_1_ctd_var_2_ef:
                # Variable 1 is of length ctd and variable 2 is of length ef.
                nans_ef = nans_var_2 | nans_tef
                nans_ctd = nans_var_1 | nans_t
                var_1, t = var_1[~nans_ctd], t[~nans_ctd]
                var_2, tef = var_2[~nans_ef], tef[~nans_ef]
                var_2_sorted, idxs = np.unique(var_2, return_index=True)
                t_interp = np.interp(var_2_vals, var_2_sorted, tef[idxs])
                t_sorted, idxs = np.unique(t, return_index=True)
                var_1_vals = np.interp(t_interp, t_sorted, var_1[idxs])

            else:
                raise RuntimeError('Cannot match time array and/or variable'
                                   ' array sizes.')

        except ValueError:
#            print('Warning: Interpolation problem with half profile {}.'
#                  ''.format(self.hpid))
            var_1_vals = np.NaN*np.zeros_like(var_2_vals)

        return var_1_vals

    def update(self, parent_float):
        """Checks for new attributes of parent float and adds them as
        attributes of itself.

        Potentially useful if additional processing is done after float
        initialisation and you want to add this as profile data.
        """

        if parent_float.floatID != self.floatID:
            raise RuntimeError('You are attempting to update this profile '
                               'with the wrong parent float.')

        data = vars(parent_float)

        indx = np.where(data['hpid'] == self.hpid)

        for key in data.keys():

            dat = data[key]
            d = np.ndim(dat)

            if d < 1 or d > 2 or '__' in key:
                continue
            elif d == 1:
                setattr(self, key, dat[indx])
            elif d == 2:
                setattr(self, key, np.squeeze(dat[:, indx]))

        pass


class EMApexFloat:
    """Initialise this class with the path to an 'allprofs##.mat' file where
    ## denotes the last two digits of the year, e.g. 11 for 2011. Also provide
    the ID number of the particular float to extract from the file.

    Notes
    -----

    Some variables have '__' in their name to stop them being transfered into
    profile data e.g. __ddist.

    Some variables have '_ca' in their names and this means that they are on a
    centered array which has length one less than the normal ctd array. These
    variables are not put on to a regular grid because they are normally also
    contained on a ctd array

    TODO: need to do something with Wp now it has been renamed. - deleted it.
    """
    def __init__(self, filepath, floatID):

        self.floatID = floatID
        print(
            "Initialising EM-APEX float: {}\n"
            "Attmpting to load data...".format(floatID)
        )

        # Loaded data is a dictionary.
        data = io.loadmat(filepath, squeeze_me=True)

        isFloat = data.pop('flid') == floatID
        del data['ar']
        self.hpid = data.pop('hpid')[isFloat]
#        flip_indx = up_down_indices(self.hpid, 'up')

        # Figure out some useful times.
        t_ctd = data['UTC'][:, isFloat]
        self.UTC_start = t_ctd[0, :]
        self.UTC_end = np.nanmax(t_ctd, axis=0)

        for key in data.keys():

            d = np.ndim(data[key])

            if d < 1 or d > 2 or '__' in key:
                print "* Skipping: {}.".format(key)
                continue
            elif d == 1:
                setattr(self, key, data[key][isFloat])
            elif d == 2:
                # This flips ascent data first, then adds it as an attribute.
#                setattr(self, key, flip_cols(data[key][:, isFloat], flip_indx))
                setattr(self, key, data[key][:, isFloat])
            else:
                print("* Don't know what to do with {}, skipping".format(key))

            print("  Loaded: {}.".format(key))

        print(
            "All numerical data appears to have been loaded successfully.\n"
            "Interpolating GPS positions and calculating thermodynamic\n"
            "variables."
        )

        # GPS interpolation to the start and end time of each half profile.
        idxs = ~np.isnan(self.lon_gps) & ~np.isnan(self.lat_gps)
        self.lon_start = np.interp(self.UTC_start, self.utc_gps[idxs],
                                   self.lon_gps[idxs])
        self.lat_start = np.interp(self.UTC_start, self.utc_gps[idxs],
                                   self.lat_gps[idxs])
        self.lon_end = np.interp(self.UTC_end, self.utc_gps[idxs],
                                 self.lon_gps[idxs])
        self.lat_end = np.interp(self.UTC_end, self.utc_gps[idxs],
                                 self.lat_gps[idxs])

        # Distance along track from first half profile.
        self.__ddist = mptls.lldist(self.lon_start, self.lat_start)
        self.dist = np.hstack((0., np.cumsum(self.__ddist)))

        # Distances, velocities and speeds of each half profile.
        self.profile_ddist = np.zeros_like(self.lon_start)
        self.profile_dt = np.zeros_like(self.lon_start)
        self.profile_bearing = np.zeros_like(self.lon_start)
        lons = np.zeros((len(self.lon_start), 2))
        lats = lons.copy()
        times = lons.copy()
        lons[:, 0], lons[:, 1] = self.lon_start, self.lon_end
        lats[:, 0], lats[:, 1] = self.lat_start, self.lat_end
        times[:, 0], times[:, 1] = self.UTC_start, self.UTC_end
        self.dist_ctd_data = self.UTC.copy()
        nans = np.isnan(self.dist_ctd_data)
        for i, (lon, lat, time) in enumerate(zip(lons, lats, times)):
            self.profile_ddist[i] = mptls.lldist(lon, lat)
            # Convert time from days to seconds.
            self.profile_dt[i] = np.diff(time)*86400.

            d = np.array([self.dist[i], self.dist[i] + self.profile_ddist[i]])
            idxs = ~nans[:, i]
            self.dist_ctd_data[idxs, i] = np.interp(self.UTC[idxs, i], time, d)

        # Pythagorian approximation (?) of bearing.
        self.profile_bearing = np.arctan2(self.lon_end - self.lon_start,
                                          self.lat_end - self.lat_start)

        # Convert to m s-1 calculate meridional and zonal velocities.
        self.sub_surf_speed = self.profile_ddist*1000./self.profile_dt
        self.sub_surf_u = self.sub_surf_speed*np.sin(self.profile_bearing)
        self.sub_surf_v = self.sub_surf_speed*np.cos(self.profile_bearing)

        # Absolute velocity.
        self.U1_abs = self.U1 + self.sub_surf_u
        self.U2_abs = self.U2 + self.sub_surf_u
        self.V1_abs = self.V1 + self.sub_surf_v
        self.V2_abs = self.V2 + self.sub_surf_v

        # Derive some important thermodynamics variables.
        # Depth.
        self.z = gsw.z_from_p(self.P, self.lat_start)
        self.z_ca = gsw.z_from_p(self.P_ca, self.lat_start)
        # Absolute salinity.
        self.SA = gsw.SA_from_SP(self.S, self.P, self.lon_start,
                                 self.lat_start)
        # Conservative temperature.
        self.CT = gsw.CT_from_t(self.SA, self.T, self.P)
        # Potential temperature with respect to 0 dbar.
        self.PT = gsw.pt_from_CT(self.SA, self.CT)
        # In-situ density.
        self.rho = gsw.rho(self.SA, self.CT, self.P)
        # Buoyancy frequency. (Bizarrely output is a tuple)
#        self.Nsquared = gsw.Nsquared(self.SA, self.CT, self.P,
#                                     self.lat_start)[0]
        # Potential density with respect to 1000 dbar.
        self.rho_1 = gsw.pot_rho_t_exact(self.SA, self.T, self.P, p_ref=1000.)

        # Vertical velocity interpolated back onto a ctd grid.
        dt = 86400.*np.diff(self.UTC, axis=0)  # [s]
        self.Wz_ca = np.diff(self.z, axis=0)/dt
        Wz_ca_flat = self.Wz_ca.flatten(order='F')
        UTC_flat = self.UTC.flatten(order='F')
        UTC_ca_flat = \
            ((self.UTC[1:, :] + self.UTC[:-1, :])/2.).flatten(order='F')
        nnans = ~np.isnan(UTC_ca_flat) & ~np.isnan(Wz_ca_flat)
        Wz_flat = UTC_flat.copy()
        unnans = ~np.isnan(Wz_flat)
        Wz_flat[unnans] = np.interp(UTC_flat[unnans], UTC_ca_flat[nnans],
                                    Wz_ca_flat[nnans])
        self.Wz = Wz_flat.reshape(self.UTC.shape, order='F')

        # Vertical water velocity.
        self.Wpef = self.Wp.copy()
        del self.Wp

        print("Creating array of half profiles.")

        self.Profiles = np.array([Profile(self, h) for h in self.hpid])

        print("Interpolating some variables onto regular grids.")

        z_vals = np.arange(-1400., -40., 6.)
        self_dict = self.__dict__
        for key in self_dict.keys():

            d = np.ndim(self_dict[key])

            if d < 2 or d > 2 or '__' in key or '_ca' in key:
                continue
            elif d == 2:
                name = 'r' + key
                __, __, var_grid = self.get_interp_grid(self.hpid, z_vals,
                                                        'z', key)
                setattr(self, name, var_grid)

            print("  Added: {}.".format(name))

        self.update_profiles()

        print("That appears to have worked.\n")

    def get_profiles(self, hpids, ret_idxs=False):
        """Will return profiles requested. Can also return indices of those
        profiles."""

        if np.ndim(hpids) == 0:
            idx = np.unique(self.hpid.searchsorted(hpids))
            if ret_idxs:
                return self.Profiles[idx[0]], idx
            else:
                return self.Profiles[idx[0]]
        elif np.ndim(hpids) == 1:
            hpids = hpids[(np.min(self.hpid) <= hpids) &
                          (hpids <= np.max(self.hpid))]
            idxs = np.unique(self.hpid.searchsorted(hpids))
            if ret_idxs:
                return self.Profiles[idxs], idxs
            else:
                return self.Profiles[idxs]
        else:
            raise RuntimeError('Check arguments.')

    def update_profiles(self):
        print("Updating half profiles.")
        [profile.update(self) for profile in self.Profiles]

    def get_interp_grid(self, hpids, var_2_vals, var_2_name, var_1_name):
        """Grid data from multiple profiles into a matrix. Linear interpolation
        is performed on, but not across profiles.

          Parameters
          ----------
          hpids : 1-D numpy.ndarray of integers or floats.
              The profile ID numbers at for which to construct grid.
          var_2_vals : 1-D numpy.ndarray of floats
              The values at which to return interpolant.
          var_2_name : string
              The variable to which var_2_vals correspond.
          var_1_name : string
              The variable to be interpolated.

          Returns
          -------
          number_grid : 2-D numpy.ndarray of integers.
              A meshgrid of numbers from 0 to len(hpids). (May be less if some
              hpids are missing.)
          var_2_grid : 2-D numpy.ndarray of floats.
              A meshed grid of var_2_vals.
          var_1_grid : 2-D numpy.ndarray of floats.
              The interpolated grid of variable 1 at given positions of
          variable 2.

          Notes
          -----
          Uses the Profile.interp function.

          Examples
          --------
          The following example returns and interpolated temperature grid.
          >>> import emapex
          >>> Float = emapex.EMApexFloat('allprofs11.mat',4977)
          >>> hpids = np.arange(10,40)
          >>> z_vals = np.arange(-1000., -100., 10)
          >>> T10 = Float.get_interp_grid(hpids, z_vals, 'z', 'T')
        """

        profiles = self.get_profiles(hpids)
        number = np.arange(len(profiles))

        number_grid, var_2_grid = np.meshgrid(number, var_2_vals)
        var_1_grid = np.zeros_like(number_grid, dtype=np.float64)

        for i, profile in enumerate(profiles):
            var_1_grid[:, i] = profile.interp(var_2_vals, var_2_name,
                                              var_1_name)

        return number_grid, var_2_grid, var_1_grid

    def get_griddata_grid(self, hpids, var_1_name, var_2_name, var_3_name,
                          var_1_vals=None, var_2_vals=None):
        """

          Parameters
          ----------

          Returns
          -------

          Notes
          -----

          Examples
          --------

        """

#        profiles = self.get_profiles(hpids)

        if var_1_vals is None:
            pass

        pass

    def get_timeseries(self, hpids, var_name):
        """TODO: Docstring..."""

        def make_timeseries(t, v):
            times = t[:, idxs].flatten(order='F')
            nnans = ~np.isnan(times)
            times = times[nnans]
            vals = v[:, idxs].flatten(order='F')[nnans]
            times, jdxs = np.unique(times, return_index=True)
            vals = vals[jdxs]
            times = m2p.datenum_to_datetime(times)
            return times, vals

        # Shorten some names.
        var = getattr(self, var_name)
        t = self.UTC
        tef = self.UTCef
        rt = self.rUTC
        __, idxs = self.get_profiles(hpids, ret_idxs=True)
        var_1_ef = var.size == tef.size
        var_1_ctd = var.size == t.size
        var_1_r = var.size == rt.size

        if var_1_ef:

            times, vals = make_timeseries(tef, var)

        elif var_1_ctd:

            times, vals = make_timeseries(t, var)

        elif var_1_r:

            times, vals = make_timeseries(rt, var)

        else:
            raise RuntimeError('Cannot match time array and/or variable'
                               ' array sizes.')

        return times, vals

#    def apply_vertical_velocity_model(self, float_velocity_model, params,
#                                      data_var_names):
#        """TODO: Docstring...
#        Must update float with regular grids for this to work.
#        """
#
#        data = [getattr(self, var_name) for var_name in data_var_names]
#
#        self.rWs = float_velocity_model(params, data)
#
#        self.rWw = self.rWz - self.rWs
#
#        self.update_profiles()
#
#    def fit_vertical_velocity_model(self, hpids):
#        """"""
#        reload(vvm)
#        params = vvm.fit_model(self, hpids)
#        self.apply_vertical_velocity_model(vvm.still_water_model_2,
#                                           params, ['rppos', 'rP', 'rrho'])
#        t1, Ws = self.get_timeseries(hpids, 'rWs')
#        __, Wz = self.get_timeseries(hpids, 'rWz')
#        __, Wf = self.get_timeseries(hpids, 'rWf')
#        __, Ww = self.get_timeseries(hpids, 'rWw')
#        pl.figure()
#        pl.plot(t1, Ws, 'r', t1, Wz, 'g', t1, Wf, 'r--', t1, Ww, 'b')


def up_down_indices(hpid_array, up_or_down='up'):
    """Given an array of hpid numbers return the indices of numbers that
    correspond to either ascents or decents.
    up_or_down can either be 'up' or 'down' corresponding to ascent and
    descent. It is 'up' by default.

    """

    # Weirdly np.where in this case returns tuples so the [0] is needed.
    if up_or_down == 'up':
        return np.where(hpid_array % 2 == 0)[0]
    elif up_or_down == 'down':
        return np.where(hpid_array % 2 == 1)[0]
    else:
        raise RuntimeError('Inputs are probably wrong.')


def flip_cols(data, cols=None):
    """Input an array of data. Receive flipped array. If array is two
    dimensional then a list of columns should be provided else the whole matrix
    will be flipped. This is different from the numpy.flipud function because
    any end chunk of data won't be flipped, e.g.:

    [1, 2, 3, nan, nan] -> [3, 2, 1, nan, nan]

    also

    [1, 2, nan, 4, nan, nan] -> [4, nan, 2, 1, nan, nan]

    If data is 2D, assumes that each column contains a list of data otherwise
    please transpose the input.

    Sometimes data sets combine combine columns of different lengths and pad
    out the empty space with nans. This is useful for flipping in that
    situation.

    """
    out_data = data.copy()
    d = np.ndim(data)

    def flip(arr):
        flip_arr = np.flipud(arr)
        nnans = ~np.isnan(flip_arr)
        idx = nnans.searchsorted(True)
#        try:
#            indx = next(i for i, el, in enumerate(flip_arr) if ~np.isnan(el))
#        except StopIteration:
#            return flip_arr
        return np.concatenate((flip_arr[idx:], flip_arr[:idx]))

    if d == 1 and cols is None:

        out_data = flip(data)

    elif d == 2 and cols is not None:

        for col_indx, col in zip(cols, data[:, cols].T):

            out_data[:, col_indx] = flip(col)

    elif d == 2 and cols is None:

        for col_indx, col in enumerate(data.T):

            out_data[:, col_indx] = flip(col)

    else:
        raise RuntimeError('Inputs are probably wrong.')

    return out_data


def what_floats_are_in_here(fname):
    """Finds all unique float ID numbers from a given allprofs##.mat file."""
    fs = io.loadmat(fname, squeeze_me=True, variable_names='flid')['flid']
    return np.unique(fs[~np.isnan(fs)])


if __name__ == "__main__":

    E76 = EMApexFloat('../../data/EM-APEX/allprofs11.mat', 4976)
