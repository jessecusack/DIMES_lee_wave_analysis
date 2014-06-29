# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 22:49:23 2013

@author: Jesse

Contains functions and classes for investigating and manipulating EM-APEX float
data.
"""

import numpy as np
import scipy.io as io
from scipy.interpolate import griddata
import gsw
import mapping_tools as mptls
import mat2py as m2p
import pickle
import copy


class Profile(object):
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
              TODO: Make this work properly without exiting program.

          Notes
          -----
          There may be issues with the orientation of the returned array
          because numpy.interp insists on increasing points a certain amount of
          sorting has to be done so that interpolation points are monotonically
          increasing.

          TODO: interp was recently changed in SciPy 0.14 and this sorting may
          no longer be necessary.

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
        t = getattr(self, 'UTC')
        tef = getattr(self, 'UTCef')
        rt = getattr(self, 'rUTC', None)

        equal_size = False

        # If sizes are not equal, check for corresponding time arrays.
        if var_2.size == var_1.size:
            equal_size = True
        else:
            for time in [t, tef, rt]:
                if time is None:
                    continue
                elif time.size == var_1.size:
                    t_1 = time
                elif time.size == var_2.size:
                    t_2 = time

        # Find NaN values.
        nans_var_1 = np.isnan(var_1)
        nans_var_2 = np.isnan(var_2)

        try:

            if equal_size:
                # Both arrays are same length.
                nans = nans_var_1 | nans_var_2
                var_1, var_2 = var_1[~nans], var_2[~nans]
                var_2_sorted, idxs = np.unique(var_2, return_index=True)
                var_1_vals = np.interp(var_2_vals, var_2_sorted, var_1[idxs])

            elif not equal_size:
                nans_1 = nans_var_1 | np.isnan(t_1)
                nans_2 = nans_var_2 | np.isnan(t_2)
                var_1, t_1 = var_1[~nans_1], t_1[~nans_1]
                var_2, t_2 = var_2[~nans_2], t_2[~nans_2]
                # np.unique is necessary to make sure inputs to interp are
                # monotonically increasing!
                var_2_sorted, idxs = np.unique(var_2, return_index=True)
                t_2_interp = np.interp(var_2_vals, var_2_sorted, t_2[idxs])
                t_1_sorted, idxs = np.unique(t_1, return_index=True)
                var_1_vals = np.interp(t_2_interp, t_1_sorted, var_1[idxs])

            else:
                raise RuntimeError('Cannot match time array and/or variable'
                                   ' array sizes.')

        except ValueError:

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


class EMApexFloat(object):
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

        # Load the data!
        for key in data.keys():

            d = np.ndim(data[key])

            if d < 1 or d > 2 or '__' in key:
                print("* Skipping: {}.".format(key))
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

        self.dist_ctd = self.UTC.copy()
        nans = np.isnan(self.dist_ctd)
        for i, (lon, lat, time) in enumerate(zip(lons, lats, times)):
            self.profile_ddist[i] = mptls.lldist(lon, lat)
            # Convert time from days to seconds.
            self.profile_dt[i] = np.diff(time)*86400.

            d = np.array([self.dist[i], self.dist[i] + self.profile_ddist[i]])
            idxs = ~nans[:, i]
            self.dist_ctd[idxs, i] = np.interp(self.UTC[idxs, i], time, d)

        self.dist_ef = self.__regrid('ctd', 'ef', self.dist_ctd)

        # Pythagorian approximation (?) of bearing.
        self.profile_bearing = np.arctan2(self.lon_end - self.lon_start,
                                          self.lat_end - self.lat_start)

        # Convert to m s-1 calculate meridional and zonal velocities.
        self.sub_surf_speed = self.profile_ddist*1000./self.profile_dt
        self.sub_surf_u = self.sub_surf_speed*np.sin(self.profile_bearing)
        self.sub_surf_v = self.sub_surf_speed*np.cos(self.profile_bearing)

        # Absolute velocity.
        self.U_abs = self.U + self.sub_surf_u
        self.V_abs = self.V + self.sub_surf_v

        # Derive some important thermodynamics variables.

        # Depth.
        self.z = gsw.z_from_p(self.P, self.lat_start)
        self.z_ca = gsw.z_from_p(self.P_ca, self.lat_start)
        self.zef = gsw.z_from_p(self.Pef, self.lat_start)

        # Absolute salinity.
        self.SA = gsw.SA_from_SP(self.S, self.P, self.lon_start,
                                 self.lat_start)
        # Conservative temperature.
        self.CT = gsw.CT_from_t(self.SA, self.T, self.P)

        # Potential temperature with respect to 0 dbar.
        self.PT = gsw.pt_from_CT(self.SA, self.CT)

        # In-situ density.
        self.rho = gsw.rho(self.SA, self.CT, self.P)

        # Potential density with respect to 1000 dbar.
        self.rho_1 = gsw.pot_rho_t_exact(self.SA, self.T, self.P, p_ref=1000.)

        # Buoyancy frequency regridded onto ctd grid.
        N2_ca, __ = gsw.Nsquared(self.SA, self.CT, self.P, self.lat_start)
        self.N2 = self.__regrid('ctd_ca', 'ctd', N2_ca)

        # Vertical velocity regridded onto ctd grid.
        dt = 86400.*np.diff(self.UTC, axis=0)  # [s]
        Wz_ca = np.diff(self.z, axis=0)/dt
        self.Wz = self.__regrid('ctd_ca', 'ctd', Wz_ca)

        # Shear calculations.
        dUdz_ca = np.diff(self.U, axis=0)/np.diff(self.zef, axis=0)
        dVdz_ca = np.diff(self.V, axis=0)/np.diff(self.zef, axis=0)
        self.dUdz = self.__regrid('ef_ca', 'ef', dUdz_ca)
        self.dVdz = self.__regrid('ef_ca', 'ef', dVdz_ca)

        # Vertical water velocity.
        self.Wpef = self.Wp.copy()
        del self.Wp

        # Regrid piston position.
        self.ppos_ctd = self.__regrid('ef', 'ctd', self.ppos)

        print("Creating array of half profiles.")

        self.Profiles = np.array([Profile(self, h) for h in self.hpid])

        print("Interpolating some variables onto regular grids.")

        z_vals = np.arange(-1400., 0., 5.)
        self.__r_z_vals = z_vals
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

    def __regrid(self, grid_from, grid_to, v):
        """Take a gridded variable and put it on a ctd/ef grid using time as
        the interpolant.

          Parameters
          ----------
          grid_from : string.
              Grid type which variable v is currently on: 'ctd', 'ef', 'ctd_ca'
              or 'ef_ca'. Where _ca suffic indicates a centered array.
          grid_to : string.
              Grid type that you want that output to be on: 'ctd' or 'ef'.
          v : 2-D numpy.ndarry.
              The variable values that need regridding.

          Returns
          -------
          x : 2-D numpy.ndarray.
              The values of v at the interpolation times.

          Notes
          -----
          This will not work if flattened time is not monotonically increasing.

        """
        if grid_from == grid_to:
            return v

        if grid_from == 'ctd':
            pUTC = self.UTC.flatten(order='F')
        elif grid_from == 'ef':
            pUTC = self.UTCef.flatten(order='F')
        elif grid_from == 'ctd_ca':
            pUTC = ((self.UTC[1:, :] + self.UTC[:-1, :])/2.).flatten(order='F')
        elif grid_from == 'ef_ca':
            pUTC = ((self.UTCef[1:, :] +
                    self.UTCef[:-1, :])/2.).flatten(order='F')
        else:
            raise ValueError("Can only grid from 'ctd', 'ef', 'ctd_ca' or "
                             "'ef_ca'.")

        if grid_to == 'ctd':
            xUTC = self.UTC.flatten(order='F')
            shape = self.UTC.shape
        elif grid_to == 'ef':
            xUTC = self.UTCef.flatten(order='F')
            shape = self.UTCef.shape
        else:
            raise ValueError("Can only grid to 'ctd' or 'ef'.")

        v_flat = v.flatten(order='F')
        nans = np.isnan(pUTC) | np.isnan(v_flat)
        x = xUTC.copy()
        xnans = np.isnan(x)
        x[~xnans] = np.interp(xUTC[~xnans], pUTC[~nans], v_flat[~nans])
        return x.reshape(shape, order='F')

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
                          var_1_vals=None, var_2_vals=None, method='linear'):
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

        __, idxs = self.get_profiles(hpids, ret_idxs=True)

        if var_1_vals is None:
            var_1_arr = getattr(self, var_1_name)[:, idxs]
            start = np.nanmin(var_1_arr)
            stop = np.nanmax(var_1_arr)
            var_1_vals = np.linspace(start, stop, 3*idxs.size)

        if var_2_vals is None:
            var_2_arr = getattr(self, var_2_name)[:, idxs]
            start = np.nanmin(np.nanmin(var_2_arr, axis=0))
            stop = np.nanmax(np.nanmax(var_2_arr, axis=0))
            var_2_vals = np.linspace(start, stop, 400)

        x_grid, y_grid = np.meshgrid(var_1_vals, var_2_vals)

        x = getattr(self, var_1_name)[:, idxs].flatten()
        y = getattr(self, var_2_name)[:, idxs].flatten()
        z = getattr(self, var_3_name)[:, idxs].flatten()

        nans = np.isnan(x) | np.isnan(y) | np.isnan(z)

        x, y, z = x[~nans], y[~nans], z[~nans]

        z_grid = griddata((x, y), z, (x_grid, y_grid), method=method)

        return x_grid, y_grid, z_grid

    def get_timeseries(self, hpids, var_name):
        """TODO: Docstring..."""

        def make_timeseries(t, v):
            times = t[:, idxs].flatten(order='F')
            nnans = ~np.isnan(times)
            times = times[nnans]
            vals = v[:, idxs].flatten(order='F')[nnans]
            times, jdxs = np.unique(times, return_index=True)
            vals = vals[jdxs]
            # Convert to datetime objects.
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

    def apply_w_model(self, fit_info):
        """Give a W_fit_info object or path to pickle of such object."""
        # Initially assume path to fit.
        try:
            with open(fit_info, 'rb') as f:
                print("Unpickling fit info.")
                setattr(self, '__wfi', pickle.load(f))
        except TypeError:
            print('Copying fit info.')
            setattr(self, '__wfi', copy.copy(fit_info))

        wfi = getattr(self, '__wfi')

        # Initialise arrays.
        self.rWs = np.empty_like(self.rUTC)
        self.rWw = np.empty_like(self.rUTC)

        if wfi.profiles == 'all':

            data = [getattr(self, 'r' + data_name) for
                    data_name in wfi.data_names]
            self.rWs = wfi.model_func(wfi.p, data, wfi.fixed)
            print("  Added: rWs.")
            self.rWw = self.rWz - self.rWs
            print("  Added: rWw.")

            # Hack to get Ww on ctd grid.
            if ('ppos' in wfi.data_names and 'P' in wfi.data_names
                and 'rho' in wfi.data_names):

                data = [getattr(self, data_name) for
                        data_name in ['ppos_ctd', 'P', 'rho']]
                self.Ws = wfi.model_func(wfi.p, data, wfi.fixed)
                self.Ww = self.Wz - self.Ws

        elif wfi.profiles == 'updown':

            up = up_down_indices(self.hpid, 'up')
            data = [getattr(self, 'r' + data_name)[:, up] for
                    data_name in wfi.data_names]
            self.rWs[:, up] = wfi.model_func(wfi.p[0], data, wfi.fixed)
            print("  Added: rWs. (ascents)")

            down = up_down_indices(self.hpid, 'down')
            data = [getattr(self, 'r' + data_name)[:, down] for
                    data_name in wfi.data_names]
            self.rWs[:, down] = wfi.model_func(wfi.p[1], data, wfi.fixed)
            print("  Added: rWs. (descents)")

            self.rWw = self.rWz - self.rWs
            print("  Added: rWw.")

        self.update_profiles()

    def apply_strain(self, N2_ref_file):
        """Input the path to file that contains grid of adiabatically levelled
        N2."""

        with open(N2_ref_file) as f:

            N2_ref = pickle.load(f)
            setattr(self, 'N2_ref', N2_ref)
            print("  Added: N2_ref.")
            setattr(self, 'strain_z', (self.N2 - N2_ref)/N2_ref)
            print("  Added: strain_z.")

        self.update_profiles()

        for key in ['N2_ref', 'strain_z']:

            name = 'r' + key
            __, __, var_grid = self.get_interp_grid(self.hpid,
                                                    self.__r_z_vals,
                                                    'z', key)
            setattr(self, name, var_grid)
            print("  Added: {}.".format(name))

        self.update_profiles()


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
