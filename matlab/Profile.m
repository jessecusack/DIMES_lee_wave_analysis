classdef Profile < handle
% The Profile class will hold all the data associated with a single
% profile and maybe some simple plotting functions.
%
% TOOLBOX DEPENDENCIES: 
%   gsw -- TEOS-10
%
% Variables in new EM-APEX Matlab velocity files
% John Dunlap, dunlap@apl.washington.edu
% August 5, 2009
% 
% A typical file name is "ema-3763a-0065-vel.mat" where,
% 	3763a -- float serial number with a suffix to indicate run number.
% 	0065  -- "hpid", half-profile ID, odd numbers descend, even numbers 
%   ascend.
% 
% On CTD grid:
% 	ctd_mlt -- UTC of CTD samples (Matlab datenum)
% 	Pctd -- pressure (dbar)
% 	T -- temperature (degrees C)
% 	S -- salinity (PSU)
% 
% Note: The CTD data is better quality ascending than descending because
% the CTD inlet is on the top.
% 
% On Processed electric field (EFP) grid:
% 	efp_mlt -- UTC of EFP estimates (Matlab datenum)
% 	Pef   -- pressure (dbar)
% 	U1,U2 -- east component of velocity from electrode sets 1,2 (m/s)
% 	V1,V2 -- north component of velocity from electrode sets 1,2 (m/s)
% 	Verr1 -- standard error of U1,V1 (m/s)
% 	Verr2 -- standard error of U2,V2 (m/s)
% 	Wef   -- vertical velocity derived from Pctd (m/s)
% 	RotP  -- rotation period (s)
% 	e1mean, e2mean -- mean value of EF (microvolts)
% 	e1sdev, e2sdev -- standard deviation of residuals of EF fit 
%                     (microvolts)
% 
% Note: When RotP becomes more than about EmaProcessNvals / 2 the
% velocity processing becomes noisy.  Also, the variables Verr1 and
% Verr2 can be used to screen out noisy U1,V1 and U2,V2.
% 
% GPS data as sampled just before phoning home with data.
% 	MLT_GPS -- UTC of GPS sample (Matlab datenum)
% 	NSAT    -- number of satellites used
% 	LAT     -- latitude of this profile (degrees)
% 	LON     -- longitude of this profile (degrees)
% 	STAT    -- status from GPGGA NMEA message, 1: OK, 0: no fix
% 	HDOP    -- horizontal dilution of precision
% 
% Typically GPS data is sampled once every several profiles so the following
% is interpolated position:
% 	lat     -- latitude interpolated from other profiles (degrees)
% 	lon     -- longitude interpolated from other profiles (degrees)
% 
% More:
% 	MLT_ref -- UTC of beginning of half-profile. (Matlab datenum)
% 	           For ascents the reference is at the deepest point.
% 	magvar  -- magnetic variation at lat/lon (degrees)
% 	fh,fz   -- earth's horizontal and vertical magnetic field (nT)
% 	
% 	hpid    -- half-profile number:  odd numbers descend, even numbers 
%              ascend.
    properties
        % Basic properties
        EmaProcessNvals     
        HDOP
        LAT
        LON
        MLT_GPS
        MLT_ref
        NSAT
        Pctd
        Pef
        RotP
        S
        STAT
        T
        U1
        U2
        V1
        V1woW
        V2
        V2woW
        Verr1
        Verr2
        Wef
        alpha1
        alpha2
        ctd_mlt
        e1mean
        e1sdev
        e2mean
        e2sdev
        efp_mlt
        esep
        fh
        fz
        hpid
        lat
        lon
        magvar
        mlt_surface
        pc_ctd
        pc_efp
        sfv1
        sfv2
        sfw
        
        % Derived properties
        P           % Sea pressure (dbar).
        CT          % Conservative temperature (deg C).
        SA          % Absolute salinity (g kg-1).
        pdens       % Potential density relative to surface (kg m-3).
        isdens      % In situ density (kg m-3).
        z           % Height (m).
        V1abs       % Absolute northward velocity (m s-1).
        V2abs       % Absolute northward velocity (m s-1).
        U1abs       % Absolute eastward velocity (m s-1).
        U2abs       % Absolute eastward velocity (m s-1).
        Wrel        % Relative water velocity (modelled) (m s-1).
        Ww          % Absolute water velocity (modelled) (m s-1).

    end % properties
    
    methods
        function obj = Profile(S)
            % Takes a EM-APEX data structure using the load function and
            % copies the data into the class.
            for fn = fieldnames(S)'
                obj.(fn{1}) = S.(fn{1})';
            end % for
            % Change the sign of the effective velocity so that it is in a
            % coordinate system where z = 0 at the sea surface and
            % decreases going down.
            obj.Wef = -obj.Wef;

            obj.P = obj.Pctd; %- 10.1325;
            % TODO: The lat lon coordinates from the individual half
            % profile files may be incorrect. Better to interpolate lat lon
            % position from GPS data file.
            obj.SA = gsw_SA_from_SP(obj.S, obj.P, obj.lon, obj.lat);
            obj.CT = gsw_CT_from_t(obj.SA, obj.T, obj.P);
            obj.pdens = gsw_rho(obj.SA, obj.CT, 1000);
            obj.isdens = gsw_rho(obj.SA, obj.CT, obj.P);
            obj.z = gsw_z_from_p(obj.P, obj.lat);
            
        end % function 
        
        function [h] = plot(obj, func, propx, propy, varargin)
            % USAGE:
            %  plot(func, propx, propy, varargin)
            %  
            % DESCRIPTION:
            %  This isn't just a plot function! The idea is that it applies
            %  a function to two data properties of the profile. The user
            %  needs to figure out what are appropriate properties to use
            %  with the function. Additional function arguments can also be
            %  passed i.e. keyword value pairs. The function must be able
            %  to return a handle.
            %
            %  If chosen properties are different sizes then some
            %  interpolation is attempted but this is not fully developed
            %  and may fail. TODO: combine with interp_var.
            %
            % INPUT:
            %  func = function handle, e.g. @plot or @line
            %  propx = Profile property and first argument of func e.g. 'T'.
            %  propy = a second Profile property and second argument of func.
            %  varargin = any additional inputs to func.
            %
            % OUTPUT:
            %  h = handle returned from func.
            %
            % EXAMPLE:
            %  % Plot a temperature profile:
            %  Profile.plot(@plot, 'T', 'z', 'Color', 'r')
            
            x = obj.(propx);
            y = obj.(propy);
            
            if length(x) == length(y)
                h = func(x, y, varargin{:});
            elseif length(x) < length(y) && ...
                    length(x) == length(obj.efp_mlt) && ...
                    length(y) == length(obj.ctd_mlt)
                % interpolate y variable onto times of measurements.
                Y = interp1(obj.ctd_mlt, y, obj.efp_mlt);
                h = func(x, Y, varargin{:});
            else
                error('Have you chosen appropriate arguments?')
            end % end if
                        
        end % function
        
        function [] = plot_TD(obj)
            % Use the SA - CT plotting function from the GSW toolbox.
            figure;
            gsw_SA_CT_plot(obj.SA, obj.CT)
            
        end % function
        
        function [] = calc_vert_vel(obj, params)%, model)
            % The idea here is that the parameters are an array that can by
            % understood by the function 'model'. The input 'model' should
            % be a function handle to the model function. The model
            % function should accept the array of parameters and also a
            % profile object. The model function then uses attribute of the
            % profile object (i.e. pressure, density) to calculate the
            % water vertical velocity which is returned here.
            %
            % Nice idea, but simpler to just use one model hard-coded.
            
            try
                pdenef = obj.interp_var('pdens', 'efp_mlt', obj.efp_mlt);
            catch
                pdenef = nan*ones(size(obj.efp_mlt));
            end % try
            
            obj.Wrel = Wrel_model_lin(params, obj.pc_efp, obj.Pef, ... 
                                          pdenef.^(-1), obj.hpid);
            obj.Ww = obj.Wef - obj.Wrel;
            
        end % function
        
        function [] = calc_abs_vel(obj, Udi, Vdi)
            % The input arguments are depth averaged eastward and northward
            % velocities. They are simply added to the relative velocites
            % U1, V1, U2 and V2 to generate the absolute velocites. 
            
            obj.U1abs = obj.U1 + Udi;
            obj.U2abs = obj.U2 + Udi;
            obj.V1abs = obj.V1 + Vdi;
            obj.V2abs = obj.V2 + Vdi;
            
        end % function
        
        function [var1vals] = interp_var(obj, var1, var2, var2vals)
            % Return the values of variable 1 (e.g. 'U1') at the same
            % position as values of variable 2 (e.g. 'z') using linear 
            % interpolation.
            % It is up to the user to put in sensible values for variable
            % 2.
            % This function has only strictly been tested for the case
            % where both variables are measured at the same time as the ctd
            % array. 
            
            % Simplify input.
            v1 = obj.(var1);
            v2 = obj.(var2);
            tctd = obj.ctd_mlt;
            tefp = obj.efp_mlt;
            % Make sure strictly monotonic by sorting and removing repeated
            % values. 
            [v2sorted indxs] = unique(v2);
            
            if length(v2) == length(v1) % Both arrays are the same length.
                
                var1vals = interp1(v2sorted, v1(indxs), var2vals);
                
            elseif length(v1) < length(v2) && ... % Array 1 is efp.
                    length(v1) == length(tefp) && ... % Array 2 is ctd.
                    length(v2) == length(tctd)
                
                t = interp1(v2sorted, tctd(indxs), var2vals);
                var1vals = interp1(tefp, v1, t);                
                
            elseif length(v1) > length(v2) && ... % Array 1 is ctd.
                    length(v2) == length(tefp) && ... % Array 2 is efp.
                    length(v1) == length(tctd)
                
                t = interp1(v2sorted, tefp(indxs), var2vals);
                var1vals = interp1(tctd, v1, t);                
                
            else
                error('Have you chosen appropriate variables?')
            end
                
        end % function
        
    end % methods
    
end % classdef