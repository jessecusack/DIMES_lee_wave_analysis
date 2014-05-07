function [u v] = velFromTrack(t, lon, lat)
% USAGE:
%  function [u v] = velFromTrack(t, lon, lat)
%
% DESCRIPTION:
%  Given arrays of time, longitude and latitude this function calculates
%  the zonal and meridional velocities estimated for each point in time. It
%  calculates distances assuming a spherical earth. The equations used are
%  only approximate so that if the input lon and lat positions are far
%  apart the error will be large. 
%
%  Requires: m_map toolbox.
%
% INPUT:
%  t = time in days
%  lon = longitude in degrees, meridian at Greenwich, west of 0 is
%  negative.
%  lat = latitude in degrees. Negative is south of the equator.
%
% OUTPUT:
%  u = zonal velocity component (m s-1)
%  v = meridional velocity component (m s-1)

% Earth radius.
re = 6371e+3;

%dd = m_lldist(lon, lat)*1000;            

dlon = diff(lon);
dlat = diff(lat);
dt = diff(t)*24*60*60;

midlat = (lat(2:end) + lat(1:end-1))./2;
midt = (t(2:end) + t(1:end-1))./2;

dx = re.*dlon.*cos(midlat*pi/180)./360;
dy = re.*dlat./180;

midu = dx./dt;
midv = dy./dt;

u = interp1(midt, midu, t, 'linear', 'extrap');
v = interp1(midt, midv, t, 'linear', 'extrap');