function [interped_ele] = interp_bathymetry(lons, lats)
% The input arguments must be between -180 to 180 longitude and
% -90 to 90 latitude. Interpolation is done using interp2. This uses Smith
% and Sandwell bathymetry.
margin = 0.5;
latlim = [min(lats)-margin, max(lats)+margin];
lonlim = [min(lons)-margin, max(lons)+margin];
[ele, elelat, elelon] = mygrid_sand([latlim lonlim], 1);
elelon(elelon > 180) = elelon(elelon > 180) - 360;
[elelonm, elelatm] = meshgrid(elelon, elelat);
interped_ele = interp2(elelonm, elelatm, ele, lons, lats);