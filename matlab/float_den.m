function [fden] = float_den(flabel,T,P,ppos)
% function float_den(fnum,T,P,ppos) determines the density of an EM-APEX float
%  Need to know which float, in-situ T and P, and piston position.
% Ballasting gave weight and neutral buoyancy at 1000m at a given T, S and
% piston position. Also compressibility and thermal expansion
% coefficients. Some of these will need verification...
% 

beta = 32e-6;  % all floats use the same temperature coefficient

[pres_bp,temp_bp,sal_bp,pc_bp,mass,alpha,pccpc,addwt] = float_ballast(flabel);

% sw density at ballast point, g/cc
rho_bp = sw_dens(sal_bp,temp_bp,pres_bp) / 1000;

% volume (cc) of float and seawater at ballast point
vbp = mass / rho_bp;

% volume (cc) of float ref'd to ballast point
% with piston at ballast point
vfp = mass * (-alpha) .* (P-pres_bp);
vft = mass * beta .* (T-temp_bp);
vfpc = (ppos - pc_bp) * pccpc;

% float volume at this p,T,pc
vf = vbp + vfp + vft + vfpc;

fden = (mass + addwt)./vf * 1000; % kg/m^3

return
