function [Wrel2] = Wrel2_model_lin(p, k, P, rho)
% USAGE:
%  Wrel2_model_lin(p, k, P, rho)
%  
% DESCRIPTION:
%  Calculate relative velocity squared between float and water.
%
% INPUT:
%  p = Array of parameters.
%  k = Piston position.
%  P = Pressure (dbar).
%  rho = Potential density (kg m-3) or its inverse. 
%
% OUTPUT:
%  Wrel2 = Relative water velocity squared(m2 s-2)

Wrel2 = p(1) + p(2).*k + p(3).*P + p(4).*rho;