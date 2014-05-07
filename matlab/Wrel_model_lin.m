function [Wrel] = Wrel_model_lin(p, k, P, rho, hpids)
% USAGE:
%  Wrel_model_lin(p, k, P, rho, hpids)
%  
% DESCRIPTION:
%  Calculate relative velocity between float and water. 
%
% INPUT:
%  p = Array of parameters.
%  k = Piston position.
%  P = Pressure (dbar).
%  rho = Potential density (kg m-3) !!-- > or its inverse. <--!!
%  hpids = hpid number or array of such numbers.
%
% OUTPUT:
%  Wrel = Relative water velocity (m s-1)

if ~ (size(k,2) == size(P, 2) == size(rho, 2) == length(hpids))
    error(['Number of half profile numbers does not match number of' ...
          ' columns of data.'])
end


if length(hpids) > 1

    evenPfls = rem(hpids, 2) == 0;

    Wrel = sqrt(Wrel2_model_lin(p, k, P, rho));

    Wrel(:, evenPfls) = -sqrt(-Wrel2_model_lin(p, k(:, evenPfls), ...
        P(:, evenPfls), rho(:, evenPfls)));

else
      
    if rem(hpids, 2) == 1
        Wrel = sqrt(Wrel2_model_lin(p, k, P, rho));
    elseif rem(hpids, 2) == 0
        Wrel = -sqrt(-Wrel2_model_lin(p, k, P, rho));  
    end
    
end