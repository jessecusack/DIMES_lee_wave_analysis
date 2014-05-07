function [p,T,S,pc,M,alpha,pccpc,M_add] = float_ballast(flabel);
  % master ballasting database for EM-APEX floats
  %
  % input (flabel) can be a single string or cell array of M strings, but
  % must contain both float number and deployment string (letter)
  %
  % Returns the following output vectors:
  % p = ballast pressure
  % T = ballast temperature (in-situ)
  % S = ballast salinity
  % M = float mass (used for ballast point)
  % alpha = compressibility
  % M_add = added weight (in water) for this deployment

  % Project list (used)
  %   CBLAST/Frances: 1633,1634,1636
  %   EDDIES: 1632a,1633a,1633b,1635a,1635b,1636a,1636b,1636c (b2)
  %   CLIMODE: 1633c,1636d
  %   PhilEx Exploratory: 1633d,1636e,1636f
  %     PhilEx 2007-2008 IOP: 1633e,1633f,1636g,3304a,3304b,3304c,3305a
  %     PhilEx 2009 IOP: 3305c,4088a,4089a,4090a
  %   Canyon: 3305b
  %   DIMES: 3767a,4086a,4087a

  % Project list (ownership)
  %   CBLAST/Frances (Sanford): 1632,1633,1634,1635,1636
  %   PhilEx (Girton): 3304,3305,4088,4089,4090
  %   DIMES (Girton): 3767,4086,4087
  %   Kerguelen (Phillips): 3760,3761,3762,3764,3950,3951,3952,4051
  %   Gustav (D'Asaro): 3763,3765,3766
  
  X = [...
  	1632 1 1000 5.102 34.878 31 2.72e-6   1.156 27885 0; ...
  	1633 1 1000 5.102 34.878 31 2.72e-6   1.156 28008 0; ...
  	1634 1 1000 5.102 34.878 31 2.72e-6   1.156 27948 0; ...
  	1635 1 1000 5.102 34.878 31 2.72e-6   1.156 27964 0; ...
  	1636 1 1000 5.102 34.878 31 2.72e-6   1.156 27948 0; ...

  	3304 1 1500 3.0   34.6   16 3.67e-6   1.156 26988 0; ...
  	3305 1 1500 3.0   34.6   16 3.67e-6   1.156 27016 0; ...

	3760 1 2000 1.773 34.766 16 3.6394E-6 1.156 27146 45.5; ...
	3761 1 2000 1.915 34.770 16 3.7314E-6 1.156 27201 45.5; ...
	3762 1 2000 2.017 34.766 16 3.7600E-6 1.156 27193 45.5; ...
	3764 1 2000 2.147 34.756 16 3.8227E-6 1.156 27201 45.5; ...
	3950 1 2000 1.773 34.766 16 3.7249E-6 1.156 27226 45.5; ...
	3951 1 2000 1.915 34.770 16 3.6002E-6 1.156 27239 45.5; ...
	3952 1 2000 2.017 34.766 16 3.6316E-6 1.156 27210 45.5; ...
	4051 1 2000 2.147 34.756 16 3.6737E-6 1.156 27171 45.5; ...

  	3763 1 300  13.5  36.0   16 3.80e-6   1.156 27900 30; ...
  	3765 1 300  13.5  36.0   16 3.69e-6   1.156 27900 30; ...
  	3766 1 300  13.5  36.0   16 3.73e-6   1.156 27900 30; ...

  	3767 1 2000 0.400 34.71  16 3.6736e-6 1.156 27178   49; ...
	4086 1 2000 0.400 34.71  16 3.7100e-6 1.156 27177   46; ...
	4087 1 2000 0.400 34.71  16 3.6700e-6 1.156 27175   33; ...
	4088 1 2000 0.400 34.71  16 3.6700e-6 1.156 27154   0; ...
	4089 1 2000 0.400 34.71  16 3.6700e-6 1.156 27136   0; ...
	4090 1 2000 0.400 34.71  16 3.6700e-6 1.156 27197   0; ...
	4812 1 2000 0.400 34.71  16 3.6700e-6 1.156 27197   0; ...
	4813 1 2000 0.400 34.71  16 3.6700e-6 1.156 27197   0; ...
	4814 1 2000 0.400 34.71  16 3.6700e-6 1.156 27197   0; ...
	4815 1 2000 0.400 34.71  16 3.6700e-6 1.156 27197   0; ...
	4976 1 2000 0.400 34.71  16 3.6700e-6 1.156 27197   0; ...
	4977 1 2000 0.400 34.71  16 3.6700e-6 1.156 27197   0; ...
	4594 1 2000 0.400 34.71  16 3.6700e-6 1.156 27197   0; ...
	4595 1 2000 0.400 34.71  16 3.6700e-6 1.156 27197   -50; ...
	4596 1 2000 0.400 34.71  16 3.6700e-6 1.156 27197   0; ...
	4597 1 2000 0.400 34.71  16 3.6700e-6 1.156 27197   -70; ...
	6478 1 2000 0.400 34.71  39 3.8100e-6 0.8 27262.6 0; ...
	6480 1 2000 0.400 34.71  39 3.6700e-6 0.8 27197   0; ...
	6481 1 2000 0.400 34.71  39 3.6700e-6 0.8 27197   0; ...
	6625 1 2000 0.400 34.71  39 3.6700e-6 0.8 27197   0; ...
	6626 1 2000 0.400 34.71  39 3.6700e-6 0.8 27197   0; ...
	];
      
 % Notes: - initial weight adjustments for 3304,3305,4088,4089,4090 were only to account for
 % recovery aids, not to change ballasting. So ACTUAL mass of package is
 % different, but ballast point should apply. (could have interpreted this differently)
 % - Not sure about ballast values for 1632-1636. Going from notes in float_dens.m
 
 % Are the following constants across all floats?
 %
 % (1) Piston CC per count
 % pccpc (or pcvol) = 252/(227-9) = 1.156 from Dana Swift, July 2004
 %   alternatives are 260/239 (May 2003)
 %   or 1.1 (Oct 2008)
 %
 % (2) Thermal expansion coefficient
 %     beta = 32e-6;
 %
 % (3) Fin rotation efficiency factor (certainly decreases with the addition
 % of recovery aids)

    
 floatno_list = X(:,1)';
 depl_list = X(:,2)';
 pres_bp_list = X(:,3)';
 temp_bp_list = X(:,4)';
 sal_bp_list  = X(:,5)';
 piston_bp_list  = X(:,6)';
 alpha_list   = X(:,7)';
 pccpc_list   = X(:,8)'; % piston CC per count
 mass_list    = X(:,9)';
 addwt_list    = X(:,10)';
 
 if ~iscell(flabel),
   flabel = {flabel};
 end
 
 M = length(flabel);
 for ii = 1:M,
   fltid = flabel{ii}; % float/deployment string
   fid = str2num(fltid(1:4)); % float number
   if length(fltid)>4
     depchar = fltid(5); % deployment char
   else
     depchar = 'a';
   end
   switch depchar,
     case 'a',
       deploy = 1;
     case 'b',
       deploy = 2;
     case 'c',
       deploy = 3;
     case 'd',
       deploy = 4;
     case 'e',
       deploy = 5;
     case 'f',
       deploy = 6;
     case 'g',
       deploy = 7;
     otherwise
       error('depchar should be a, b, c, d, e, f, or g');
   end

   iflt = find(floatno_list==fid);
   idepl = find(depl_list(iflt)==deploy);
   if length(iflt)==1,
     iline = iflt;
   elseif ~isempty(depl_list(iflt(idepl)))
     iline = iflt(idepl);
   else
     ddiff = abs(depl_list(iflt)-deploy);
     [md,imd] = min(ddiff);
     iline = iflt(imd(1));
     disp(['Couldn''t find float ' num2str(fid) ' deployment ' ...
	   num2str(deploy) '.']);
     disp(['Instead using float ' num2str(fid) ' deployment ' ...
	   num2str(depl_list(iline)) '.']);
   end
   
   p(ii) = pres_bp_list(iline);
   T(ii) = temp_bp_list(iline);
   S(ii) = sal_bp_list(iline);
   pc(ii) = piston_bp_list(iline);
   M(ii) = mass_list(iline);
   alpha(ii) = alpha_list(iline);
   pccpc(ii) = pccpc_list(iline);
   M_add(ii) = addwt_list(iline);
   
 end
