# -*- coding: utf-8 -*-
"""
Created on Tue May 20 15:45:36 2014

A place for finescale parameterisation functions.

@author: jc3e13
"""


def adiabatic_level():
    """ """
    pass

#% Calculate a section of vertical strain variance from WOCE SR1b CTD data using
#% the method of Mauritzen et al. (2002)
#
#clear all
#close all
#
#istn_centre=[3:20 23:54]; % choose stations to include in calculation
#stn_overlap=0;
#
#switch_deep=1; % to get rid of uppermost 150 m
#switch_etaz=1; % to skip strain calculation
#
#% set pressure range (in db) of adiabatic levelling
#
#plev=400;
#
#% set wavelength limits of strain spectrum integration in K calculation
#
#lzmin_fixed = 50;
#lzmax_fixed = 320;
#
#load ksection_sr1b97_out lzmin lzmax
#
#if switch_etaz~=1
#
#% load CTD data
#
#load /home/rock3/e035/dobox/raw/drake/sr1_97.mat bin_sal bin_temp bin_press lat lon stn % CTD data
#[c,ia,ib]=intersect(stn,istn_centre);
#salin=bin_sal(:,ia);
#temp=bin_temp(:,ia);
#clear bin_sal bin_temp
#for i=1:length(ia)
#  buf_press(:,i)=bin_press;
#end
#press=buf_press; clear buf_press bin_press
#lon=lon(ia);
#lat=lat(ia);
#stn=stn(ia);
#
#if switch_deep==1
#  icyc=[150/2+1:size(press,1)];
#  salin=salin(icyc,:);
#  temp=temp(icyc,:);
#  press=press(icyc,:);
#end
#
#% calculate N in 2 db data
#
#[n2,dum,p_ave] = sw_bfrq(salin,temp,press,lat);
#
#% apply adiabatic levelling method of Bray and Fofonoff (1981)
#% to calculate a reference N^2 profile at each station
#
#pgrid = p_ave(:,1);
#
#for ii=1:length(istn_centre)
#
#  for jj=1:length(pgrid)
#
#    pmin_lev = max([pgrid(jj)-plev/2 press(1,1)]);
#    pmax_lev = min([pgrid(jj)+plev/2 press(end,1)]);
#
#    icyc = find(press(:,ii)>=pmin_lev & press(:,ii)<=pmax_lev);
#
#    if isempty(icyc)
#
#      n2_bray(jj,ii) = nan;
#
#    else
#
#      pbar = nanmean(press(icyc,ii)); % nominal reference press
#      tbar = nanmean(temp(icyc,ii));
#      sbar = nanmean(salin(icyc,ii));
#      rhobar = sw_pden(sbar,tbar,pbar,pbar);
#
#      sv = 1./sw_pden(salin(icyc,ii),temp(icyc,ii),press(icyc,ii),pbar);
#
#% regress press against de-meaned sv and store coefficients
#
#      order = 1; % use linear fit for now
#
#      [s dum] = polyfit(press(icyc,ii),sv-nanmean(sv),order);
#
#      alpha(:,jj,ii) = s';
#      alpha(2,jj,ii) = alpha(2,jj,ii)+nanmean(sv);
#
#      clear s dum
#
#      g = sw_g(lat(ii),sw_dpth(pbar,lat(ii)));
#
#      n2_bray(jj,ii) = -1e-4 * rhobar^2 * g^2 * alpha(1,jj,ii);
#
#    end
#
#  end
#
#end
#
#% calculate vertical strain
#
#etaz = (n2-n2_bray)./n2_bray;
#
#else
#
#load strainsec_sr1b97_out
#
#end
#
#% station loop
#
#for iloc=1:length(istn_centre)
#
#  ilocmin=max(1,iloc-stn_overlap);
#  ilocmax=min(istn_centre(end),iloc+stn_overlap);
#
#% press loop
#
#  dp_centre=[100:100:4800];
#  dp_overlap=256;
#
#  for idp=1:length(dp_centre)
#
#    dzg=diff(p_ave(1:2,1)); % assumes a regular p_ave grid
#
#    imin=max(1,(dp_centre(idp)-dp_overlap)/dzg+1);  % upper bound of N^2 profile segment
#    imax=min((dp_centre(idp)+dp_overlap)/dzg,size(etaz,1)); % lower bound of N^2 profile segment
#
#% isolate desired depth range and calculate vertical strain spectrum
#
#    fftpt = 256;
#    spec = nan*ones(fftpt/2,length(ilocmin:ilocmax));
#
#    strain = etaz(imin:imax,ilocmin:ilocmax);
#
#    dz = sw_dpth(nanmean(diff(p_ave(imin:imax,iloc))),lat(iloc)); % assumes depth bin is uniform over depth range of transform
#    Fs = 1/dz;
#
#    kzax = (2*pi)*(0:(fftpt/2-1))/fftpt*Fs; % wavenumber axis
#
#    for j=1:size(strain,2)
#
#      straintemp = strain(~isnan(strain(:,j)),j);
#      pts(ilocmin+j-1,idp) = length(straintemp);
#      ratio(ilocmin+j-1,idp)=pts(ilocmin+j-1,idp)/size(strain,1);
#
#      pcthresh=0.5; % define threshold for calculation
#
#      if ratio(ilocmin+j-1,idp) >= pcthresh
#
#% calculate power spectral density of vertical strain
#
#        ffttemp=fft(straintemp,fftpt);
#        spectemp=(ffttemp.*conj(ffttemp))/pts(ilocmin+j-1,idp)/(pi*Fs);
#        spec(:,j)=spectemp(1:fftpt/2); % use only meaningful points
#      end
#    end
#
#    clear meanspec*
#    meanspec=meannan(spec,2)'; % mean spectrum for selected depth range and station group
#
#    stnthres=0.5; % define no. of stations threshold
#
#    [cii,cjj]=find(isfinite(spec));
#    if isempty(cjj) | length(cjj([diff(cjj')~=0 logical(1)]))/size(spec,2) < stnthres
#      meanspec=nan*ones(size(meanspec));
#    end
#
#% integrate vertical strain power spectral density
#
#    kmin_fixed=max(find(kzax<(2*pi/lzmax_fixed)));
#    kmax_fixed=min(find(kzax>(2*pi/lzmin_fixed)));
#%   kmin=max(find(kzax<(2*pi/lzmax(iloc,idp))));
#%   kmax=min(find(kzax>(2*pi/lzmin(iloc,idp))));
#    kmin=max(find(kzax<(2*pi/300)));
#    kmax=min(find(kzax>(2*pi/80)));
#
#    num_fixed=polyarea([kzax(kmin_fixed) kzax(kmin_fixed:kmax_fixed) kzax(kmax_fixed) kzax(kmin_fixed)],[0 meanspec(kmin_fixed:kmax_fixed)' 0 0]);
#    num=polyarea([kzax(kmin) kzax(kmin:kmax) kzax(kmax) kzax(kmin)],[0 meanspec(kmin:kmax)' 0 0]);
#
#% power spectral density of vertical strain for the GM76 model
#
#    epsilon0=7.8e-10;
#    N0=5.24e-3;
#    N2mean = nanmean(nanmean(n2_bray(imin:imax,ilocmin:ilocmax)));
#
#    E=6.3e-5; % dimensionless energy level
#    b=1300; % scale depth of thermocline
#    jstar=3; % mode number
#
#    betastar=pi*jstar/b*sqrt(N2mean)/N0;
#
#    phi_zeta = E*b^3*N0^2/(2*jstar*pi*N2mean) ./(1+kzax/betastar).^2; % power spectral density of vertical displacement
#    phi_eta = kzax.^2.*phi_zeta; % power spectral density of vertical strain
#
#% integrate GM76 power spectral density
#
#    den_fixed=polyarea([kzax(kmin_fixed) kzax(kmin_fixed:kmax_fixed) kzax(kmax_fixed) kzax(kmin_fixed)],[0 phi_eta(kmin_fixed:kmax_fixed) 0 0]);
#    den=polyarea([kzax(kmin) kzax(kmin:kmax) kzax(kmax) kzax(kmin)],[0 phi_eta(kmin:kmax) 0 0]); % not sure why Mauritzen et al. (2002) use factor of 7/3 for this
#
#    strain_spec(idp,iloc)=num; % strain spectrum for shear-to-strain ratio calculation
#    strain_gmspec(idp,iloc)=den; % GM76 strain spectrum for shear-to-strain ratio calculation
#    strain_n2mean(idp,iloc)=N2mean; % N^2 for shear-to-strain ratio calculation
#    norm_strain_var(idp,iloc)=num_fixed^2/den_fixed^2; % strain variance normalized to GM76
#    epsilon=epsilon0*N2mean/N0^2*num_fixed^2/den_fixed^2;
#
#% calculate diapycnal turbulent eddy diffusivity
#
#    gamma=0.2; % mixing efficiency gamma <= 0.2
#
#    strain_K(idp,iloc)=gamma*epsilon/N2mean; % diapycnal turbulent diffusivity
#
#  end
#end
#
#figure
#pcolor(istn_centre,dp_centre,norm_strain_var)
#axis ij
#dum=cat(2,'Norm. strain variance section for stations ',num2str(istn_centre(1)),'-',num2str(istn_centre(end)));
#title(dum)
#xlabel 'Station'
#ylabel 'Depth [m]'
#shading flat
#colorbar
#
#save strainsec_sr1b97_out.mat
