function [] = wave_analysis(EMF, hpid_range, varargin)
% Basic docstring: Input an EMApexFloat and an hpid range and this function
% will calculate the dominant wavelengths present in the the velocity and
% maybe the potential density. 

pfls = EMF.getPflSubset(hpid_range);
N = length(pfls);

% Interpolate velocities onto regular 4m grid.
D = 4; % Sampling distance (m).
Ks = 1/D;
Z = -1400:D:-40;
WNS = 1./(0:50:1400);
i = 1;

figure
for pfl = pfls'

    subplot(N,1,i)
    W = pfl.interp_var('Ww', 'z', Z);
    h = spectrum.welch;
    pwelch(W)%, 'Fs', Ks, 'ConfLevel', 0.95);%, 'Frequencies', WNS);
    %plot(Hpsd.Frequencies, Hpsd.Data)
    i = i + 1;
    
end