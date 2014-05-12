% Model parameter file:

mp.savePath = 'W_model_results';
mp.saveFName = 'odd_4976_10-110.mat';

% My idea: list model parameters here. The model reads this file and does
% its thing. 

mp.hpidIndx = 10:110;
mp.PLvls = 150:12.5:1350;
mp.dens0 = 1027;
mp.useOdd = true;
mp.useEven = false;

% Number of random groups of profiles and number of profiles in each group.
mp.NGrps = 50;
mp.NPfls = 10;