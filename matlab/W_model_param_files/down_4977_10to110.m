% Model parameter file:

mp.savePath = 'W_model_results\pdens1';
mp.saveFName = 'down_4977_10to110.mat';

% My idea: list model parameters here. The model reads this file and does
% its thing. 

mp.hpidIndx = 10:110;
mp.PLvls = 150:12.5:1350;
mp.dens0 = 1031;
mp.useOdd = true;
mp.useEven = false;

% Number of random groups of profiles and number of profiles in each group.
mp.NGrps = 50;
mp.NPfls = 20;