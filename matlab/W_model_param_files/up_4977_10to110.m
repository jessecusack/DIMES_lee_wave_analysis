% Model parameter file:

saveFName = 'W_model_results\pdens1\up_4977_10to110';

% My idea: list model parameters here. The model reads this file and does
% its thing. 

hpidIndx = 10:110;
PLvls = 150:12.5:1350;
dens0 = 1031;
useOdd = false;
useEven = true;

% Number of random groups of profiles and number of profiles in each group.
NGrps = 50;
NPfls = 20;