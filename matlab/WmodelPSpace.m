% Explore the model parameter space.

%%
EM1 = load_floats({'4976a'});
%%
hpidIndx = 1:300;
PLvls = 100:12.5:1450;
dens0 = 1027;

% Am only doing it for profiles 1:300 for the moment because I do not have
% the correct data quality controls in place. (Some profiles at hpid > 300
% do not have enough data to produce the grid and need to be discarded.)
[ig, Pg, Wf] = EM1.gridVar(hpidIndx, 'P', PLvls, 'Wef');
[~, ~, pc] = EM1.gridVar(hpidIndx, 'P', PLvls, 'pc_ctd');
[~, ~, dens] = EM1.gridVar(hpidIndx, 'P', PLvls, 'isdens');
dens = dens.^(-1);

[~, validHpidIndx] = EM1.getPflSubset(hpidIndx);

oddp = rem(validHpidIndx, 2) == 1;
evenp = rem(validHpidIndx, 2) == 0;
allp =  oddp | evenp;

% Calculate w^2 - will use thsi for drag calculation
% Change sign of w**2 to be negative on the up profiles
Wf2 = Wf.^2;
Wf2(:,oddp) = -Wf2(:,oddp);

%%
% Number of random parameters to generate.
N = 100000;
% Parameter ranges/
amin = 0; amax = 200;
bmin = 0; bmax = 1e-3;
cmin = 0; cmax = 1e-3;
dmin = 0; dmax = 1e-3;

a = amin + (amax - amin).*rand(N,1);
b = bmin + (bmax - bmin).*rand(N,1);
c = cmin + (cmax - cmin).*rand(N,1);
d = dmin + (dmax - dmin).*rand(N,1);

err = zeros(N,1);

parfor i = 1:N
    
    Wrel2 = a(i) + b(i).*pc + c(i).*Pg + d(i).*dens;
    err(i) = mean(mean(Wf2(:,allp) - Wrel2(:,allp)));
    
end

