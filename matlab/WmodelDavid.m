% The idea of this script is to derive absolute vertical velocities. 

% Using the same method as David.

EM1 = load_floats({'4976a'});
%%
clearvars -except EM1
hpidIndx = 10:400;
PLvls = 200:12.5:1300;
dens0 = 1027;

% Am only doing it for profiles 1:300 for the moment because I do not have
% the correct data quality controls in place. (Some profiles at hpid > 300
% do not have enough data to produce the grid and need to be discarded.)
[ig, Pg, Wf] = EM1.gridVar(hpidIndx, 'P', PLvls, 'Wef');
[~, ~, pc] = EM1.gridVar(hpidIndx, 'P', PLvls, 'pc_ctd');
[~, ~, dens] = EM1.gridVar(hpidIndx, 'P', PLvls, 'isdens');

[~, validHpidIndx] = EM1.getPflSubset(hpidIndx);

oddPfls = rem(validHpidIndx, 2) == 1;
evenPfls = rem(validHpidIndx, 2) == 0;
allPfls =  oddPfls | evenPfls;

% Calculate w^2 - will use thsi for drag calculation
% Change sign of w**2 to be negative on the up profiles
Wf2 = Wf.^2;
Wf2(:,evenPfls) = -Wf2(:,evenPfls);

% What profiles to choose? What depth levels to use?
usePfls = allPfls;
useLvls = 1:size(pc,1);

Wf2_ = Wf2(useLvls,usePfls);
pc_ = pc(useLvls,usePfls);
Pg_ = Pg(useLvls,usePfls);
dens_ = dens(useLvls,usePfls);

fprintf(1,'Wref\t\tkappa\t\talpha\t\tk0\t\tscc\n');

% Number of random groups of profiles and number of profiles in each group.
NGrps = 20; NPfls = 50; 

x = zeros(NGrps,4); Wref = zeros(NGrps,1); alpha = zeros(NGrps,1); 
kappa = zeros(NGrps,1); k0 = zeros(NGrps,1); scc = zeros(NGrps,1);
for i = 1:NGrps
    
    [~, idx] = datasample(Wf2_, NPfls, 2);
    
    % Solve equation Ax = B for x using lscov.
    B = reshape(Wf2_(:,idx),[],1);
    A(:,2) = reshape(pc_(:,idx),[],1);
    A(:,3) = reshape(Pg_(:,idx),[],1);
    A(:,4) = reshape(dens_(:,idx),[],1);
    A(:,1) = 1;

    % Exclude very large velocities 
    ix = abs(B) < 0.03;

    % Now find parameters by least squares minimisation of residual 
    [x(i,:), sx] = lscov(A(ix,:),B(ix));

    % Stupid matlab... 
    a = x(i,1); b = x(i,2); c = x(i,3); d = x(i,4);

    Wref(i) = sqrt(d*dens0);
    alpha(i) = b/(dens0*d);
    kappa(i) = -c/(dens0*d);
    k0(i) = -a/b;

    % Given these parameters this is hte modeled value of squared velocity
    Wrel2 = a + b.*pc + c.*Pg + d.*dens;

    % Difference between actual data and modelled data (both are w^2)
    % Form thsi we coudl sovle for water vertical velocity
    err = Wf2 - Wrel2;
    scc(i) = std(reshape(err,[],1));
    fprintf(1,'%3.3f\t\t%6.2e\t%6.2e\t%5.1f\t\t%6.2e\n',...
        Wref(i), kappa(i), alpha(i), k0(i), scc(i));
end