function [] = plot_W_model_output(output_file)
% Figures are plotted in a directory in ../figures/ of the same name as the
% function argument. FOR SOME REASON UNIX MESSES UP FIGURES SO RUN ON
% WINDOWS.

output = load(output_file);

% Input and output stuff.
[~, name, ~] = fileparts(output_file);
plot_dir = ['../figures/' name];
if ~exist(plot_dir, 'dir')
    mkdir(plot_dir)
end
file_type = '-djpeg';

x = output.x;
% mabcd = output.mabcd;
dens0 = output.mp.dens0;

Wref = sqrt(x(:,4)./dens0);
alpha = -x(:,2).*dens0./x(:,4);
kappa = x(:,3).*dens0./x(:,4);
k0 = -x(:,1)./x(:,2) - x(:,4)./(x(:,2).*dens0);

paramTitles = {'a', 'b', 'c', 'd'};
paramUnits = {'m^2 s^{-2}', 'm^2 s^{-2} pp^{-1}', 'm^2 s^{-2} dbar^{-1}', 'm^2 s^{-2} kg m^{-3}'};

for i = 1:4
    h = figure;
    hist(x(:,i),20)
    title(paramTitles{i})
    xlabel([paramTitles{i} ' (' paramUnits{i} ')'])
    print(h, file_type, fullfile(plot_dir, ['hist_' paramTitles{i}]))
end

for i = 1:4
    for j = 1:i
        if i == j; continue; end 
        h = figure;
        plot(x(:,i), x(:,j), '.')
        xlabel([paramTitles{i} ' (' paramUnits{i} ')'])
        ylabel([paramTitles{j} ' (' paramUnits{j} ')'])
        grid
        print(h, file_type, fullfile(plot_dir, ['scatter_' paramTitles{i} '_vs_' paramTitles{j}]))
    end
end

% Rearrange physical parameters into a matrix.
physTitles = {'Wref', 'kappa', 'alpha', 'k0'};
physSymbs = {'W_{ref}', '\kappa', '\alpha_k', 'k_0'};
physUnits = {'m s^{-1}', 'dbar^{-1}', 'pp^{-1}', 'pp'};
physParam = [Wref kappa alpha k0];

for i = 1:4
    h = figure;
    hist(physParam(:,i),20)
    title(physTitles{i})
    xlabel([physSymbs{i} ' (' physUnits{i} ')'])
    ylabel('Counts')
    print(h, file_type, fullfile(plot_dir, ['hist_' physTitles{i}]))
end

for i = 1:4
    for j = 1:i
        if i == j; continue; end 
        h = figure;
        plot(physParam(:,i), physParam(:,j), '.')
        xlabel([physSymbs{i} ' (' physUnits{i} ')'])
        ylabel([physSymbs{j} ' (' physUnits{j} ')'])
        grid
        print(h, file_type, fullfile(plot_dir, ['scatter_' physTitles{i} '_vs_' physTitles{j}]))
    end
end