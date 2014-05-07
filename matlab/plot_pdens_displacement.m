function [] = plot_pdens_displacement(float, hpid_range, varargin)
% Plot the vertical displacement of potential density surfaces. 

print_flag = false;
file_type = '-djpeg';

if ~isempty(varargin)
    if strcmpi(varargin{1}, 'print')
        plot_dir = fullfile('../figures/', ['W_' float.ID]);
        if ~exist(plot_dir, 'dir')
            mkdir(plot_dir)
        end
        print_flag = true;
    end
end

% Vertical position of potential density surfaces.
pdens_contours = 1031.5:0.1:1032.2;
[ig, pdens, Z] = float.gridVar(hpid_range, 'pdens', pdens_contours, 'z');
h = figure;
plot(ig', Z')
xlabel('hpid')
ylabel('z (m)')
xlim([ig(1,1) ig(1,end)])
x_tick_labels = get(gca, 'XTick');
new_x_tick_labels = hpid_range(ismember(ig(1,:), x_tick_labels));
set(gca, 'XTickLabel', new_x_tick_labels)
legend(num2str(pdens(:, 1), '%.1f'))

if print_flag
    print(h, file_type, fullfile(plot_dir, [ float.ID '_pdens_displacement']))
end