function [] = plot_tracks(float, hpid_range, varargin)
% Plot EM-Apex track on bathymetry and save.

plot_dir = '../figures/';
print_flag = false;
file_type = '-djpeg';
N = length(hpid_range);

if ~isempty(varargin)
    if strcmpi(varargin{1}, 'print')
        print_flag = true;
    end
end

[~, ~, indxBool] = float.getPflSubset(hpid_range);
LATS = float.LATS(indxBool); LONS = float.LONS(indxBool);
h = float.plot_track(hpid_range, true);
if print_flag
    print(h, file_type, fullfile(plot_dir, [ float.ID '_track']))
end

% Bathymetry along track:
ele = interp_bathymetry(LONS, LATS);

h = figure;
plot(ele)
xlim([1 N])
x_tick_labels = get(gca, 'XTick');
new_x_tick_labels = hpid_range(ismember(1:N, x_tick_labels));
set(gca, 'XTickLabel', new_x_tick_labels)
if print_flag
    print(h, file_type, fullfile(plot_dir, [ float.ID '_track_bathy']))
end