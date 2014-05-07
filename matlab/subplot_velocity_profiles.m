function [] = subplot_velocity_profiles(float, hpid_range, varargin)

% Plots of vertical velocity with depth.

% Plots of horizontal velocity too.

N = length(hpid_range);
[pfls, hpids, ~] = float.getPflSubset(hpid_range);

% Input and output stuff.
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

% Velocity component vertical profiles.
figure;
i = 1;
for pfl = pfls'
    subplot(1,N,i)

    pfl.plot(@line, 'Ww', 'z', 'Color', 'b');
    ax1 = gca;
    axis([-0.2, 0.2 -1600, 0])
    set(ax1, 'XColor', 'b')
    xlabel('w_{water} (m s^{-1})', 'Color', 'b')
    text(-0.15, -100, num2str(hpids(i)))
    if i == 1
        ylabel('Depth (m)')
    else
        set(ax1, 'YTickLabel', [])
    end
%     title(num2str(hpids(i)))

    ax2 = axes('Position', get(ax1, 'Position'), ...
               'XAxisLocation', 'top', ...
               'YAxisLocation', 'right', ...
               'Color', 'none', ...
               'XColor', 'r', ...
               'YTickLabel', []);
           
    pfl.plot(@line, 'U1', 'z', 'Color', 'r', ...
        'LineStyle', '-');
    pfl.plot(@line, 'V1', 'z', 'Color', 'r', ...
        'LineStyle', '--');
    axis([-0.5, 0.5 -1600, 0])
    xlabel('u, v (m s^{-1})', 'Color', 'r')
    i = i + 1;
end

if print_flag
end
