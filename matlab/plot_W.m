function [] = plot_W(float, hpid_range, varargin)
% Input last argument as 'print' if you want to save the files. Check code
% for default save directory and edit as appropriate. 
%
% TODO: Keep figure handles and print at end to simplify code. 
%
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

figure; 
[a, b, c] = float.gridVar(hpid_range, 'z', -1400:10:-50, 'Ww');
[~, hpids, ~] = float.getPflSubset(hpid_range);
pcolor2(a, b, c);
set(gca(), 'XTick', a(1,:))
set(gca(), 'XTickLabel', hpids)
cb = colorbar;
ylabel(cb, 'w_w (m s^{-1})')
if print_flag
    print(gcf, file_type, fullfile(plot_dir, 'grid'))
end

figure;
float.section(hpid_range, 'DIST', 'z', 'Ww', @pcolor2, false)
if print_flag
print(gcf, file_type, fullfile(plot_dir, 'dist_interp'))
end 

% figure;
% float.section(hpid_range, 'MLT_refs', 'z', 'Ww', @pcolor2, false)
% if print_flag
% print(gcf, file_type, fullfile(plot_dir, 'time_interp'))
% end
% 
% float.plot_subset(@plot, hpid_range, 'Ww', 'z')
% if print_flag
% print(gcf, file_type, fullfile(plot_dir, 'multi_profile'))
% end