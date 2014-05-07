function [h] = pcolor2(varargin)
h = pcolor(varargin{:});
set(h, 'EdgeColor', 'none')