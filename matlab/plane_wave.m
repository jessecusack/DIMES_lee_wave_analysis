function [y] = plane_wave(p, x)
y = p(1) +  p(2)*sin(p(3)*(x + p(4)));