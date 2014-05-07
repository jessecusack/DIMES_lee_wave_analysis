function [] = plot_everything(float, hpid_range, varargin)
% Does what it says!

plot_tracks(float, hpid_range, varargin{:})

plot_pdens_displacement(float, hpid_range, varargin{:})

plot_W(float, hpid_range, varargin{:})

subplot_velocity_profiles(float, hpid_range, varargin{:})

wave_analysis(float, hpid_range, varargin{:})