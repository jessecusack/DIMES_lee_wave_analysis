[E76, E77] = load_floats();
upE77 = load('W_model_results\pdens1\up_4977_10to110.mat');
downE77 = load('W_model_results\pdens1\down_4977_10to110.mat');
downE76 = load('W_model_results\pdens1\down_4976_10to110.mat');
upE76 = load('W_model_results\pdens1\up_4976_10to110.mat');
E76.calc_vert_vel(upE76.mabcd, downE76.mabcd)
E77.calc_vert_vel(upE77.mabcd, downE77.mabcd)