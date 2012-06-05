global params;

%Model
params.dt = 1;
params.N = 100;
params.proc_var = 1;
params.x0 = [-100; 50; 10; 0];
% params.bng_var    Set for each test
% params.rng_var    Set for each test

% Prior
params.prior_var = diag([0.0005; 0.0005; 0.001; 0.001]);

% Algorithm
params.Np = 100;
params.S = 100;
% params.M          Set of each test
