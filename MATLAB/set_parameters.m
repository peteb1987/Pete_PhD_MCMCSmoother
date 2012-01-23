global params;

%Model
params.dt = 1;
params.N = 100;
params.proc_var = 1;
params.bng_var = (pi/180)^2;    % Try (pi/720)^2; (pi/180)^2;
params.rng_var = 10;           % Try 0.1;        10;
params.x0 = [-100; 50; 10; 0];

% Prior
params.prior_var = diag([0.0005; 0.0005; 0.001; 0.001]);

% Algorithm
params.Np = 100;
params.S = 10;
params.M = 5;                   % Try 5;           50;
