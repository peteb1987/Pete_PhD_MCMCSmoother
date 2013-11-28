global params;

%Model
params.dt = 1;
params.N = 100;
params.proc_var = 1;
params.bng_var = (pi/36)^2;    % Try (pi/720)^2; (pi/180)^2;
params.rng_var = 0.1;           % Try 0.1;        10;
params.x0 = [-100; 50; 10; 0];

% Prior
params.prior_var = diag([0.0005; 0.0005; 0.001; 0.001]);

% Algorithm
params.Np = 1000;
params.S = 100;
params.M = 10;                   % Try 5;           50;
