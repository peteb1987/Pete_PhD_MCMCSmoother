function [ t, x, y ] = generate_radar_data
%GENERATE_RADAR_DATA Generates some radar-style tracking data for inference
%algorithm testing

global params;

[A, Q] = build_AQ(params.dt, params.proc_var);
R = [params.bng_var 0; 0 params.rng_var];
N = params.N;
x0 = params.x0;

% Initialise arrays
x = zeros(4, N);
y = zeros(2, N);
t = params.dt:params.dt:params.dt*N;

last_x = x0;

% Loop
for n = 1:N
    
    % Simulate new state
    x(:,n) = mvnrnd((A*last_x)', Q)';
    last_x = x(:,n);
    
    % Simulate new observation
    [b, r] = cart2pol(x(1,n), x(2,n));
    y(:,n) = mvnrnd([b, r], R)';
    
end

end

