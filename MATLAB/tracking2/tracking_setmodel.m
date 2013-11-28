function model = tracking_setmodel(test)

% Model parameters

% Using a 2D near constant velocity model and a
% bearing-range observation model.

%%%%%%%%%%%%%%%%

% General things
model.K = 250;                  % Number of time points
model.ds = 4;                   % Dimension of the state
model.do = 2;                   % Dimension of the observations

% Parameters
sigx = 1;
sigtheta = ( 5*(pi/180) )^2; % Bearing covariance (0.25/5)
sigr     = 100;                 % Range covariance (0.1/100)

% Matrixes
T = 1;                          % Sampling period
a = 0.02;                        % Ensures stability
model.A = [exp(-a*T) 0 T 0; 
           0 exp(-a*T) 0 T;
           0 0 exp(-a*T) 0;
           0 0 0 exp(-a*T)];
model.Q = sigx * ...
    [T^3/3  0      T^2/2  0    ;
     0      T^3/3  0      T^2/2;
     T^2/2  0      T      0    ;
     0      T^2/2  0      T    ];
model.R = diag([sigtheta sigr]);

% x1 distribution
model.m1 = [-100; 50; 10; 0];
model.P1 = diag([0.0005; 0.0005; 0.001; 0.001]);

end