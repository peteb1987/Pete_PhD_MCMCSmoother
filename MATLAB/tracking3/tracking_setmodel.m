function model = tracking_setmodel(test)

% Model parameters

% Using a 3D near constant velocity model and a
% bearing-elevation-range-rangerate observation model.

%%%%%%%%%%%%%%%%

% General things
model.K = 100;                  % Number of time points
model.ds = 6;                   % Dimension of the state
model.do = 3;                   % Dimension of the observations

% Parameters
sigx = 1;
sigtheta = ( 0.1*(pi/180) )^2;   % Bearing covariance %0.25;5;
sigphi   = ( 0.1*(pi/180) )^2;   % Elevation covariance
sigr     = 0.01;                  % Range covariance %0.1;100

% Matrixes
T = 1;                          % Sampling period
a = 0.02;                        % Ensures stability
model.A = [exp(-a*T) 0 0 T 0 0; 
           0 exp(-a*T) 0 0 T 0;
           0 0 exp(-a*T) 0 0 T;
           0 0 0 exp(-a*T) 0 0;
           0 0 0 0 exp(-a*T) 0;
           0 0 0 0 0 exp(-a*T)];
model.Q = sigx * ...
    [T^3/3  0      0      T^2/2  0      0    ;
     0      T^3/3  0      0      T^2/2  0    ;
     0      0      T^3/3  0      0      T^2/2;
     T^2/2  0      0      T      0      0    ;
     0      T^2/2  0      0      T      0    ;
     0      0      T^2/2  0      0      T    ];
model.R = diag([sigtheta sigphi sigr]);

% x1 distribution
model.m1 = [-100 50 50 10 0 0]';
model.P1 = diag([10 10 10 1 1 1]);

end