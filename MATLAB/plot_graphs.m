clup
dbstop if error

% Function handles for the model
addpath('tracking3');
fh.setmodel = @tracking_setmodel;
fh.setalgo = @tracking_setalgo;
fh.generatedata = @tracking_generatedata;
fh.transition = @tracking_transition;
fh.observation = @tracking_observation;
fh.stateprior = @tracking_stateprior;
fh.ekfproposal = @tracking_ekfproposal;
fh.ukfproposal = @tracking_ukfproposal;

%%% SETTINGS %%%

% DEFINE RANDOM SEED
rand_seed = 1;

test_case = 3;

% Set model parameters
[model] = feval(fh.setmodel, []);

switch test_case
    case 1
        sigtheta = ( 0.1*(pi/180) )^2;  % Bearing covariance
        sigphi   = ( 0.1*(pi/180) )^2;  % Elevation covariance
        sigr     = 0.01;                % Range covariance
    case 2
        sigtheta = ( 0.1*(pi/180) )^2;  % Bearing covariance
        sigphi   = ( 0.1*(pi/180) )^2;  % Elevation covariance
        sigr     = 100;                 % Range covariance
    case 3
        sigtheta = ( 5*(pi/180) )^2;    % Bearing covariance
        sigphi   = ( 5*(pi/180) )^2;    % Elevation covariance
        sigr     = 0.01;                % Range covariance
    otherwise
        error('Invalid test case');
end
model.R = diag([sigtheta sigphi sigr]);


%% Setup
fprintf('   Random seed: %u.\n', rand_seed);

% Set random seed
rng(rand_seed);

% Generate data
[time, state, observ] = feval(fh.generatedata, model);

close all

colours = 'kbgrm';

if model.ds == 4
    
    [ox, oy] = pol2cart(observ(1,:), observ(2,:));
    figure, hold on
    plot(state(1,:), state(2,:), 'k--');
    plot(ox, oy, 'r:*')
    plot(0, 0, 'xk', 'markersize', 3);
    
elseif model.ds == 6
    
    [ox, oy, oz] = sph2cart(observ(1,:), observ(2,:), observ(3,:));
    figure, hold on
    plot3(state(1,:), state(2,:), state(3,:), 'linewidth', 2);
    plot3(ox, oy, oz, 'r:*')
    plot3(0, 0, 0, 'xk', 'markersize', 5);
    
    view(45,30);
    
end

export_pdf(gcf, ['trajectory_case_' num2str(test_case) '.pdf']);