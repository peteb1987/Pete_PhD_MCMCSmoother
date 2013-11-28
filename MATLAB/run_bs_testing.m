% Base script for backward simulation smoothing comparisons on a tracking model

%% Preliminaries

if ~exist('test', 'var') || ~isfield(test,'flag_batch') || (~test.flag_batch)
    
    clup
    dbstop if error
    %     dbstop if warning
    
    % Set flag to non-batch
    test.flag_batch = false;
    
    %%% SETTINGS %%%
    
    % DEFINE RANDOM SEED
    rand_seed = 0;
    
    % Set display options
    display.text = true;
    display.plot_during = false;
    display.plot_after = true;
    
    % Tests to run
    test.samp_types = [0 1 2 3 4];
    test.chain_lengths = [NaN NaN 10 10 10];
    
end

%% Setup
fprintf('   Random seed: %u.\n', rand_seed);

% Function handles for the model
addpath('tracking2');
fh.setmodel = @tracking_setmodel;
fh.setalgo = @tracking_setalgo;
fh.generatedata = @tracking_generatedata;
fh.transition = @tracking_transition;
fh.observation = @tracking_observation;
fh.stateprior = @tracking_stateprior;
fh.ekfproposal = @tracking_ekfproposal;
fh.ukfproposal = @tracking_ukfproposal;

% Set random seed
rng(rand_seed);

% Set model parameters
[model] = feval(fh.setmodel, test);

% Set algorithm parameters
[algo] = feval(fh.setalgo, test);

% Generate data
[time, state, observ] = feval(fh.generatedata, model);

%% Filtering

% Run the particle filter
ppsl_type = 3; % 1,2,3 = bootstrap,EKF,UKF resp.
[pf, pf_diagnostics] = particle_filter(display, algo, model, fh, observ, ppsl_type, state);

%% Smoothing

num_alg = length(test.samp_types);
ps = cell(num_alg,1);
ps_diagnostics = cell(num_alg,1);
rt = zeros(num_alg,1);

% Run the smoothers
for aa = 1:num_alg
    type = test.samp_types(aa);
    algo.M = test.chain_lengths(aa);
    tic;
    if type == 4
        [ps{aa}, ps_diagnostics{aa}] = refreshed_backward_simulation_smoother(display, algo, model, fh, pf, observ, state);
    else
        [ps{aa}, ps_diagnostics{aa}] = backward_simulation_smoother(display, algo, model, fh, pf, observ, type, state);
    end
    rt(aa) = toc;
end

%% Analysis
pos_rmse = cell(num_alg,1);
vel_rmse = cell(num_alg,1);
tnees = cell(num_alg,1);
nus = cell(num_alg,1);
for aa = 1:num_alg
    [pos_rmse{aa}, vel_rmse{aa}, tnees{aa}, nus{aa}] = particle_smoother_analysis(model, state, ps{aa});
end

%% Display
if display.plot_after
    close all
    
    colours = 'kbgrm';
    
    if model.ds == 4
        
        [ox, oy] = pol2cart(observ(1,:), observ(2,:));
        figure, hold on
        plot(state(1,:), state(2,:), 'k--');
        plot(ox, oy, 'r-*')
        
        for aa = 1:num_alg
            figure, hold on
            plot(state(1,:), state(2,:), 'k--');
            for ii = 1:length(ps{aa}), plot(ps{aa}(ii).traj_state(1,:), ps{aa}(ii).traj_state(2,:), colours(aa)); end
        end
        
    elseif model.ds == 6
        
        [ox, oy, oz] = sph2cart(observ(1,:), observ(2,:), observ(3,:));
        figure, hold on
        plot3(state(1,:), state(2,:), state(3,:));
        plot3(ox, oy, oz, 'r-*')
        
        for aa = 1:num_alg
            figure, hold on
            plot3(state(1,:), state(2,:), state(3,:));
            for ii = 1:length(ps{aa}), plot3(ps{aa}(ii).traj_state(1,:), ps{aa}(ii).traj_state(2,:), ps{aa}(ii).traj_state(3,:), colours(aa)); end
        end
        
    end
    
    figure, hold on
    for aa = 1:num_alg
        plot(pos_rmse{aa}, colours(aa));
    end
    
    figure, hold on
    for aa = 1:num_alg
        plot(vel_rmse{aa}, colours(aa));
    end
    
    figure, hold on
    for aa = 1:num_alg
        plot(tnees{aa}, colours(aa));
    end
    
    figure, hold on
    for aa = 1:num_alg
        plot(nus{aa}, colours(aa));
    end
    
end
