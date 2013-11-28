clup
dbstop if error

% DEFINE RANDOM SEED
rand_seed = 0;

% Set random seed
s = RandStream('mt19937ar', 'seed', rand_seed);
RandStream.setDefaultStream(s);

% Rejection sampling? Set to false - its too slow
rej_samp = true;

% Set parameters
set_parameters;

%% Generate some Bearings only tracking data
[ t, x, y ] = generate_radar_data;

% Plot it
xmax = max(max(abs(x(1:2,:))))+5;
figure(1); hold on; xlim([-xmax, xmax]); ylim([-xmax, xmax]);
plot(x(1,:), x(2,:), 'b', 'linewidth', 3);
[x_obs, y_obs] = pol2cart(y(1,:), y(2,:)); plot(x_obs, y_obs, 'r');
plot(0, 0, 'xk', 'markersize', 10)
legend('Target Position', 'Observations')
% figure(2); plot(t, y);

%% Run a PF
init_pts = num2cell(mvnrnd(params.x0', params.prior_var, params.Np)', 1);
[pts_array, wts_array, filter_pts] = particle_filter( init_pts, t, y, @tracking_ppsl, @tracking_trans, @tracking_obs, 0.5 );

% Resample final frame particles
parents = systematic_resample(wts_array{end});
kiti_pts = pts_array{end}(parents);

% Plot filter
figure(3); hold on; xlim([-xmax, xmax]); ylim([-xmax, xmax]);
plot(x(1,:), x(2,:), 'b', 'linewidth', 3)
cellfun(@(x) plot(x(1,:), x(2,:), 'color', [rand rand rand]), filter_pts)

% Plot Kitigawa smoother
figure(4); hold on; xlim([-xmax, xmax]); ylim([-xmax, xmax]);
plot(x(1,:), x(2,:), 'b', 'linewidth', 3)
cellfun(@(x) plot(x(1,:), x(2,:), 'color', [rand rand rand]), kiti_pts)

%% Run an ordinary smoother
tic;
bs_smooth_pts = backward_sampling_smoother( params.S, t, pts_array, wts_array, @tracking_trans );
bs_time = toc;

% Plot it
figure(5); hold on; xlim([-xmax, xmax]); ylim([-xmax, xmax]);
plot(x(1,:), x(2,:), 'b', 'linewidth', 3)
cellfun(@(x) plot(x(1,:), x(2,:), 'color', [rand rand rand]), bs_smooth_pts)

%% Run a rejection-sampling smoother
if rej_samp
    tic;
    brs_smooth_pts = backward_rejectionsampling_smoother( params.S, t, pts_array, wts_array, @tracking_trans );
    brs_time = toc;
    
    % Plot it
    figure(6); hold on; xlim([-xmax, xmax]); ylim([-xmax, xmax]);
    plot(x(1,:), x(2,:), 'b', 'linewidth', 3)
    cellfun(@(x) plot(x(1,:), x(2,:), 'color', [rand rand rand]), brs_smooth_pts)
end

%% Run the MCMC smoother
tic;
mcmc_smooth_pts = mcmc_smoother( params.S, params.M, t, pts_array, wts_array, @tracking_trans );
mcmc_time = toc;

% Plot it
figure(7); hold on; xlim([-xmax, xmax]); ylim([-xmax, xmax]);
plot(x(1,:), x(2,:), 'b', 'linewidth', 3)
cellfun(@(x) plot(x(1,:), x(2,:), 'color', [rand rand rand]), mcmc_smooth_pts)

%% Run the MCMC new-state smoother
tic;
mcmc_newstate_smooth_pts = mcmc_newstate_smoother( params.S, params.M, t, pts_array, wts_array, y, @tracking_trans, @tracking_obs, @tracking_bidirec_ppsl );
mcmc_ns_time = toc;

% Plot it
figure(8); hold on; xlim([-xmax, xmax]); ylim([-xmax, xmax]);
plot(x(1,:), x(2,:), 'b', 'linewidth', 3)
cellfun(@(x) plot(x(1,:), x(2,:), 'color', [rand rand rand]), mcmc_newstate_smooth_pts)

%% Analysis

filt_rmse = RMSE(x, filter_pts);
[filt_nees, filt_nees_over_time] = NEES(x, filter_pts);

kiti_rmse = RMSE(x, kiti_pts);
[kiti_nees, kiti_nees_over_time] = NEES(x, kiti_pts);
[kiti_Nup, kiti_Nuh] = count_unique_particles(kiti_pts);

bs_rmse = RMSE(x, bs_smooth_pts);
[bs_nees, bs_nees_over_time] = NEES(x, bs_smooth_pts);
[bs_Nup, bs_Nuh] = count_unique_particles(bs_smooth_pts);

if rej_samp
    brs_rmse = RMSE(x, brs_smooth_pts);
    [brs_nees, brs_nees_over_time] = NEES(x, brs_smooth_pts);
    [brs_Nup, brs_Nuh] = count_unique_particles(brs_smooth_pts);
else
    brs_rmse = bs_rmse;
    brs_nees = bs_nees;
    brs_nees_over_time = bs_nees_over_time;
    brs_Nup = bs_Nup;
    brs_Nuh = bs_Nuh;
    brs_time = 0;
end

mcmc_rmse = RMSE(x, mcmc_smooth_pts);
[mcmc_nees, mcmc_nees_over_time] = NEES(x, mcmc_smooth_pts);
[mcmc_Nup, mcmc_Nuh] = count_unique_particles(mcmc_smooth_pts);

mcmc_ns_rmse = RMSE(x, mcmc_newstate_smooth_pts);
[mcmc_ns_nees, mcmc_ns_nees_over_time] = NEES(x, mcmc_newstate_smooth_pts);
[mcmc_ns_Nup, mcmc_ns_Nuh] = count_unique_particles(mcmc_newstate_smooth_pts);

%% Output results
figure(9), clf, hold on
plot(t, filt_rmse.pos, 'r'), plot(t, kiti_rmse.pos, 'g'), plot(t, bs_rmse.pos, 'b'), plot(t, brs_rmse.pos, 'k'), plot(t, mcmc_rmse.pos, 'c'), plot(t, mcmc_ns_rmse.pos, 'm')
legend('filter', 'Kitigawa smoother', 'direct smoother', 'direct rejection smoother', 'MCMC smoother', 'MCMC new-state smoother');
xlabel('time'), ylabel('position error')

figure(10), clf, hold on
plot(t, filt_rmse.vel, 'r'), plot(t, kiti_rmse.vel, 'g'), plot(t, bs_rmse.vel, 'b'), plot(t, brs_rmse.vel, 'k'), plot(t, mcmc_rmse.vel, 'c'), plot(t, mcmc_ns_rmse.vel, 'm')
legend('filter', 'Kitigawa smoother', 'direct smoother', 'direct rejection smoother', 'MCMC smoother', 'MCMC news-state smoother');
xlabel('time'), ylabel('velocity error')

figure(11), clf, hold on
plot(t, filt_nees_over_time, 'r'), plot(t, kiti_nees_over_time, 'g'), plot(t, bs_nees_over_time, 'b'), plot(t, brs_nees_over_time, 'k'), plot(t, mcmc_nees_over_time, 'c'), plot(t, mcmc_ns_nees_over_time, 'm')
legend('filter', 'Kitigawa smoother', 'direct smoother', 'direct rejection smoother', 'MCMC smoother', 'MCMC news-state smoother');
xlabel('time'), ylabel('NEES')

figure(12), clf, hold on
plot(t, kiti_Nup, 'g'), plot(t, bs_Nup, 'b'), plot(t, brs_Nup, 'k'), plot(t, mcmc_Nup, 'c'), plot(t, mcmc_ns_Nup, 'm')
legend('Kitigawa smoother', 'direct smoother', 'direct rejection smoother', 'MCMC smoother', 'MCMC news-state smoother');
xlabel('time'), ylabel('num. of unique particles')

figure(13), clf, hold on
plot(t, kiti_Nuh, 'g'), plot(t, bs_Nuh, 'b'), plot(t, brs_Nuh, 'k'), plot(t, mcmc_Nuh, 'c'), plot(t, mcmc_ns_Nuh, 'm')
legend('Kitigawa smoother', 'direct smoother', 'direct rejection smoother', 'MCMC smoother', 'MCMC news-state smoother');
xlabel('time'), ylabel('num. of unique histories')

fprintf(1, '\n');
fprintf(1, 'RMSEs:  |  Filter  |  KitiSm  |  DireSm  |  DiReSm  |  MCMCSm  |  MCMCnsSm|\n');
fprintf(1, 'Position:  %5.2f   |  %5.2f   |  %5.2f   | %5.2f    | %5.2f    | %5.2f    |\n', filt_rmse.mean_pos, kiti_rmse.mean_pos, bs_rmse.mean_pos, brs_rmse.mean_pos, mcmc_rmse.mean_pos, mcmc_ns_rmse.mean_pos);
fprintf(1, 'Velocity:  %5.2f   |  %5.2f   |  %5.2f   | %5.2f    | %5.2f    | %5.2f    |\n', filt_rmse.mean_vel, kiti_rmse.mean_vel, bs_rmse.mean_vel, brs_rmse.mean_vel, mcmc_rmse.mean_vel, mcmc_ns_rmse.mean_vel);
fprintf(1, '\n');
fprintf(1, 'NEESs:  |  KitiSm  |  DireSm  |  DiReSm  |  MCMCSm  |  MCMCnsSm|\n');
fprintf(1, 'Position:  %5.3f   |  %5.3f   | %5.3f    | %5.3f    | %5.3f    |\n', kiti_nees, bs_nees, brs_nees, mcmc_nees, mcmc_ns_nees);
fprintf(1, '\n');
fprintf(1, 'Direct smoother took:           %f seconds.\n', bs_time);
fprintf(1, 'Rejection smoother took:        %f seconds.\n', brs_time);
fprintf(1, 'MCMC smoother took:             %f seconds.\n', mcmc_time);
fprintf(1, 'MCMC new-state smoother took:   %f seconds.\n', mcmc_ns_time);
fprintf(1, '\n');
