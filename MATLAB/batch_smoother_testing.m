clup
dbstop if error

% Get test flag
% test_flag = str2double(getenv('SGE_TASK_ID'));
test_flag = 2;

% Batch test parameters
num_seed = 10;
set_parameters;

% Results arrays
bs_times = zeros(num_seed,1);
mcmc1_times = zeros(num_seed,1);
mcmc3_times = zeros(num_seed,1);
mcmc10_times = zeros(num_seed,1);
mcmc30_times = zeros(num_seed,1);
mcmc100_times = zeros(num_seed,1);

kiti_unique_pts = zeros(num_seed,params.N);
bs_unique_pts = zeros(num_seed,params.N);
mcmc1_unique_pts = zeros(num_seed,params.N);
mcmc3_unique_pts = zeros(num_seed,params.N);
mcmc10_unique_pts = zeros(num_seed,params.N);
mcmc30_unique_pts = zeros(num_seed,params.N);
mcmc100_unique_pts = zeros(num_seed,params.N);

filt_mean_pos_rmse = zeros(num_seed,1);
kiti_mean_pos_rmse = zeros(num_seed,1);
bs_mean_pos_rmse = zeros(num_seed,1);
mcmc1_mean_pos_rmse = zeros(num_seed,1);
mcmc3_mean_pos_rmse = zeros(num_seed,1);
mcmc10_mean_pos_rmse = zeros(num_seed,1);
mcmc30_mean_pos_rmse = zeros(num_seed,1);
mcmc100_mean_pos_rmse = zeros(num_seed,1);

filt_mean_vel_rmse = zeros(num_seed,1);
kiti_mean_vel_rmse = zeros(num_seed,1);
bs_mean_vel_rmse = zeros(num_seed,1);
mcmc1_mean_vel_rmse = zeros(num_seed,1);
mcmc3_mean_vel_rmse = zeros(num_seed,1);
mcmc10_mean_vel_rmse = zeros(num_seed,1);
mcmc30_mean_vel_rmse = zeros(num_seed,1);
mcmc100_mean_vel_rmse = zeros(num_seed,1);


% Loop through random seend
for rs = 1:num_seed;
    
    % Set random seed
    s = RandStream('mt19937ar', 'seed', rs);
    RandStream.setDefaultStream(s);
    
    % Set parameters
    set_parameters;
    
    %% Generate some Bearings only tracking data
    [ t, x, y ] = generate_radar_data;
    
    %% Run a PF
    init_pts = num2cell(mvnrnd(params.x0', params.prior_var, params.Np)', 1);
    [pts_array, wts_array, filter_pts] = particle_filter( init_pts, t, y, @tracking_ppsl, @tracking_trans, @tracking_obs, 0.5 );
    
    % Resample final frame particles
    parents = systematic_resample(wts_array{end});
    kiti_pts = pts_array{end}(parents);
    
    % Analyse and store
    filt_rmse = RMSE(x, filter_pts);
    
    kiti_rmse = RMSE(x, kiti_pts);
    [kiti_Nup, kiti_Nuh] = count_unique_particles(kiti_pts);
    
    kiti_unique_pts(rs,:) = kiti_Nup;
    filt_mean_pos_rmse = filt_rmse.mean_pos;
    kiti_mean_pos_rmse = kiti_rmse.mean_pos;
    filt_mean_vel_rmse = filt_rmse.mean_vel;
    kiti_mean_vel_rmse = kiti_rmse.mean_vel;
    
    %% Run an ordinary smoother
    tic;
    bs_smooth_pts = backard_sampling_smoother( params.S, t, pts_array, wts_array, @tracking_trans );
    bs_time = toc;
    
    % Analyse and store
    bs_rmse = RMSE(x, bs_smooth_pts);
    [bs_Nup, bs_Nuh] = count_unique_particles(bs_smooth_pts);
    
    bs_times(rs) = bs_time;
    bs_unique_pts(rs,:) = bs_Nup;
    bs_mean_pos_rmse = bs_rmse.mean_pos;
    bs_mean_vel_rmse = bs_rmse.mean_vel;
    
    %% Run the MCMC smoothers
    
    M = 1;
    % Run
    tic;
    mcmc_smooth_pts = mcmc_smoother( params.S, M, t, pts_array, wts_array, @tracking_trans );
    mcmc_time = toc;
    % Analyse and store
    mcmc_rmse = RMSE(x, mcmc_smooth_pts);
    [mcmc_Nup, mcmc_Nuh] = count_unique_particles(mcmc_smooth_pts);
    mcmc1_times(rs) = mcmc_time;
    mcmc1_unique_pts(rs,:) = mcmc_Nup;
    mcmc1_mean_pos_rmse = mcmc_rmse.mean_pos;
    mcmc1_mean_vel_rmse = mcmc_rmse.mean_vel;
    
    M = 3;
    % Run
    tic;
    mcmc_smooth_pts = mcmc_smoother( params.S, M, t, pts_array, wts_array, @tracking_trans );
    mcmc_time = toc;
    % Analyse and store
    mcmc_rmse = RMSE(x, mcmc_smooth_pts);
    [mcmc_Nup, mcmc_Nuh] = count_unique_particles(mcmc_smooth_pts);
    mcmc3_times(rs) = mcmc_time;
    mcmc3_unique_pts(rs,:) = mcmc_Nup;
    mcmc3_mean_pos_rmse = mcmc_rmse.mean_pos;
    mcmc3_mean_vel_rmse = mcmc_rmse.mean_vel;
    
    M = 10;
    % Run
    tic;
    mcmc_smooth_pts = mcmc_smoother( params.S, M, t, pts_array, wts_array, @tracking_trans );
    mcmc_time = toc;
    % Analyse and store
    mcmc_rmse = RMSE(x, mcmc_smooth_pts);
    [mcmc_Nup, mcmc_Nuh] = count_unique_particles(mcmc_smooth_pts);
    mcmc10_times(rs) = mcmc_time;
    mcmc10_unique_pts(rs,:) = mcmc_Nup;
    mcmc10_mean_pos_rmse = mcmc_rmse.mean_pos;
    mcmc10_mean_vel_rmse = mcmc_rmse.mean_vel;
    
    M = 30;
    % Run
    tic;
    mcmc_smooth_pts = mcmc_smoother( params.S, M, t, pts_array, wts_array, @tracking_trans );
    mcmc_time = toc;
    % Analyse and store
    mcmc_rmse = RMSE(x, mcmc_smooth_pts);
    [mcmc_Nup, mcmc_Nuh] = count_unique_particles(mcmc_smooth_pts);
    mcmc30_times(rs) = mcmc_time;
    mcmc30_unique_pts(rs,:) = mcmc_Nup;
    mcmc30_mean_pos_rmse = mcmc_rmse.mean_pos;
    mcmc30_mean_vel_rmse = mcmc_rmse.mean_vel;
    
    M = 100;
    % Run
    tic;
    mcmc_smooth_pts = mcmc_smoother( params.S, M, t, pts_array, wts_array, @tracking_trans );
    mcmc_time = toc;
    % Analyse and store
    mcmc_rmse = RMSE(x, mcmc_smooth_pts);
    [mcmc_Nup, mcmc_Nuh] = count_unique_particles(mcmc_smooth_pts);
    mcmc100_times(rs) = mcmc_time;
    mcmc100_unique_pts(rs,:) = mcmc_Nup;
    mcmc100_mean_pos_rmse = mcmc_rmse.mean_pos;
    mcmc100_mean_vel_rmse = mcmc_rmse.mean_vel;
    
end

% Average and store

results.filt.mean_pos_rmse = mean(filt_mean_pos_rmse, 1);
results.kiti.mean_pos_rmse = mean(kiti_mean_pos_rmse, 1);
results.bs.mean_pos_rmse = mean(bs_mean_pos_rmse, 1);
results.mcmc1.mean_pos_rmse = mean(mcmc1_mean_pos_rmse, 1);
results.mcmc3.mean_pos_rmse = mean(mcmc3_mean_pos_rmse, 1);
results.mcmc10.mean_pos_rmse = mean(mcmc10_mean_pos_rmse, 1);
results.mcmc30.mean_pos_rmse = mean(mcmc30_mean_pos_rmse, 1);
results.mcmc100.mean_pos_rmse = mean(mcmc100_mean_pos_rmse, 1);

results.filt.mean_vel_rmse = mean(filt_mean_vel_rmse, 1);
results.kiti.mean_vel_rmse = mean(kiti_mean_vel_rmse, 1);
results.bs.mean_vel_rmse = mean(bs_mean_vel_rmse, 1);
results.mcmc1.mean_vel_rmse = mean(mcmc1_mean_vel_rmse, 1);
results.mcmc3.mean_vel_rmse = mean(mcmc3_mean_vel_rmse, 1);
results.mcmc10.mean_vel_rmse = mean(mcmc10_mean_vel_rmse, 1);
results.mcmc30.mean_vel_rmse = mean(mcmc30_mean_vel_rmse, 1);
results.mcmc100.mean_vel_rmse = mean(mcmc100_mean_vel_rmse, 1);

results.kiti.unique_pts = mean(kiti_unique_pts, 1);
results.bs.unique_pts = mean(bs_unique_pts, 1);
results.mcmc1.unique_pts = mean(mcmc1_unique_pts, 1);
results.mcmc3.unique_pts = mean(mcmc3_unique_pts, 1);
results.mcmc10.unique_pts = mean(mcmc10_unique_pts, 1);
results.mcmc30.unique_pts = mean(mcmc30_unique_pts, 1);
results.mcmc100.unique_pts = mean(mcmc100_unique_pts, 1);

results.bs.time = mean(bs_times);
results.mcmc1.time = mean(mcmc1_times);
results.mcmc3.time = mean(mcmc3_times);
results.mcmc10.time = mean(mcmc10_times);
results.mcmc30.time = mean(mcmc30_times);
results.mcmc100.time = mean(mcmc100_times);

save(['smoother_test' num2str(test_flag) '.mat'], 'results', 'params');