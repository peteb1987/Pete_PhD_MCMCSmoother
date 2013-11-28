clup
dbstop if error

% Get test flag
% test_flag = str2double(getenv('SGE_TASK_ID'));
test_flag = 3;

% Batch test parameters
num_seed = 100;
set_batch_parameters;

% Set test-specific parameters
if test_flag == 1
    params.bng_var = (pi/720)^2;
    params.rng_var = 0.1;
elseif test_flag == 2
    params.bng_var = (pi/36)^2;
    params.rng_var = 0.1;
elseif test_flag == 3
    params.bng_var = (pi/36)^2;
    params.rng_var = 100;
end

% Number of smoothing algorithms
NA = 13; % filter, kitigawa, godsill, 5 MCMC-resampling, 5 MCMC-Fearnhead

% Results arrays
times = cell(NA,1);
unique_pts = cell(NA,1);
mean_pos_rmse = cell(NA,1);
mean_vel_rmse = cell(NA,1);
mean_nees = cell(NA,1);
for alg = 1:NA
    times{alg} = zeros(num_seed,1);
    unique_pts{alg} = zeros(num_seed,params.N);
    mean_pos_rmse{alg} = zeros(num_seed,1);
    mean_vel_rmse{alg} = zeros(num_seed,1);
    mean_nees{alg} = zeros(num_seed,1);
end

% Loop through random seed
for rs = 1:num_seed;
    
    rs
    
    % Set random seed
    s = RandStream('mt19937ar', 'seed', rs);
    RandStream.setDefaultStream(s);
    
    %% Generate some Bearings only tracking data
    [ t, x, y ] = generate_radar_data;
    
    %% Run a PF
    init_pts = num2cell(mvnrnd(params.x0', params.prior_var, params.Np)', 1);
    [pts_array, wts_array, filter_pts] = particle_filter( init_pts, t, y, @tracking_ppsl, @tracking_trans, @tracking_obs, 0.5 );
    
    % Resample final frame particles
    parents = systematic_resample(wts_array{end});
    kiti_pts = pts_array{end}(parents);
    
    % Analyse
    filt_rmse = RMSE(x, filter_pts);
    
    kiti_rmse = RMSE(x, kiti_pts);
    [kiti_Nup, kiti_Nuh] = count_unique_particles(kiti_pts);
    [kiti_nees, kiti_nees_over_time] = NEES(x, kiti_pts);
    
    %Store
    mean_pos_rmse{1}(rs) = filt_rmse.mean_pos;
    mean_vel_rmse{1}(rs) = filt_rmse.mean_vel;
    unique_pts{2}(rs,:) = kiti_Nup;
    mean_pos_rmse{2}(rs) = kiti_rmse.mean_pos;
    mean_vel_rmse{2}(rs) = kiti_rmse.mean_vel;
    mean_nees{2}(rs) = kiti_nees;
    
    %% Run smoothers
    
    for alg = 3:NA
        
        % Reset random seed for each smoother
        s = RandStream('mt19937ar', 'seed', 0);
        RandStream.setDefaultStream(s);
        
        switch alg
            case 3
                
                % Run an ordinary smoother
                tic;
                smooth_pts = backward_sampling_smoother( params.S, t, pts_array, wts_array, @tracking_trans );
                time = toc;
                
            case {4,5,6,7,8}
                
                if alg == 4
                    M = 1;
                elseif alg == 5
                    M = 3;
                elseif alg == 6
                    M = 10;
                elseif alg == 7
                    M = 30;
                elseif alg == 8
                    M = 100;
                else
                    error('Invalid branch')
                end

                tic;
                smooth_pts = mcmc_smoother( params.S, M, t, pts_array, wts_array, @tracking_trans );
                time = toc;
                
            case {9,10,11,12,13}
                
                if alg == 9
                    M = 1;
                elseif alg == 10
                    M = 3;
                elseif alg == 11
                    M = 10;
                elseif alg == 12
                    M = 30;
                elseif alg == 13
                    M = 100;
                else
                    error('Invalid branch')
                end
                
                tic;
                smooth_pts = mcmc_newstate_smoother( params.S, M, t, pts_array, wts_array, y, @tracking_trans, @tracking_obs, @tracking_bidirec_ppsl );
                time = toc;
                
        end
        
        % Analyse and store
        rmse = RMSE(x, smooth_pts);
        [Nup, Nuh] = count_unique_particles(smooth_pts);
        [nees, nees_over_time] = NEES(x, smooth_pts);
        
        % Store
        times{alg}(rs) = time;
        unique_pts{alg}(rs,:) = Nup;
        mean_pos_rmse{alg}(rs) = rmse.mean_pos;
        mean_vel_rmse{alg}(rs) = rmse.mean_vel;
        mean_nees{alg}(rs) = nees;
        
        
    end
    
end

%% Average and store
results = cell(NA,1);
for alg = 1:NA
    results{alg}.mean_pos_rmse = mean(mean_pos_rmse{alg});
    results{alg}.mean_vel_rmse = mean(mean_vel_rmse{alg});
    results{alg}.mean_nees = mean(mean_nees{alg});
    results{alg}.unique_pts = mean(unique_pts{alg}, 1);
    results{alg}.times = mean(times{alg});
end

%%
save(['smoother_test' num2str(test_flag) '.mat'], 'results', 'times', 'unique_pts', 'mean_pos_rmse', 'mean_vel_rmse', 'mean_nees', 'params');