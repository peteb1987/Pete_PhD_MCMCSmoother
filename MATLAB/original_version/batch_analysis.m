%% Preliminaries

close all

% Set test flag
test_flag = 3;

% Load data
load(['smoother_test' num2str(test_flag) '.mat']);

% 
NA = 13;
t = params.dt:params.dt:params.dt*params.N;
col = {[1 0 0], [1 .5 .5], [0 1 0], [.5 1 .5], [0 0 1], [.5 .5 1]};

%% Output results

% Unique particles over time
figure(1), hold on
to_plot = [2,3,4,6,8,10];
for ii = 1:length(to_plot)
    alg = to_plot(ii);
    plot(t, results{alg}.unique_pts, 'color', col{ii});
end
legend('filter-smoother', 'direct resampling smoother', 'MCMC resampling smoother (M=1)', 'MCMC resampling smoother (M=10)', 'MCMC sampling smoother (M=1)', 'MCMC sampling smoother (M=10)');
xlabel('time'), ylabel('num. of unique particles')

figure(2), hold on
bar( cellfun(@(x) mean(x.unique_pts), results) );

figure(3), hold on
bar( cellfun(@(x) x.mean_pos_rmse, results) );

figure(4), hold on
bar( cellfun(@(x) x.mean_vel_rmse, results) );

figure(5), hold on
bar( cellfun(@(x) x.times, results) );

figure(6), hold on
inds = 2; plot( cellfun(@(x) x.times, results(inds)), cellfun(@(x) x.mean_pos_rmse, results(inds)), '*k' )
inds = 3; plot( cellfun(@(x) x.times, results(inds)), cellfun(@(x) x.mean_pos_rmse, results(inds)), '*r' )
inds = 4:8; plot( cellfun(@(x) x.times, results(inds)), cellfun(@(x) x.mean_pos_rmse, results(inds)), '*-b' )
inds = 9:13; plot( cellfun(@(x) x.times, results(inds)), cellfun(@(x) x.mean_pos_rmse, results(inds)), '*-g' )
legend('FS', 'DBRS', 'MCMC-BRS (M=1-100)', 'MCMC-BSS (M=1-100)');
xlabel('Running Time'); ylabel('Position RMSE');
% saveas(6, 'smoother_comparison_posRMSE_time.pdf', 'pdf');

% fprintf(1, '\n');
% fprintf(1, 'Mean num of unique pts: |  KitiSm  |  DireSm  |  MCMCS1  |  MCMCS3  |  MCMCS10 |  MCMCS30 | MCMCS100 |\n');
% fprintf(1, '                           %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |\n', mean(results.kiti.unique_pts), mean(results.bs.unique_pts), mean(results.mcmc1.unique_pts), mean(results.mcmc3.unique_pts), mean(results.mcmc10.unique_pts), mean(results.mcmc30.unique_pts), mean(results.mcmc100.unique_pts));
% fprintf(1, '\n');
% fprintf(1, 'RMSEs:  |  Filter  |  KitiSm  |  DireSm  |  MCMCS1  |  MCMCS3  |  MCMCS10 |  MCMCS30 | MCMCS100 |\n');
% fprintf(1, 'Position:  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |\n', results.filt.mean_pos_rmse, results.kiti.mean_pos_rmse, results.bs.mean_pos_rmse, results.mcmc1.mean_pos_rmse, results.mcmc3.mean_pos_rmse, results.mcmc10.mean_pos_rmse, results.mcmc30.mean_pos_rmse, results.mcmc100.mean_pos_rmse);
% fprintf(1, 'Velocity:  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |\n', results.filt.mean_vel_rmse, results.kiti.mean_vel_rmse, results.bs.mean_vel_rmse, results.mcmc1.mean_vel_rmse, results.mcmc3.mean_vel_rmse, results.mcmc10.mean_vel_rmse, results.mcmc30.mean_vel_rmse, results.mcmc100.mean_vel_rmse);
% fprintf(1, '\n');
% fprintf(1, 'Times:  |  DireSm  |  MCMCS1  |  MCMCS3  |  MCMCS10 |  MCMCS30 | MCMCS100 |\n');
% fprintf(1, '           %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |\n', results.bs.time, results.mcmc1.time, results.mcmc3.time, results.mcmc10.time, results.mcmc30.time, results.mcmc100.time);
% fprintf(1, '\n');