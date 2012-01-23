%% Preliminaries

% Set test flag
test_flag = 2;

% Load data
load(['smoother_test' num2str(test_flag) '.mat']);

% 
t = params.dt:params.dt:params.dt*params.N;

%% Output results

figure(3), clf, hold on
plot(t, results.kiti.unique_pts, 'r'), plot(t, results.bs.unique_pts, 'b')
plot(t, results.mcmc1.unique_pts, 'g'), plot(t, results.mcmc3.unique_pts, 'g'), plot(t, results.mcmc10.unique_pts, 'g'), plot(t, results.mcmc30.unique_pts, 'g'), plot(t, results.mcmc100.unique_pts, 'g')
legend('Kitigawa smoother', 'direct smoother', 'MCMC smoother');
xlabel('time'), ylabel('num. of unique particles')

fprintf(1, '\n');
fprintf(1, 'Mean num of unique pts: |  KitiSm  |  DireSm  |  MCMCS1  |  MCMCS3  |  MCMCS10 |  MCMCS30 | MCMCS100 |\n');
fprintf(1, '                           %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |\n', mean(results.kiti.unique_pts), mean(results.bs.unique_pts), mean(results.mcmc1.unique_pts), mean(results.mcmc3.unique_pts), mean(results.mcmc10.unique_pts), mean(results.mcmc30.unique_pts), mean(results.mcmc100.unique_pts));
fprintf(1, '\n');
fprintf(1, 'RMSEs:  |  Filter  |  KitiSm  |  DireSm  |  MCMCS1  |  MCMCS3  |  MCMCS10 |  MCMCS30 | MCMCS100 |\n');
fprintf(1, 'Position:  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |\n', results.filt.mean_pos_rmse, results.kiti.mean_pos_rmse, results.bs.mean_pos_rmse, results.mcmc1.mean_pos_rmse, results.mcmc3.mean_pos_rmse, results.mcmc10.mean_pos_rmse, results.mcmc30.mean_pos_rmse, results.mcmc100.mean_pos_rmse);
fprintf(1, 'Velocity:  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |\n', results.filt.mean_vel_rmse, results.kiti.mean_vel_rmse, results.bs.mean_vel_rmse, results.mcmc1.mean_vel_rmse, results.mcmc3.mean_vel_rmse, results.mcmc10.mean_vel_rmse, results.mcmc30.mean_vel_rmse, results.mcmc100.mean_vel_rmse);
fprintf(1, '\n');
fprintf(1, 'Times:  |  DireSm  |  MCMCS1  |  MCMCS3  |  MCMCS10 |  MCMCS30 | MCMCS100 |\n');
fprintf(1, '           %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |  %5.2f   |\n', results.bs.time, results.mcmc1.time, results.mcmc3.time, results.mcmc10.time, results.mcmc30.time, results.mcmc100.time);
fprintf(1, '\n');