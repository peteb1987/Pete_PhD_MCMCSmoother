clup;

test_name = 'model_case_3'
num_tests = 100;
num_algs = 16;
num_traj = 100;
result_type = 'nus'; % 'pos_rmse' 'vel_rmse' 'tnees' 'nus'

mean_pos_rmse = zeros(num_algs,num_tests);
mean_vel_rmse = zeros(num_algs,num_tests);
mean_tnees = zeros(num_algs,num_tests);
mean_nus = zeros(num_algs,num_tests);
mean_rt = zeros(num_algs,num_tests);

for tt = num_tests:-1:1
    
    try
        load(['batch_tests/' test_name '/bs_test_' num2str(tt) '.mat']);
        
        for aa = 1:num_algs
            mean_pos_rmse(aa,tt) = mean(results.pos_rmse{aa});
            mean_vel_rmse(aa,tt) = mean(results.vel_rmse{aa});
            mean_tnees(aa,tt) = mean(results.tnees{aa});
            mean_nus(aa,tt) = mean(results.nus{aa});
            mean_rt(aa,tt) = mean(results.rt(aa));
        end
        
    catch
        
        mean_pos_rmse(:,tt) = [];
        mean_vel_rmse(:,tt) = [];
        mean_tnees(:,tt) = [];
        mean_nus(:,tt) = [];
        mean_rt(:,tt) = [];
        
    end
    
end

mean_mean_pos_rmse = mean(mean_pos_rmse,2);
mean_mean_vel_rmse = mean(mean_vel_rmse,2);
mean_mean_tnees = mean(mean_tnees,2);
mean_mean_nus = mean(mean_nus,2);
mean_mean_rt = mean(mean_rt,2);

%% Plot some graphs

% inds = {1, 2, [3:7], [8:12], [13:16]};
% colours = 'kgcbr';

inds = {1, 2, [8:12], [13:16]};
colours = 'krbc';

figure, hold on
for aa = 1:length(inds)
    
    switch result_type
        case 'pos_rmse'
            plot(mean_mean_rt(inds{aa})/num_traj, mean_mean_pos_rmse(inds{aa}), 'x-', 'color', colours(aa), 'markersize', 10);
        case 'vel_rmse'
            plot(mean_mean_rt(inds{aa})/num_traj, mean_mean_vel_rmse(inds{aa}), 'x-', 'color', colours(aa), 'markersize', 10);
        case 'tnees'
            plot(mean_mean_rt(inds{aa})/num_traj, mean_mean_tnees(inds{aa}), 'x-', 'color', colours(aa), 'markersize', 10);
        case 'nus'
            plot(mean_mean_rt(inds{aa})/num_traj, mean_mean_nus(inds{aa}), 'x-', 'color', colours(aa), 'markersize', 10);
    end
    
end

xl = xlim;
yl = ylim;
plot([0 0], yl, ':k');
xlim([-10 xl(2)]);
ylim(yl);

xlabel('Average time per smoothing particle (s)');
switch result_type
    case 'pos_rmse'
        ylabel('Average Position RMSE')
        legend('Filter', 'Direct Backward Simulation Smoother', 'MH Backward Simulation Smoother', 'MH Refreshed Backward Simulation Smoother')
    case 'vel_rmse'
        ylabel('Average Velocity RMSE')
        legend('Filter', 'Direct Backward Simulation Smoother', 'MH Backward Simulation Smoother', 'MH Refreshed Backward Simulation Smoother')
    case 'tnees'
        ylabel('Average TNEES')
        legend('Filter', 'Direct Backward Simulation Smoother', 'MH Backward Simulation Smoother', 'MH Refreshed Backward Simulation Smoother')
    case 'nus'
        ylabel('Average Number of Unique States')
        legend('Filter', 'Direct Backward Simulation Smoother', 'MH Backward Simulation Smoother', 'MH Refreshed Backward Simulation Smoother', 'Location','SouthEast')
end

export_pdf(gcf, [test_name '_results_' result_type '.pdf'], 6, 4);