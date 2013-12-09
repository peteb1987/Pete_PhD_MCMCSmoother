clup;

test_name = 'model_case_1'
num_tests = 100;
num_algs = 16;

mean_pos_rmse = zeros(num_algs,num_tests);
mean_vel_rmse = zeros(num_algs,num_tests);
mean_tnees = zeros(num_algs,num_tests);
mean_nus = zeros(num_algs,num_tests);
mean_rt = zeros(num_algs,num_tests);

for tt = num_tests:-1:1
    
    try
        load(['3D_batch_tests/' test_name '/bs_test_' num2str(tt) '.mat']);
        
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

inds = {1, 2, [3:7], [8:12], [13:16]};
colours = 'kgcbr';

figure, hold on
for aa = 1:length(inds)
    
    plot(mean_mean_rt(inds{aa}), mean_mean_pos_rmse(inds{aa}), 'x-', 'color', colours(aa));
    
end