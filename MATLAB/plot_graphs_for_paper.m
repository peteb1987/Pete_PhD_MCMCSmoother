clup

% Set test flag
for test_flag = [1,3];
    
    test_case = test_flag;
    if test_flag == 3
        test_case = 2;
    end
    
    % Load data
    load(['smoother_test' num2str(test_flag) '.mat']);
    
    fig = figure; hold on
    inds = 2; plot( cellfun(@(x) x.times, results(inds)), cellfun(@(x) x.mean_pos_rmse, results(inds)), 'ok', 'markersize', 7 )
    inds = 3; plot( cellfun(@(x) x.times, results(inds)), cellfun(@(x) x.mean_pos_rmse, results(inds)), 'r^', 'markersize', 7 )
    inds = 4:8; plot( cellfun(@(x) x.times, results(inds)), cellfun(@(x) x.mean_pos_rmse, results(inds)), '*-b', 'markersize', 7 )
    inds = 9:13; plot( cellfun(@(x) x.times, results(inds)), cellfun(@(x) x.mean_pos_rmse, results(inds)), 'x-.g', 'markersize', 7 )
    
    legend('FS', 'DBRS', 'MCMC-BRS (M=1-100)', 'MCMC-BSS (M=1-100)');
    xlabel('Running Time (s)'); ylabel('Position RMSE');
    xlimits = get(gca, 'XLim');
    ylimits = get(gca, 'YLim');
    xlim([-50, xlimits(2)]);
    plot([0 0], ylimits, ':k');
    ylim(ylimits);
    
    wid = 4; hei = 3;
    format_graph_for_pdf;
    print(fig, '-dpdf', ['case' num2str(test_case) '_smoother_comparison_posRMSE_time.pdf']);
    
    %%% Trajectories
    
    % DEFINE RANDOM SEED
    rand_seed = 1;
    
    % Set random seed
    s = RandStream('mt19937ar', 'seed', rand_seed);
    RandStream.setDefaultStream(s);
    
    % Set parameters
    set_parameters;
    
    if test_flag == 1
        params.bng_var = (pi/720)^2;
        params.rng_var = 0.1;
    elseif test_flag == 2
        params.bng_var = (pi/180)^2;
        params.rng_var = 10;
    elseif test_flag == 3
        params.bng_var = (pi/36)^2;
        params.rng_var = 100;
    end
    
    %% Generate some Bearings only tracking data
    [ t, x, y ] = generate_radar_data;
    
    % Plot it
    xmax = max(max(abs(x(1:2,:))))+5;
    fig = figure; hold on; xlim([-xmax, xmax]); ylim([-xmax, xmax]);
    plot(x(1,:), x(2,:), 'b', 'linewidth', 3);
    [x_obs, y_obs] = pol2cart(y(1,:), y(2,:)); plot(x_obs, y_obs, 'r');
    plot(0, 0, 'xk', 'markersize', 10)
    legend('Target Position', 'Observations', 'Location', 'SouthWest')
    xlabel('x coordinate'); ylabel('y coordinate');
    
    xlim([-300, 50]);
    ylim([-50, 300]);
    
    wid = 4; hei = 4;
    format_graph_for_pdf;
    pos = get(gca, 'Position');
    set(gca, 'Position', [pos(1), pos(2), pos(3), pos(3)]);
    print(fig, '-dpdf', ['case' num2str(test_case) '_trajectory.pdf']);
    
end