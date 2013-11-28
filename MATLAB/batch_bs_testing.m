% Batch testing script for BS algorithms

% Add toolboxes to path
addpath('../toolbox/user');
addpath('../toolbox/lightspeed');
addpath('../toolbox/ekfukf');

% Clean up
clup

% Get environment variable specifying test number
test_num = str2double(getenv('SGE_TASK_ID'));
if isnan(test_num)
    test_num = 0;
end
fprintf(1, 'Test number: %u.\n', test_num);

% Set flags
test.flag_batch = true;

% Define random seed
rand_seed = test_num;

% Set display options
display.text = true;
display.plot_during = false;
display.plot_after = false;

% Tests to run
test.name = 'model_case_2'
test.samp_types = [0 1 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4];
test.chain_lengths = [NaN NaN 3 10 32 100 320 3 10 32 100 320 3 10 32 100 320];
  
% Run the script
run_bs_testing;

% Get the numbers we want
results.rt =   rt;
results.pos_rmse = pos_rmse;
results.vel_rmse = vel_rmse;
results.tnees = tnees;
results.nus = nus;

% Save
save(['batch_tests/' test.name '/bs_test_' num2str(test_num) '.mat'], 'results', 'model', 'test')

% SHOW WE'VE FINISHED
disp(['Test:' num2str(test_num) ' - DONE!']);
