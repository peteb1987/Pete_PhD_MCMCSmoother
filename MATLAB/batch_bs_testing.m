% Batch testing script for BS algorithms

% Add toolboxes to path
addpath('../toolbox/user');
addpath('../toolbox/lightspeed');
addpath('../toolbox/ekfukf');

% Clean up
clup

% Get environment variable specifying test number
sys_num = str2double(getenv('SGE_TASK_ID'));
if isnan(sys_num)
    sys_num = 1;
end

fprintf(1, 'System count number: %u.\n', sys_num);

% Function handles for the model
addpath('tracking3');
fh.setmodel = @tracking_setmodel;
fh.setalgo = @tracking_setalgo;
fh.generatedata = @tracking_generatedata;
fh.transition = @tracking_transition;
fh.observation = @tracking_observation;
fh.stateprior = @tracking_stateprior;
fh.ekfproposal = @tracking_ekfproposal;
fh.ukfproposal = @tracking_ukfproposal;

% Set flags
test.flag_batch = true;

% Set display options
display.text = true;
display.plot_during = false;
display.plot_after = false;

% Calculate test case and number
test_num = ceil(sys_num/3);
test_case = rem(sys_num-1,3)+1;

% Tests to run
test.name = ['model_case_' num2str(test_case)];
test.samp_types = [0 1 2 2 2 2 2 3 3 3 3 3 4 4 4 4];
test.chain_lengths = [NaN NaN 3 10 32 100 316 3 10 32 100 316 3 10 32 100];

% Define random seed
rand_seed = test_num;

% Set model parameters
[model] = feval(fh.setmodel, test);
switch test_case
    case 1
        sigtheta = ( 0.1*(pi/180) )^2;  % Bearing covariance
        sigphi   = ( 0.1*(pi/180) )^2;  % Elevation covariance
        sigr     = 0.01;                % Range covariance
    case 2
        sigtheta = ( 0.1*(pi/180) )^2;  % Bearing covariance
        sigphi   = ( 0.1*(pi/180) )^2;  % Elevation covariance
        sigr     = 100;                 % Range covariance
    case 3
        sigtheta = ( 5*(pi/180) )^2;    % Bearing covariance
        sigphi   = ( 5*(pi/180) )^2;    % Elevation covariance
        sigr     = 0.01;                % Range covariance
    otherwise
        error('Invalid test case');
end
model.R = diag([sigtheta sigphi sigr]);

% Run the script
run_bs_testing;

% Get the numbers we want
results.rt = rt;
results.pos_rmse = pos_rmse;
results.vel_rmse = vel_rmse;
results.tnees = tnees;
results.nus = nus;

% Save
save(['batch_tests/' test.name '/bs_test_' num2str(test_num) '.mat'], 'results', 'model', 'test')

% SHOW WE'VE FINISHED
disp(['Test:' num2str(test_num) ' for case ' num2str(test_case) ' - DONE!']);
