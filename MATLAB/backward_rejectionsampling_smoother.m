function [ smooth_pts ] = backward_rejectionsampling_smoother( S, times, pts_array, wts_array, h_trans )
%BACKARD_SAMPLING_SMOOTHER Godsill2004-type backward sampling smoother
%using rejection sampling
%
%
%
global params; %UGLY

fprintf(1, '\n\n*** Running direct backward sampling smoother. ***\n');

% Initialise variables
K = size(pts_array,1);

% Initialise arrays
smooth_pts = cell(S,1);

% Sample the final filter weights
parents = systematic_resample(wts_array{K}, S);

% Loop through trajectories
for ii = 1:S
        
    fprintf(1, 'Now simulating trajectory %u.\n', ii);
    
    % Copy filtering particle
    smooth_pts{ii} = pts_array{K}{parents(ii)};
    
    % Loop backwards in time
    for kk = K-1:-1:1
        
        dt = times(kk+1) - times(kk);
        
        % Get weights
        wts = wts_array{kk};
        
        % Rejection sample
        [~, Q] = build_AQ(dt, params.proc_var);  % This is ugly and should be encapsulated.
        max_trans_prb = log(det(2*pi*Q));
        u = inf; thresh = 1;
        count = 0;
        while u > thresh
            
            % Sample an index and a uniform variable
            ind = sample_weights(exp(wts));
            u = log(rand);
            
            % Calculate the acceptance threshold
            [~, trans_prb] = feval(h_trans, dt, pts_array{kk}{ind}(:,kk), smooth_pts{ii}(:,kk+1) );
            thresh = trans_prb - max_trans_prb;
            
            count = count + 1;
            
        end
        
        count
        
        % Update trajectory
        smooth_pts{ii} = [pts_array{kk}{ind}(:,1:kk) smooth_pts{ii}(:,kk+1:end)];
        
    end
    
end

end

