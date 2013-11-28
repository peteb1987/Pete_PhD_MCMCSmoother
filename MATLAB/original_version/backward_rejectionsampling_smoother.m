function [ smooth_pts ] = backward_rejectionsampling_smoother( S, times, pts_array, wts_array, h_trans )
%BACKARD_SAMPLING_SMOOTHER Godsill2004-type backward sampling smoother
%using rejection sampling
%
%
%
global params; %UGLY

fprintf(1, '\n\n*** Running rejection sampling smoother. ***\n');

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
        max_trans_prb = -0.5*log(det(2*pi*Q));
        u = inf; thresh = 1;
        count = 0;
        while (u > thresh)
            
            % Sample an index and a uniform variable
            ind = sample_weights(exp(wts));
            u = log(rand);
            
            % Calculate the acceptance threshold
            [~, trans_prb] = feval(h_trans, dt, pts_array{kk}{ind}(:,kk), smooth_pts{ii}(:,kk+1) );
            thresh = trans_prb - max_trans_prb;
            
            assert(thresh < 0)
            
            count = count + 1;
            
            if count == params.M
                % Do it normally
%                 fprintf(1,'No sample - direct evaluation.');
                
                for jj = 1:length(pts_array{kk})
                    % Calculate weight
                    [~, trans_prb] = feval(h_trans, dt, pts_array{kk}{jj}(:,kk), smooth_pts{ii}(:,kk+1) );
                    wts(jj) = wts(jj) + trans_prb;
                end
                
                % Normalise weights
                wts = wts - max(wts);
                lin_wts = exp(wts); lin_wts = lin_wts/sum(lin_wts);
                wts = log(lin_wts);
                
                % Sample
                ind = sample_weights(exp(wts));
                
                break
            end
            
        end
        
%         count
        
        % Update trajectory
        smooth_pts{ii} = [pts_array{kk}{ind}(:,1:kk) smooth_pts{ii}(:,kk+1:end)];
        
    end
    
end

end

