function [ smooth_pts ] = backard_sampling_smoother( S, times, pts_array, wts_array, h_trans )
%BACKARD_SAMPLING_SMOOTHER Normal Godsill2004-type backward sampling
%smoother.
%
%
%

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
        
        % Loop through filter particles
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
%         fprintf(1, '     ESS in frame %u: %f\n', kk, ESS(wts));
        ind = sample_weights(exp(wts));
        
        % Update trajectory
        smooth_pts{ii} = [pts_array{kk}{ind}(:,1:kk) smooth_pts{ii}(:,kk+1:end)];
        
    end
    
end

end

