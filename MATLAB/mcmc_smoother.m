function [ smooth_pts ] = mcmc_smoother( S, M, times, pts_array, wts_array, h_trans )
%MCMC_SMOOTHER New MCMC-based backward sampling smoother
%
%
%

fprintf(1, '\n\n*** Running MCMC-based backward resampling smoother. ***\n');

% Initialise variables
K = size(pts_array,1);
Np = size(pts_array{K},1);

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
        
        % Initialise probability
        [~, old_trans_prb] = feval(h_trans, dt, smooth_pts{ii}(:,kk), smooth_pts{ii}(:,kk+1) );
%         ind = 1;
%         while smooth_pts{ii}(1,kk)~=pts_array{kk}{ind}(1,kk)
%             ind = ind + 1;
%         end
%         old_weight = wts_array{kk}(ind);
        
        for jj = 1:M
            
            % Propose a new particle
%             ind = randsample(Np, 1);
            ind = sample_weights(exp(wts_array{kk}));
            
            % Calculate probability
            [~, new_trans_prb] = feval(h_trans, dt, pts_array{kk}{ind}(:,kk), smooth_pts{ii}(:,kk+1) );
            new_weight = wts_array{kk}(ind);
            
            % Construct acceptance probability
            ap = (new_trans_prb - old_trans_prb);%...
%                 +(new_weight - old_weight);
            
            if log(rand)<ap
                % Update trajectory
                smooth_pts{ii} = [pts_array{kk}{ind}(:,1:kk) smooth_pts{ii}(:,kk+1:end)];
                old_weight = new_weight;
                old_trans_prb = new_trans_prb;
            end
            
        end
        
    end
    
end

end

