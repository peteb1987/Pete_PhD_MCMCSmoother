function [ smooth_pts ] = mcmc_newstate_smoother( S, M, times, pts_array, wts_array, observs, h_trans, h_obs, h_bidirec_ppsl )
%MCMC_NEWSTATE_SMOOTHER 
%
%
%

fprintf(1, '\n\n*** Running MCMC-based backward sampling smoother. ***\n');

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
    for kk = K-1:-1:2
        
        dt_next = times(kk+1) - times(kk);
        dt_prev = times(kk) - times(kk-1);
        
        old_x = smooth_pts{ii}(:,kk);
        
        % Initialise probability
        [~, old_trans_prb_prev] = feval(h_trans, dt_prev, smooth_pts{ii}(:,kk-1), old_x );
        [~, old_trans_prb_next] = feval(h_trans, dt_next, old_x, smooth_pts{ii}(:,kk+1) );
        [~, old_obs_prb] = feval(h_obs, smooth_pts{ii}(:,kk), observs(:,kk));
        [~, old_ppsl_prb] = feval(h_bidirec_ppsl, smooth_pts{ii}(:,kk-1), smooth_pts{ii}(:,kk+1), observs(:,kk), dt_prev, dt_next, h_trans, h_obs, smooth_pts{ii}(:,kk));
        
%         ind = 1;
%         while smooth_pts{ii}(1,kk-1)~=pts_array{kk-1}{ind}(1,kk-1)
%             ind = ind + 1;
%         end
%         old_weight = wts_array{kk-1}(ind);
        
        for jj = 1:M
            
            % Propose a new k-1 particle
%             ind = randsample(Np, 1);
            ind = randsample(Np, 1, true, exp(wts_array{kk-1}));
            
            % Propose a new k state
            [new_x, new_ppsl_prb] = feval(h_bidirec_ppsl, pts_array{kk-1}{ind}(:,kk-1), smooth_pts{ii}(:,kk+1), observs(:,kk), dt_prev, dt_next, h_trans, h_obs);
            
            % Calculate probability
            [~, new_trans_prb_prev] = feval(h_trans, dt_prev, pts_array{kk-1}{ind}(:,kk-1), new_x );
            [~, new_trans_prb_next] = feval(h_trans, dt_next, new_x, smooth_pts{ii}(:,kk+1) );
            [~, new_obs_prb] = feval(h_obs, new_x, observs(:,kk));
%             new_weight = wts_array{kk}(ind);
            
            % Construct acceptance probability
            ap = (old_ppsl_prb - new_ppsl_prb)...
                +(new_trans_prb_prev - old_trans_prb_prev)...
                +(new_trans_prb_next - old_trans_prb_next)...
                +(new_obs_prb - old_obs_prb);%...
%                 +(new_weight - old_weight);
            
            if log(rand)<ap
                % Update trajectory
                smooth_pts{ii} = [pts_array{kk-1}{ind}(:,1:kk-1), new_x, smooth_pts{ii}(:,kk+1:end)];
%                 old_weight = new_weight;
                old_trans_prb_prev = new_trans_prb_prev;
                old_trans_prb_next = new_trans_prb_next;
                old_obs_prb = new_obs_prb;
%                 old_x = new_x;
            end
            
        end
        
    end
    
end

end

