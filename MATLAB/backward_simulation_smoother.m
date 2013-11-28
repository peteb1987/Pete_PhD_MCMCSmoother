function [ ps, diagnostics ] = backward_simulation_smoother( display, algo, model, fh, pf, observ, flag_samp_type, true_state )
%BACKWARD_SIMULATION_SMOOTHER Run backward simulation smoothing with one of
%a variety of sampling methods

% flag_samp_type
% 0 = filter
% 1 = direct (exhaustive evaluation)
% 2 = rejection sampling with early stopping
% 3 = mh sampling

if display.text
    fprintf(1, 'Running backward simulation smoother, type %u.\n', flag_samp_type);
end

% Local copies
K = model.K;
N = algo.N;
S = algo.S;
M = algo.M;

% If we're just resampling the filter particles, keep N of them
if flag_samp_type == 0
    S = N;
end

% Set up diagnostics
switch flag_samp_type  
    case 2
        diagnostics.time_out = false(S,K);
    case 3
        diagnostics.acceptance_rate = zeros(S,K);
    otherwise
        diagnostics = [];
end

% Resample final frame filter particles
filt_ancestor = sample_weights(pf(K).weight, S);

% Particle loop
for ii = 1:S
    
    if display.text
        fprintf(1, '   Simulating trajectory %u.\n', ii);
    end
    
    % Initialise new particle index and state trajectory
    traj_ind = zeros(1,K);
    traj_state = zeros(model.ds,K);
    traj_ind(K) = filt_ancestor(ii);
    traj_state(:,K) = pf(K).state(:,traj_ind(K));
    
    % Loop backwards in time
    for kk = K-1:-1:1
        
        % Get next state
        next_state = traj_state(:,kk+1);
        
        % Sample by appropriate method
        switch flag_samp_type
            
            case 0
                % Reproduce the filter particle
                next_ind = traj_ind(kk+1);
                traj_ind(kk) = pf(kk+1).ancestor(next_ind);
                traj_state(:,kk) = pf(kk).state(:,traj_ind(kk));
            
            case 1
                % Direct sampling
                bs_weight = zeros(N,1);
                for jj = 1:N
                    [~, trans] = tracking_transition(model, pf(kk).state(:,jj), next_state);
                    bs_weight(jj) = pf(kk).weight(jj) + trans;
                end
                traj_ind(kk) = sample_weights(bs_weight, 1);
                traj_state(:,kk) = pf(kk).state(:,traj_ind(kk));
                
            case 2
                % Rejection sampling
                max_trans = -0.5*log(det(2*pi*model.Q));
                done = false;
                for jj = 1:M
                    ind = sample_weights(pf(kk).weight, 1);
                    [~, trans] = tracking_transition(model, pf(kk).state(:,ind), next_state);
                    assert(trans < max_trans);
                    if log(rand) < (trans-max_trans)
                        traj_ind(kk) = ind;
                        traj_state(:,kk) = pf(kk).state(:,traj_ind(kk));
                        done = true;
                        break
                    end
                end
                if ~done
                    % Timed out - use direct sampling
                    diagnostics.time_out(ii,kk) = true;
                    bs_weight = zeros(N,1);
                    for jj = 1:N
                        [~, trans] = tracking_transition(model, pf(kk).state(:,jj), next_state);
                        bs_weight(jj) = pf(kk).weight(jj) + trans;
                    end
                    traj_ind(kk) = sample_weights(bs_weight, 1);
                    traj_state(:,kk) = pf(kk).state(:,traj_ind(kk));
                end
                
            case 3
                % MH sampling
                num_acceptances = 0;
                next_ind = traj_ind(kk+1);
                ind = pf(kk+1).ancestor(next_ind);
                state = pf(kk).state(:,ind);
                [~, trans] = tracking_transition(model, state, next_state);
                for jj = 1:M
                    ppsl_ind = sample_weights(pf(kk).weight, 1);
                    ppsl_state = pf(kk).state(:,ppsl_ind);
                    [~, ppsl_trans] = tracking_transition(model, ppsl_state, next_state);
                    if (ind ~= ppsl_ind) && (log(rand) < (ppsl_trans-trans))
                        num_acceptances = num_acceptances + 1;
                        ind = ppsl_ind;
                        state = ppsl_state;
                        trans = ppsl_trans;
                    end
                end
                traj_ind(kk) = ind;
                traj_state(:,kk) = state;
                diagnostics.acceptance_rate(ii,kk) = num_acceptances / M;
                
            otherwise
                error('Not a valid sampling method')
                
        end
        
    end
    
    ps(ii).traj_ind = traj_ind;
    ps(ii).traj_state = traj_state;
    
end


end

