function [ ps, diagnostics ] = refreshed_backward_simulation_smoother( display, algo, model, fh, pf, observ, true_state )
%BACKWARD_SIMULATION_SMOOTHER Run refreshed backward simulation smoothing
%using MH sampling.

if display.text
    fprintf(1, 'Running refreshed backward simulation smoother.\n');
end

% Local copies
K = model.K;
N = algo.N;
S = algo.S;
M = algo.M;

% Set up diagnostics
diagnostics.acceptance_rate = zeros(S,K);

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
    
    % Loop backwards in time
    for kk = K:-1:1
        
        if kk < K
            next_state = traj_state(:,kk+1);
        else
            next_state = [];
        end
        
        % Initialise
        num_acceptances = 0;
        if kk > 1
            prev_ind = pf(kk).ancestor(traj_ind(kk));
        else
            prev_ind = 0;
        end
        state = pf(kk).state(:,traj_ind(kk));
        
        if kk > 1
            [import_mn, import_vr] = double_kf_update(model, pf(kk-1).state(:,prev_ind), next_state, observ(:,kk));
        else
            [import_mn, import_vr] = double_kf_update(model, [], next_state, observ(:,kk));
        end
        import = loggausspdf(state, import_mn, import_vr);
        
        if kk < K
            [~, trans1] = tracking_transition(model, state, next_state);
        else
            trans1 = 0;
        end
        if kk > 1
            [~, trans2] = tracking_transition(model, pf(kk-1).state(:,prev_ind), state);
        else
            [~, trans2] = tracking_stateprior(model, state);
        end
        [ ~, lhood ] = tracking_observation( model, state, observ(:,kk) );
        
        % Chain
        for jj = 1:M
            
            if kk > 1
                ppsl_prev_ind = sample_weights(pf(kk-1).weight, 1);
                [ppsl_mn, ppsl_vr] = double_kf_update(model, pf(kk-1).state(:,ppsl_prev_ind), next_state, observ(:,kk));
            else
                ppsl_prev_ind = 0;
                [ppsl_mn, ppsl_vr] = double_kf_update(model, [], next_state, observ(:,kk));
            end
            ppsl_state = mvnrnd(ppsl_mn', ppsl_vr)';
            ppsl_import = loggausspdf(ppsl_state, ppsl_mn, ppsl_vr);
            
            if kk < K
                [~, ppsl_trans1] = tracking_transition(model, ppsl_state, next_state);
            else
                ppsl_trans1 = 0;
            end
            if kk > 1
                [~, ppsl_trans2] = tracking_transition(model, pf(kk-1).state(:,ppsl_prev_ind), ppsl_state);
            else
                [~, ppsl_trans2] = tracking_stateprior(model, ppsl_state);
            end
            [ ~, ppsl_lhood ] = tracking_observation( model, ppsl_state, observ(:,kk) );
            
            ap = (ppsl_trans1 + ppsl_trans2 + ppsl_lhood - ppsl_import) ...
                - (trans1 + trans2 + lhood - import);
            if log(rand) < ap
                num_acceptances = num_acceptances + 1;
                prev_ind = ppsl_prev_ind;
                state = ppsl_state;
                trans1 = ppsl_trans1;
                trans2 = ppsl_trans2;
                lhood = ppsl_lhood;
                import = ppsl_import;
            end
            
        end
        
        if kk > 1
            traj_ind(kk-1) = prev_ind;
        end
        traj_state(:,kk) = state;
        diagnostics.acceptance_rate(ii,kk) = num_acceptances / M;
        
    end
    
    ps(ii).traj_ind = traj_ind;
    ps(ii).traj_state = traj_state;
    
end


end

function [mn, vr] = double_kf_update(model, prev_state, next_state, obs)

% KF Predict - from previous state
if ~isempty(prev_state)
    trans_mn = model.A * prev_state;
    trans_vr = model.Q;
else
    trans_mn = model.m1;
    trans_vr = model.P1;
end

% KF Update
if ~isempty(next_state)
    [bi_trans_mn, bi_trans_vr] = kf_update(trans_mn, trans_vr, next_state, model.A, model.Q);
else
    bi_trans_mn = trans_mn;
    bi_trans_vr = trans_vr;
end

% UKF update
h = @(x, par)tracking_h(model, x);
[mn, vr] = ukf_update1(bi_trans_mn, bi_trans_vr, obs, h, model.R);
vr = (vr+vr')/2;

end

