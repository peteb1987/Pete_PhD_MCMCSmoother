function [ pos_rmse, vel_rmse, tnees, nus ] = particle_smoother_analysis( model, true_state, ps )
%PARTICLE_SMOOTHER_ANALYSIS Calculate RMSE, TNEES and NUS over time for
%a particle smoother approximation

K = model.K;
S = length(ps);

% Initialise arrays
nees = zeros(1,K);
nus = zeros(1,K);

% Cut up state
true_pos = true_state(1:model.ds/2,:);
true_vel = true_state(model.ds/2+1:model.ds,:);

% Make a great big array of all particle states
state = cat(3,ps.traj_state);

% Take the mean to give us a point estimate
state_est = mean(state, 3);
pos_est = state_est(1:model.ds/2,:);
vel_est = state_est(model.ds/2+1:model.ds,:);

% RMSEs
pos_rmse = sqrt(sum( (true_pos-pos_est).^2, 1 ));
vel_rmse = sqrt(sum( (true_vel-vel_est).^2, 1 ));

% NUS
for kk = 1:K
    nus(kk) = size(unique(squeeze(state(:,kk,:))', 'rows'),1);
end

% TNEES
for kk = 1:K
    if nus(kk)>model.ds
        state_diff = squeeze(bsxfun(@minus, state(:,kk,:), state_est(:,kk)));
        vr_est = state_diff*state_diff'/S;
        err = state_est(:,kk) - true_state(:,kk);
        nees(kk) = err'*(vr_est\err);
    else
        nees(kk) = Inf;
    end
end
tnees = nees./(1+nees);
tnees(isnan(tnees))=1;



end

