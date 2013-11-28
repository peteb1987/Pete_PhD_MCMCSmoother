function [ time, state, observ ] = tracking_generatedata( model )
%tracking_generatedata Generate a data set for thye tracking model.

% Initialise arrays
time = 1:model.K;
state = zeros(model.ds, model.K);
observ = zeros(model.do, model.K);

% First state
[state(:,1), ~] = tracking_stateprior(model);

% Loop through time
for kk = 1:model.K
    
    if kk > 1
        [state(:,kk), ~] = tracking_transition(model, state(:,kk-1));
    end
    
    [observ(:,kk), ~] = tracking_observation(model, state(:,kk));
    
end

end

