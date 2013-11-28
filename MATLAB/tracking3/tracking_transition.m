function [ new_state, prob ] = tracking_transition( model, state, new_state )
%tracking_transition Sample and/or evaluate observation density for the
%tracking model.

% prob is a log-probability.

% Calculate new_state mean
mn = model.A * state;

% Sample state if not provided
if (nargin<3)||isempty(new_state)
    new_state = mvnrnd(mn', model.Q)';
end

% Calculate probability if required
if nargout>1
    prob = loggausspdf(new_state, mn, model.Q);
else
    prob = [];
end

end

