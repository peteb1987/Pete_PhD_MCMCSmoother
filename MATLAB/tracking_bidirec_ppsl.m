function [x, ppsl_prb] = tracking_bidirec_ppsl(x_prev, x_next, y, dt_prev, dt_next, h_trans, h_obs, x)
%TRACKING_BIDIREC_PPSL Propose a new state for a 2D tracking model when the
% next state is already known (for smoother).

global params;

[A_prev, Q_prev] = build_AQ(dt_prev, params.proc_var);
[A_next, Q_next] = build_AQ(dt_next, params.proc_var);
R = [params.bng_var 0; 0 params.rng_var];

% Interpolate between previous and next x
[m, P] = kf_update(A_prev*x_prev, Q_prev, x_next, A_next, Q_next);

% Now include y
H = [-m(2)/(m(1)^2+m(2)^2),    m(1)/(m(1)^2+m(2)^2),     0, 0;
     m(1)/sqrt(m(1)^2+m(2)^2), m(2)/sqrt(m(1)^2+m(2)^2), 0, 0];
[m, P] = ekf_update1(m, P, y, H, R, h_obs);

% Sample
if (nargin < 8)||isempty(x)
    x = mvnrnd(m', P)';
end
ppsl_prb = log(mvnpdf(x', m', P));

end

