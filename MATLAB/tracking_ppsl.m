function [ x, ppsl_prb ] = tracking_ppsl( last_x, y, dt, h_trans, h_obs)
%TRACKING_PPSL Propose a new state for a 2D tracking model

global params;

[A, Q] = build_AQ(dt, params.proc_var);
R = [params.bng_var 0; 0 params.rng_var];

% Optimal importance
[m, P] = kf_predict(last_x, eps*eye(4), A, Q);
H = [-m(2)/(m(1)^2+m(2)^2),    m(1)/(m(1)^2+m(2)^2),     0, 0;
     m(1)/sqrt(m(1)^2+m(2)^2), m(2)/sqrt(m(1)^2+m(2)^2), 0, 0];
[m, P] = ekf_update1(m, P, y, H, R, h_obs);
x = mvnrnd(m', P)';
ppsl_prb = log(mvnpdf(x', m', P));

% % Bootstrap
% x = mvnrnd((A*last_x), Q)';
% ppsl_prb = log(mvnpdf(x', (A*last_x)', Q));

end

