function [ x_mn, trans_prb ] = tracking_trans( dt, last_x, x )
%TRACKING_TRANS Transition probability for linear 2D tracking dynamics

global params;

[A, Q] = build_AQ(dt, params.proc_var);

x_mn = A*last_x;

if ~isempty(x)
    trans_prb = fast_log_mvnpdf(x', x_mn', Q);
end

end

