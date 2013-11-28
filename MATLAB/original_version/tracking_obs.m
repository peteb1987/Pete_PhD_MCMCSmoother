function [ y_mn, obs_prb ] = tracking_obs( x, y )
%TRACKING_OBS Observation likelihood for bearing-only 2D radar observations

global params;

R = [params.bng_var 0; 0 params.rng_var];

[b, r] = cart2pol(x(1), x(2));
y_mn = [b; r];

if ~isempty(y)
    obs_prb = fast_log_mvnpdf(y', y_mn', R);
end

end

