function [ pts_array, wts_array, filter_pts ] = particle_filter( init_pts, times, observs, h_ppsl, h_trans, h_obs, ess_thresh )
%PARTICLE_FILTER Run a particle filter
%
% *** Input ***
% init_pts      - Cell array of particles defining prior (see below)
% times         - Row vector of observation times
% observs       - Matrix of observations, observs(:,k) is the kth observation
% h_ppsl        - Function handle to the proposal function.
% h_trans       - Function handle to the transition likelihood function.
% h_obs         - Function handle to the observation likelihood function.
%
% *** Output ***
% pts_array     - Cell array of particle cell arrays, the filter output at each observation time
% wts_array     - Cell array of corresponding weight column vectors
% filter_pts    - Cell array of equal-weight, current filter points (i.e. no Kitigawa smoothing)
%
% Particle cell arrays have the following structure: pts{ii}(:,k) is the
% state of the iith particle at the kth time.
%
% All weights and probabilities are log-ed.
%

fprintf(1, '\n\n*** Running particle filter. ***\n');

% Initialise constants
Np = length(init_pts);
K = size(times, 2);
d = size(init_pts{1}, 1);

% Initialise arrays
pts_array = cell(K,1);
wts_array = cell(K,1);
filter_pts = repmat({zeros(d,K)}, Np, 1);
wts = log(ones(Np,1)/Np);
last_pts = init_pts;

% Loop through time
for kk = 1:K
    
    fprintf(1, 'Now processing frame %u.\n', kk);
    
    % Create a new particle array
    pts = repmat({zeros(d,kk)}, Np, 1);
    
    % Find time increment
    t = times(kk);
    if kk == 1
        dt = t;
    else
        dt = times(kk) - times(kk-1);
    end
    
    % Loop through particles
    for ii = 1:Np
        
        % Copy state trajectory from old particle
        pts{ii} = last_pts{ii}(:,1:kk-1);
        
        % Get previous state
        if kk == 1
            last_x = init_pts{ii}(:,1);
        else
            last_x = pts{ii}(:,kk-1);
        end
        
        % Propose a new value for the particle
        [x, ppsl_prb] = feval(h_ppsl, last_x, observs(:,kk), dt, h_trans, h_obs);
        
        % Store new value
        pts{ii}(:,kk) = x;
        
        % Evaluate probabilities
        [~, trans_prb] = feval(h_trans, dt, last_x, x);
        [~, obs_prb] = feval(h_obs, x, observs(:,kk));
        
        % Update weight
        wts(ii) = wts(ii) + trans_prb + obs_prb - ppsl_prb;
        
    end
    
    % Normalise weights
    wts = wts - max(wts);
    assert(~any(isnan(wts)), 'Invalid set of weights');
    lin_wts = exp(wts); lin_wts = lin_wts/sum(lin_wts);
    wts = log(lin_wts);
    
    % Store particles and weights
    pts_array{kk} = pts;
    wts_array{kk} = wts;
    
    % Resample
    parents = systematic_resample(wts);
    resamp_pts = pts(parents);
    if ESS(wts) < ess_thresh*Np
        pts = resamp_pts;
        wts = log(ones(Np,1)/Np);
        fprintf(1, '     Resampled in frame %u.\n', kk);
    end
    for ii = 1:Np
        filter_pts{ii}(:,kk) = resamp_pts{ii}(:,kk);
    end
    
    last_pts = pts;
    
end

end

