function [ pf, diagnostics ] = particle_filter( display, algo, model, fh, observ, flag_ppsl_type, true_state )
%particle_filter

if display.text
    fprintf(1, 'Running particle filter, type %u.\n', flag_ppsl_type);
end

% Initialise diagnostics structure
diagnostics = repmat(struct('rt', [], 'ess', [], 'se', []), 1, model.K);

% Initialise particle filter structure
pf = repmat(struct('state', [], 'ancestor', [], 'weight', []), 1, model.K);

%%% Time loop %%%
for kk = 1:model.K
    
    if display.text
        fprintf(1, '   Time step %u.\n', kk);
    end
    
    % Start timer
    tic;
    
    %%% Particle selection %%% - CUSTOMISE THIS IF NEEDED
    if kk > 1
        aux_weight = pf(kk-1).weight;
        pf(kk).ancestor = sample_weights(aux_weight, algo.N, 2);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Initialise state and weight arrays
    pf(kk).state = zeros(model.ds, algo.N);
    pf(kk).weight = zeros(1, algo.N);
    step_count = zeros(1, algo.N);
    
    switch flag_ppsl_type
        case {1, 2, 3, 4}
            %%% Ordinary particle filters %%%
            
            % Particle loop
            for ii = 1:algo.N
                
                % Ancestory
                if kk > 1
                    prev_state = pf(kk-1).state(:,pf(kk).ancestor(ii));
                    prior_weight = pf(kk-1).weight(pf(kk).ancestor(ii)) - aux_weight(pf(kk).ancestor(ii));
                else
                    prior_weight = 0;
                end
                
                % Sample from importance density
                switch flag_ppsl_type
                    case 1
                        %%% Bootstrap %%%
                        
                        if kk > 1
                            [state, ~] = feval(fh.transition, model, prev_state);
                        else
                            [state, ~] = feval(fh.stateprior, model);
                        end
                        
                        % Density
                        [~, lhood_prob] = feval(fh.observation, model, state, observ(:,kk));
                        
                        % Weight
                        weight = prior_weight + lhood_prob;
                        
                    case {2, 3, 4}
                        %%% Other importance densities %%%
                        
                        switch flag_ppsl_type
                            
                            case 2
                                %%% EKF approx to OID %%%
                                if kk > 1
                                    [state, ppsl_prob] = feval(fh.ekfproposal, model, prev_state, observ(:,kk));
                                else
                                    [state, ppsl_prob] = feval(fh.ekfproposal, model, [], observ(:,kk));
                                end
                                
                            case 3
                                %%% UKF approx to OID %%%
                                if kk > 1
                                    [state, ppsl_prob] = feval(fh.ukfproposal, model, prev_state, observ(:,kk));
                                else
                                    [state, ppsl_prob] = feval(fh.ukfproposal, model, [], observ(:,kk));
                                end
                                
                            case 4
                                %%% Linearised OID %%%
                                if kk > 1
                                    [state, ppsl_prob] = feval(fh.linearisedoidproposal, model, prev_state, observ(:,kk));
                                else
                                    [state, ppsl_prob] = feval(fh.linearisedoidproposal, model, [], observ(:,kk));
                                end
                                
                        end
                        
                        % Densities
                        if kk > 1
                            [~, trans_prob] = feval(fh.transition, model, prev_state, state);
                        else
                            [~, trans_prob] = feval(fh.stateprior, model, state);
                        end
                        [~, lhood_prob] = feval(fh.observation, model, state, observ(:,kk));
                        
                        % Weight
                        weight = prior_weight + lhood_prob + trans_prob - ppsl_prob;
                        
                        if isinf(ppsl_prob)
                            weight = -inf;
                        end
                        
                end
                
                assert(isreal(weight));
                assert(~isnan(weight));
                
                % Store new state and weight
                pf(kk).state(:,ii) = state;
                pf(kk).weight(ii) = weight;
                
            end
            
        otherwise
            error('Invalid choice of proposal type');
            
    end
    
    % Make sure we have real weights
    assert(all(isreal(pf(kk).weight)));
    
    % Handle the case where all the particles have weight 0
    if all(isinf(pf(kk).weight))
        pf(kk).weight = zeros(size(pf(kk).weight));
    end
    
    % Particle mean and variance
    norm_weight = exp(normalise_weights(pf(kk).weight));
    pf(kk).mn = pf(kk).state*norm_weight';
    err = bsxfun(@minus, pf(kk).state, pf(kk).mn);
    pf(kk).vr = err*diag(norm_weight)*err';
    
    % Diagnostics
    diagnostics(kk).rt = toc;
    diagnostics(kk).ess = calc_ESS(pf(kk).weight);
    diagnostics(kk).se = true_state(:,kk) - pf(kk).mn;
    diagnostics(kk).step_count = step_count;
    
    if display.text
        fprintf(1, '      - Took %fs.\n', diagnostics(kk).rt);
        fprintf(1, '      - ESS of %f.\n', diagnostics(kk).ess);
    end
    
end

end

