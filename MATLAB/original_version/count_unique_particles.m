function [ Nup, Nuh ] = count_unique_particles( pts )
%COUNT_UNIQUE_PARTICLES Counts number of unique particles (at time k) and
%histories (up to time k)

% Initialise variables
Np = length(pts);
[d, K] = size(pts{1});

% Initialise arrays
Nup = Np*ones(1,K);
Nuh = Np*ones(1,K);

fprintf(1, '\nCounting unique particles.\n');

% Loop through time
for k = 1:K
    
    if rem(k, 100)==0
        fprintf(1, '     Reached time step %u.\n', k);
    end
    
    %%%% Particles
    
    unique = true(Np,1);
    
    % Loop through particles
    for ii = 1:Np
        
        x_ii = pts{ii}(:,k);
        
        % Loop through later particles
        for jj = ii+1:Np
            
            % Only look at particles which we still know to be unique
            if unique(jj)
                
                x_jj = pts{jj}(:,k);
                
                % Compare number of jumps
                if all(all(x_ii==x_jj))
                    
                    % Particle is not unique, cross it off
                    unique(jj) = false;
                    Nup(k) = Nup(k) - 1;
                    
                end
                
            end
            
        end
        
    end
    
    %%%% Histories
    
    unique = true(Np,1);
    
    % Loop through particles
    for ii = 1:Np
        
        x_ii = pts{ii}(:,1:k);
        
        % Loop through later particles
        for jj = ii+1:Np
            
            % Only look at particles which we still know to be unique
            if unique(jj)
                
                x_jj = pts{jj}(:,1:k);
                
                % Compare number of jumps
                if all(all(x_ii==x_jj))
                    
                    % Particle is not unique, cross it off
                    unique(jj) = false;
                    Nuh(k) = Nuh(k) - 1;
                    
                end
                
            end
            
        end
        
    end
    
end

fprintf(1, '     Average number of unique particles: %3.1f.\n', mean(Nup));
fprintf(1, '     Number of unique trajectories:      %u.\n', Nuh(end));

end

