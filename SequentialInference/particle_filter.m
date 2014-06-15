% Jens Raaby
% ATDM Assignment 3, June 2013

function [ p,estimates ] = particle_filter( observations,...
                                            initial_noise,...
                                            R_noise,...
                                            Q_noise,...
                                            evaluation_probability_function,...
                                            dynamic_model_function,...
                                            N)
%PARTICLE_FILTER Runs particle filtering over a set of observations
%   OBSERVATIONS is a T by M matrix (M dimensions, T time steps)
%   INITIAL_NOISE is a 1 by N matrix specifying the noise for the initial
%   particles
%   R_noise is a 1 by N matrix with noise used for evaluation probability
%   Q_noise is a 1 by N matrix with noise for the dynamical model
%   EVALUATION_PROBABILITY_FUNCTION is a handle for the p(z|x) probability
%   DYNAMIC_MODEL_FUNCTION is a handle for a update function for the particles 
%   N is number of particles


% assume input is rows of observations
z_dims = size(observations,2);
x_dims = size(initial_noise,2);
T = size(observations,1);

% Generate initial particles from first observation with initial noise
% parameter. This could do with refactoring
p = zeros(N,x_dims);
x0 = [observations(1,:), 0 , 0];
for i = 1:N
   xi = dynamic_model_function(x0,initial_noise); 
   p(i,:) = xi;
end

estimates = zeros(T,x_dims);

% loop over observations
for t=1:T
    
    z = observations(t,:);
    
    % MOTION UPDATE
    % 
    p2 = zeros(N,x_dims);
    for n = 1:N
        % move every particle using supplied function
        % and noise is added!
        p2(n,:) = dynamic_model_function(p(n,:),Q_noise);
    end
    
    p = p2;

    
    % MEASUREMENT UPDATE
    w = zeros(N,1);
    for n = 1:N
        w(n) = evaluation_probability_function(z,p(n,:),R_noise);
    end
    % normalised:
    a = w/sum(w);
    
    
    % RESAMPLING
    %   using Thrun's cicular algorithm which is O(N)
    %   see Thrun 2001 and https://www.udacity.com/course/cs373
    
    p3 = zeros(N,x_dims);
    maxw = max(w);
    % random initial particle (between 1 and N):
    index = 1 + round(rand() * (N-1));

    beta = 0;
    for n = 1:N
        % increment beta by some random amount
        % the choice of 2 * maxw is arbitrary
        % (from Sebastian Thrun's Udacity course)
        beta = beta + (rand() * 2 * maxw);
        
        while beta > w(index)
            % decrease beta by current weight
            beta = beta - w(index);
            
            % increment the index by 1
            index = index + 1;
            if (index == N+1)
                index = 1; % modular doesn't work with 0 based indexing
            end
        end
        % we have sampled a particle, add it to the set:
        p3(n,:) = p(index,:);
    end
    
    % store this estimate - could also use another mean (weighted)
    estimate = mean(p3,1);
    estimates(t,:) = estimate;
    
    % replace the initial set of particles
    p = p3;
end

end
