% Jens Raaby
% ATDM Assignment 3, June 2013
function [error] = evaluation_probability(z, x, noise)
%EVALUATION_PROBABILITY evaluation probability for simple state spaces
% We assume the first N dimensions of x are the same as those in z
% noise by default is 0.5 in each dimension
% loopy is a flag to use a loop instead of matlab's vectorised code
    z_dims = length(z);
    
    % default noise
    if nargin < 3
        noise = 0.5 * ones(1,z_dims);  
    end
        
    % for each dimension, compute probability of observation given some
    % noise
    gs = normpdf(z,x(1:z_dims),noise);
    
    % could also use mvnpdf - probably better metric
%     noise2 = eye(z_dims);
%     noise2(1,1) = noise(1);
%     noise2(2,2) = noise(2);
%     g2 = mvnpdf(z,x(1:z_dims),noise2);
    

    % the overall probability is the product of these probabilities
    error = prod([1,gs]);
% error = g2;
end