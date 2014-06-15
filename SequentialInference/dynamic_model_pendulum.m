function [ x_now ] = dynamic_model_pendulum( x_prev, x_noise, dt)
%DYNAMIC_MODEL_PENDULUM A simple dynamic model for 2D motion
%   We just use the velocity and timestep to estimate next position
%   The velocity is the same as before.
%   Gaussian noise is added for each dimension

x_dims = length(x_prev);


% default to zero noise
if nargin<2
   x_noise = zeros(1,x_dims);
end

% default timestep of 1
if nargin<3
    dt = 1;
end

% compute new state without noise
% first copy old state:
x_now = x_prev;

% update positions:
x_now(1) = x_prev(1) + (dt * x_prev(3));
x_now(2) = x_prev(2) + (dt * x_prev(4));

% now add noise:
for i = 1:x_dims
    
    % generate gaussian random number (mean 0, variance specified by noise)
    noise = x_noise(i).*randn();
    
    % add the noise to the state
    x_now(i) = x_now(i) + noise;
end

end

