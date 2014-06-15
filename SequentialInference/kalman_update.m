% Jens Raaby
% ATDM Assignment 3, June 2013

function [ xest, P ] = kalman_update( x, A, B, u, P, z, H, R, Q )
%KALMAN_UPDATE Note: state and observation should be column vectors
% x = initial estimate of state
% A = state transition matrix
% B = external motion matrix
% u = external motion
% P = uncertainty covariance (for prediction)
% z = observation
% H = observation model
% R = observation noise
% Q = transition noise

% find the dimensionality 
N = length(x); % state space
M = length(z); % observation space
I = eye(N); % identity matrix in state space

% PREDICTION UPDATE
xhat = A*x + B*u;
Phat = (A * P * A') + Q;% + Q is transition noise


% MEASUREMENT UPDATE
y = z - (H*xhat); % what is error between observation and state

S = (H * Phat * H') + R; % "innovation covariance"
K = (Phat * H') * inv(S); % Kalman Gain

% Estimate
xest = xhat + (K*y);
P = (I - (K*H)) * Phat;


end

