% Jens Raaby
% ATDM Assignment 3, June 2013

clear all
clc

% Load in the data
noisy_pend = 'noisy_pendulum.csv';
true_pend =  'true_pendulum.csv';

Z        = csvread(['data/' noisy_pend]);
Zperfect = csvread(['data/' true_pend]);

measurements = size(Z,1);

% 4D state space
% State is: x,y,vx,vy
% dt is the time step used for the dynamic model:
dt = 1;

% State transition matrix
A = [1 0 dt 0;...
     0 1 0 dt;...
     0 0 1 0 ;...
     0 0 0 1];


% Measurement function
% how to transform from state to measurement
H = [1 0 0 0;...
     0 1 0 0];

P = 100 * eye(size(A,1)); % initial state covariance - makes little difference

% R is the measurement uncertainty:
R = eye(2);
R(1,1) = 0.5; % lower to reduce effect of Kalman gain
R(2,2) = 1; % lower to reduce effect of Kalman gain

% Q - transition noise:
Q = 5 * eye(size(A));

% external motion terms:
u = [0;0;0;0]; 
B = zeros(4);

x = Z(1,:)'; % initial observation is used for initial state
x = [x; 0; 0]; % assume everything else is static


% LOOP THROUGH THE OBSERVATIONS
estimates = zeros(measurements,size(A,1));
for i = 1:measurements
   
    % load observation:
   z = Z(i,:)';
   
   % compute estimate using kalman_update
   [x1, P] = kalman_update( x, A, B, u, P, z, H, R, Q );
   
   % cache the output
   estimates(i,:) = x1'; 
end

%%

h1 = plotPendulum(Z,Zperfect,estimates,301:400,'Kalman Filter, final 100 observations');
h2 = plotPendulum(Z,Zperfect,estimates,1:400,'Kalman Filter, all observations');
h3 = plotOneDim(Z,Zperfect,estimates,1,301:400,'X axis, Kalman filter, final 100 observations');
h4 = plotOneDim(Z,Zperfect,estimates,2,301:400,'Y axis, Kalman filter, final 100 observations');

print(h1,'-depsc2','Report/figures/Kalman_last100_trace.eps');
print(h2,'-depsc2','Report/figures/Kalman_all_trace.eps');
print(h3,'-depsc2','Report/figures/Kalman_last100_X.eps');
print(h4,'-depsc2','Report/figures/Kalman_last100_Y.eps');

% root mean squared error function for both dimensions:
rms1 = @(compare1,est,rg) sqrt(mean(sum((compare1(rg,1:2) - est(rg,1:2)).^2,2)));
% for each dimension
rms = @(compare1,est,rg) sqrt(mean((compare1(rg,1:2) - est(rg,1:2)).^2,1));

rg = 301:400;

rms_truth = rms(Zperfect,estimates,rg)
rms_truth2 = rms1(Zperfect,estimates,rg)

rms_observed = rms(Z,estimates,rg)
rms_observed2 = rms1(Z,estimates,rg)

% fnameprefix = sprintf('4D_R%d'

%%
% 6D state space:% position, velocity and acceleration
% state is x,y,vx,vy,ax,ay

% dt is the time step used for the dynamic model:
dt = 1;

% State transition matrix
A = [1 0 dt 0 dt^2 0;...
     0 1 0 dt 0    dt^2;...
     0 0 1 0  dt   0;...
     0 0 0 1  0    dt;...
     0 0 0 0  1    0;...
     0 0 0 0  0    1];


% Measurement function
% how to transform from state to measurement
H = [1 0 0 0 0 0;...
     0 1 0 0 0 0];

P = 100 * eye(size(A,1)); % initial uncertainty - makes little difference

% R is the measurement uncertainty:
R = eye(2);
R(1,1) = 0.5;
R(2,2) = 0.8;

Q = 5 * eye(6); % transition uncertainty
u = [0;0;0;0;0;0]; % external motion
% 
B = zeros();

x = Z(1,:)'; % initial observation is used for initial state
x = [x; 0; 0; 0; 0]; % assume everything else is static

estimates = zeros(measurements,6);
for i = 1:measurements
   
    % load observation:
   z = Z(i,:)';
   
   % compute estimate using kalman_update
   [x1, P] = kalman_update(x,A,B,u,P,z,H,R,Q);
   
   % cache the output
   estimates(i,:) = x1'; 
end

h1 = plotPendulum(Z,Zperfect,estimates,301:400,'Kalman Filter, final 100 observations');
h2 = plotPendulum(Z,Zperfect,estimates,1:400,'Kalman Filter, all observations');
h3 = plotOneDim(Z,Zperfect,estimates,1,301:400,'X axis, Kalman filter, final 100 observations');
h4 = plotOneDim(Z,Zperfect,estimates,2,301:400,'Y axis, Kalman filter, final 100 observations');

