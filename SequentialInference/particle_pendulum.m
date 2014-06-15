clear all
clc



% Pendulum tracking
% read in the data
noisy_pend = 'noisy_pendulum.csv';
true_pend =  'true_pendulum.csv';

f = noisy_pend;
Z = csvread(['data/' noisy_pend]);
Zperfect = csvread(['data/' true_pend]);
T = size(Z,1); % how many observations

% we need to create parameters
% the dynamic model is needed for estimation

% observation space:
z_dims = size(Z,2);
R = [0.5,0.5]; % noise in the measurements

% State space
% just try 2D position and 2D velocity:
% [x,y,vx,vy]
x_dims = 4;
x_noise_initial = [0.9,0.9,0.9,0.9]; 
Q = [0.5,0.5,0.5,0.5]; % dynamic model noise

% time step
dt = 1;

% Create a function to evaluate an observation conditioned on state
% p(z | x)
evalprob = @(z,x,R) evaluation_probability(z,x,R);

% Create a function to update a particle with the dynamical model
transitionprob = @(x0,noise) dynamic_model_pendulum(x0,noise,dt);

% how many particles:
N = 1000;

[particles,estimates] = particle_filter(Z, x_noise_initial, R, Q, evalprob, transitionprob, N);


h1 = plotPendulum(Z,Zperfect,estimates,301:400,'Particle Filter, final 100 observations');
h2 = plotPendulum(Z,Zperfect,estimates,1:400,'Particle Filter, all observations');
h3 = plotOneDim(Z,Zperfect,estimates,1,301:400,'X axis, Particle filter, final 100 observations');
h4 = plotOneDim(Z,Zperfect,estimates,2,301:400,'Y axis, Particle filter, final 100 observations');


% root mean squared error function for both dimensions:
rms1 = @(compare1,est,rg) sqrt(mean(sum((compare1(rg,1:2) - est(rg,1:2)).^2,2)));
% for each dimension
rms = @(compare1,est,rg) sqrt(mean((compare1(rg,1:2) - est(rg,1:2)).^2,1));

rg = 301:400;

rms_truth = rms(Zperfect,estimates,rg)
rms_truth2 = rms1(Zperfect,estimates,rg)

rms_observed = rms(Z,estimates,rg)
rms_observed2 = rms1(Z,estimates,rg)
%%
print(h1,'-depsc2','Report/figures/Particle_last100_trace.eps');
print(h2,'-depsc2','Report/figures/Particle_all_trace.eps');
print(h3,'-depsc2','Report/figures/Particle_last100_X.eps');
print(h4,'-depsc2','Report/figures/Particle_last100_Y.eps');


