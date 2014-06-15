function [ h ] = plotPendulum( Z, Z_true, estimates, range, ttl )
%PLOTPENDULUM Plots the real measurements against the estimates
% Z is observations
% Z_true is the ground truth
% ESTIMATES is the estimates
% RANGE is which observations to include
% TTL is the title for the plot

if nargin < 4
   range = 1:length(Z); 
end

h = figure();

plot(Z_true(range,1),Z_true(range,2),'g','LineWidth',3);
hold on

plot(Z(range,1),Z(range,2),'b-');
plot(estimates(range,1),estimates(range,2),'r-');

legend('Ground truth','Observed','Estimates');

if nargin > 4
    title(ttl)
end

