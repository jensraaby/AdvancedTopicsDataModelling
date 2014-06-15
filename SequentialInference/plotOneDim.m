function [ h ] = plotOneDim( Z, Z_true, estimates, dim, range, ttl )
%PLOTPENDULUM Plots the real measurements against the estimates
% Z is observations
% Z_true is the ground truth
% ESTIMATES is the estimates
% DIM is which dimension to plot
% RANGE is which observations to include
% TTL is the title for the plot


if nargin < 5
   range = 1:length(Z); 
end
if nargin < 4
    dim = 1;
end

h = figure();

plot(Z_true(range,dim),'g','LineWidth',1);
hold on

plot(Z(range,dim),'bo-');
plot(estimates(range,dim),'rx-');

legend('Ground truth','Observed','Estimates');

if nargin > 5
    title(ttl);
end

