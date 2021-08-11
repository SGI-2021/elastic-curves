close all; 
%clear all;
% cp = [0.2,0.5,0.9; ...
%       0.8,0.5,0.3];

% Import the spline curve and scale it to be about 30 cm long
% spline = SplineCurve.import('rect_spline1.txt');
% scale = 0.08;
% spline.cp = spline.cp * scale;
% num_samples = 400;
% 
% % Sample the curve, and approximate the arc-length parameter s at every
% % sample
% t_samples = linspace(0,spline.t_max,num_samples);
% gamma = spline.evaluate(t_samples);
% kappa = spline.curvature(t_samples);
% to = gamma(:,2:end)-gamma(:,1:end-1);
% seg_lens = sqrt(sum(to.^2,1));
% s = [0 cumsum(seg_lens)];

% Plot the B-Spline curve
% figure;
% hold on;
% scatter(spline.cp(1,:), spline.cp(2,:),32,'k','s','filled','MarkerEdgeColor','none');
% plot(spline.cp(1,:), spline.cp(2,:),'k--');

UI = UserInterface(spline.cp);
%plot(gamma(1,:), gamma(2,:),'LineWidth',2,'Color',[1 0 0]);
title('Spline Curve');
axis tight equal;





