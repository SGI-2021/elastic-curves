close all; 
%clear all;
% cp = [0.2,0.5,0.9; ...
%       0.8,0.5,0.3];

spline = SplineCurve.import('rect_spline1.txt');
scale = 0.08;
spline.cp = spline.cp * scale;

UI = UserInterface(spline.cp);




