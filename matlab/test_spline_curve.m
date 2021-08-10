% Import the spline curve and scale it to be about 30 cm long
spline = SplineCurve.import('rect_spline1.txt');
scale = 0.08;
spline.cp = spline.cp * scale;
num_samples = 400;

% Sample the curve, and approximate the arc-length parameter s at every
% sample
t_samples = linspace(0,spline.t_max,num_samples);
gamma = spline.evaluate(t_samples);
kappa = spline.curvature(t_samples);
to = gamma(:,2:end)-gamma(:,1:end-1);
seg_lens = sqrt(sum(to.^2,1));
s = [0 cumsum(seg_lens)];

% Plot the B-Spline curve
figure;
hold on;
scatter(spline.cp(1,:), spline.cp(2,:),32,'k','s','filled','MarkerEdgeColor','none');
plot(spline.cp(1,:), spline.cp(2,:),'k--');
plot(gamma(1,:), gamma(2,:),'LineWidth',2,'Color',[1 0 0]);
title('Spline Curve');
axis tight equal;

% Initialize stiffness optimizer
[~, gamma_infl] = spline.findInflectionPoints();
lpopt = LPStiffnessOptimizer(gamma, kappa, gamma_infl);

% TODO1: This is the first method that needs to be implemented
%K = lpopt.optimizeSimple();

% TODO2: Once the 'simple' version works, you can also implement this
% method, and compare the results on a curve with an inflection point.
K = lpopt.optimizeWithInflections();

if isempty(K)
    fprintf('Curve is infeasible!');
    return;
end


% Generate simple strip outline and plot
max_width = 0.05;
half_width = 0.5*max_width/max(K)*K;
p_simple = GeometryGenerator.generateProfileOutlineLoop(s,half_width,0.01);
SvgTools.exportCurves('spline-curve/strip_simple.svg', {p_simple}, 1e3/0.352778);

figure;
title('Simple Strip Outline');
patch('Vertices',p_simple','Faces',[1:size(p_simple,2)-1; 2:size(p_simple,2)]');
axis tight equal;


% Generate a perforated strip outline and plot
p_perforated = GeometryGenerator.generateMicroPerforationProfile(s,half_width*2,...
    max_width*ones(1,num_samples),3,5,0.2,5e-3);
SvgTools.exportForCricut('spline-curve/strip_perforated.svg', {p_perforated}, 1e3/0.352778);

figure;
title('Perforated Strip Outline');
for li=1:length(p_perforated)
    patch('Vertices',p_perforated{li}','Faces',[1:size(p_perforated{li},2)-1;2:size(p_perforated{li},2)]');
end
axis tight equal;


% Generate outlines for a composite strip and plot
hw2 = 0.5*max_width*ones(1,num_samples);
% Use materials, where the broad strip is thickness 0.1, and the narrow
% strip is thickness 0.2 (only the ratio between the materials matters!)
hw1 = GeometryGenerator.computeSubstripHalfWidth(0.2,0.1,hw2,K);
if isempty(hw1)
    fprintf('Composite Strip is infeasible.');
    return;
end

comp1 = GeometryGenerator.generateProfileOutlineLoop(s,hw1,0.01);
comp2 = GeometryGenerator.generateProfileOutlineLoop(s,hw2,0.01);

SvgTools.exportCurves('spline-curve/strip_comp_thick.svg', {comp1}, 1e3/0.352778);
SvgTools.exportCurves('spline-curve/strip_comp_thin.svg', {comp2}, 1e3/0.352778);

figure;
title('Composite Strip Outlines');
patch('Vertices',comp1','Faces',[1:size(comp1,2)-1; 2:size(comp1,2)]');
patch('Vertices',(comp2 + [0;0.1])','Faces',[1:size(comp2,2)-1; 2:size(comp2,2)]');
axis tight equal;