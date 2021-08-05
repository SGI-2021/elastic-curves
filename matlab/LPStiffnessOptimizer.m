classdef LPStiffnessOptimizer < handle
    % Implements the basic computational methods for finding the optimal
    % stiffness distribution of a beam given its equilibrium shape.
    
    properties
        gamma % a 2-by-num_samples matrix of curve samples; each column contains the x- and y-coordinates of a sample
        kappa % a 1-by-num_samples row vector of curvature samples
        gamma_infl = zeros(2,0) % a 2-by-n matrix of inflection points; each column contains the x- and y-coordinates of an inflection point
        num_samples % the number of curve samples (and curvature samples)
    end
    
    methods
        % The arguments have the same format as the corresponding
        % properties. The third argument is optional.
        function obj = LPStiffnessOptimizer(gamma, kappa, gamma_infl)
            obj.gamma = gamma;
            obj.kappa = kappa;
            if nargin > 2
                obj.gamma_infl = gamma_infl;
            end
            obj.num_samples = length(kappa);
        end
        
        % Compute the "best" stiffness profile for a curve without
        % inflections. Returns the stiffness samples as an 1-by-num_samples
        % row vector "K", and the corresponding values for the scalar "a"
        % and the 2-by-1 column vector "b".
        function [K, a, b] = optimizeSimple(obj)
            K = zeros(1, obj.num_samples);
            a = 0;
            b = [0;0];
        end
        
        % Like "optimizeSimple" but adds equality constraints for
        % inflection points to the linear program explicitly.
        function [K, a, b] = optimizeWithInflections(obj)
            K = zeros(1, obj.num_samples);
            a = 0;
            b = [0;0];
        end
    end
end