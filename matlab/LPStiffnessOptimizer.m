classdef LPStiffnessOptimizer < handle
    % Implements the basic computational methods for finding the optimal
    % stiffness distribution of a beam given its equilibrium shape.
    
    properties
        gamma % a 2-by-num_samples matrix of curve samples; each column contains the x- and y-coordinates of a sample
        kappa % a 1-by-num_samples row vector of curvature samples
        gamma_infl = zeros(2,0) % a 2-by-n matrix of inflection points; each column contains the x- and y-coordinates of an inflection point
        num_samples % the number of curve samples (and curvature samples)
        err % error message set to 0 if there are no errors and 1 if an error has occured
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

            % Initialise error to 0
            obj.err = 0;
        end
        
        % Compute the "best" stiffness profile for a curve without
        % inflections. Returns the stiffness samples as an 1-by-num_samples
        % row vector "K", and the corresponding values for the scalar "a"
        % and the 2-by-1 column vector "b".
        function [K, a, b] = optimizeSimple(obj)
            f = [0,0,0,1];
            
            % Setting the inequality matrix
            A1 = [-ones(obj.num_samples, 1) ./ obj.kappa', ...
                  -obj.gamma(1,:)' ./ obj.kappa', ...
                  -obj.gamma(2,:)' ./ obj.kappa', ... 
                  zeros(obj.num_samples, 1)];
            A2 = [ones(obj.num_samples, 1) ./ obj.kappa', ...
                  obj.gamma(1,:)' ./ obj.kappa', ...
                  obj.gamma(2,:)' ./ obj.kappa', ...
                  -ones(obj.num_samples, 1)];
            A = [A1; A2];

            b = [-ones(obj.num_samples, 1), zeros(obj.num_samples, 1)];

            abM = linprog(f,A,b); % solving the linear programme

            K = -(abM(1:3)' * A(1:obj.num_samples,1:3)');
            a = abM(1);
            b = abM(2:3);
            
        end
        
        % Like "optimizeSimple" but adds equality constraints for
        % inflection points to the linear program explicitly.
        function [K, a, b] = optimizeWithInflections(obj)
            f = [0,0,0,1];

            % Setting the inequality matrix
            A1 = [-ones(obj.num_samples, 1) ./ obj.kappa', ...
                  -obj.gamma(1,:)' ./ obj.kappa', ... 
                  -obj.gamma(2,:)' ./ obj.kappa', ... 
                  zeros(obj.num_samples, 1)];
            A2 = [ones(obj.num_samples, 1) ./ obj.kappa', ...
                  obj.gamma(1,:)' ./ obj.kappa', ...
                  obj.gamma(2,:)' ./ obj.kappa', ... 
                  -ones(obj.num_samples, 1)];
            A = [A1; A2];

            b = [-ones(obj.num_samples, 1), zeros(obj.num_samples, 1)];
            
            % Setting the equality matrix
            Aeq = [ones(size(obj.gamma_infl,2),1), obj.gamma_infl(1,:)', obj.gamma_infl(2,:)', zeros(size(obj.gamma_infl,2),1)];
            beq = zeros(size(obj.gamma_infl,2),1);
            
            % Checking whether there are any infinities in the matrix A and
            % returning an error if so
            if(sum(isinf(A)) == 0) 
                abM = linprog(f,A,b,Aeq,beq); % solving the linear programme
                obj.err = 0;
            else
                obj.err = 1;
                fprintf('Curve is infeasible \n');
                abM = [];
            end
            
            % Checking whether or not a solution exists and returning an
            % error if not
            if(isempty(abM))
                obj.err = 1;
                fprintf('Curve is infeasible \n');
                a = [];
                b = [];
                K = [];
            else
                K = -(abM(1:3)' * A(1:obj.num_samples,1:3)');
                a = abM(1);
                b = abM(2:3);
                obj.err = 0;
            end
        end
    end
end