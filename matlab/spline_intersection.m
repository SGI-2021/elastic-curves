function [a, b] = spline_intersection(spline1, spline2)
    % Generate a initial guess for the parameters a,b such that gamma1(a)
    % is close to gamma2(b). Let's call this initial guess a0,b0.
    
    % Start Newton's iteration with x0=(a0,b0), and run iterations until we
    % find an and bn, such that ||gamma1(an) - gamma2(bn)|| < epsilon.
    % (epsilon = 10^-12).
    
    % T1, T2 are the paramter spaces of both splines
    T1 = linspace(0, spline1.t_max);
    T2 = linspace(0, spline2.t_max);
    
    % Discetized approximations of gamma (100 sample points)
    gamma1 = spline1.evaluate(T1);
    gamma2 = spline2.evaluate(T2);
    
    % Initialize distance cell
    distances = cell(1, size(gamma1, 2));
    
    % Computes distance between i-th sample point on gamma1 and all other
    % 100 sample points of gamma 2
    for i = 1:size(gamma1, 2)
        distances{i} = sum((gamma2 - gamma1(:,i)).^2, 1);   
    end
    
    % Computes linear Index of minimum distance 
    [~, minInd] = min([distances{:}]);
    
    % Translates linear indicies into linear indicies for the splines
    % i.e the minimum distance is between sp1point and sp2point
    sp1point = ceil(minInd / 100);
    sp2point = mod(minInd, 100);
    if sp2point == 0
        sp2point = 100;
    end
    
    % Computes the parameters a and b
    a = T1(sp1point);
    b = T2(sp2point);
        
    % TODO: Implement Newton's Method
end