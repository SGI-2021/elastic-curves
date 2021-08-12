function [a,b] = spline_intersection(spline1, spline2)
    % Generate a initial guess for the parameters a,b such that gamma1(a)
    % is close to gamma2(b). Let's call this initial guess a0,b0.
    
    % Start Newton's iteration with x0=(a0,b0), and run iterations until we
    % find an and bn, such that ||gamma1(an) - gamma2(bn)|| < epsilon.
    % (epsilon = 10^-12).
end