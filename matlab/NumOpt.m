classdef NumOpt < handle
    methods(Static)
        % Generic Newtons's method for finding a root of the function given
        % by the handle f_fun. df_fun is a function handle for the
        % function's derivative. x is the starting parameter. A root is
        % returned if its absolute function value is below threshold.
        % Terminates after at most max_iter iterations.
        function [ret, converged] = newtonsMethod(f_fun, df_fun, x, threshold, max_iter)
            f = f_fun(x);
            iter = 0;
            while (abs(f) > threshold && iter < max_iter)
                df = df_fun(x);
                x = x - f/df;
                f = f_fun(x);
                iter = iter + 1;
            end
            ret = x;
            converged = abs(f) <= threshold;
        end
    end
end