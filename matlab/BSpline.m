classdef BSpline
    % A class that generates basis spline functions
    % Used by 'SplineCurve.m'
    
    properties
        degree % degree of the spline (2 for quadratic, 3 for cubic etc.)
        pieces % array of 'Polymonial' instances, of length degree+1
    end
    
    methods
        % Represents a BSpline basis function of the given degree
        function obj = BSpline(degree)
            obj.degree = degree;
            cur(degree+1) = Polynomial(0.0);
            cur(1) = Polynomial(1.0);
            for k=1:degree
                prev = cur;
                kr = 1 / k;
                cur(1) = Polynomial([0;kr]) * prev(1);
                for i=1:k-1
                    cur(i+1) = Polynomial([i*kr;kr]) * prev(i+1) + Polynomial([1-(i-1)*kr; -kr]) * prev(i);
                end
                cur(k+1) = Polynomial([kr;-kr]) * prev(k);
            end
            obj.pieces = cur;
        end
    end
end

