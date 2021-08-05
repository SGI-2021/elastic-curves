classdef SplineCurve < handle
    properties
        spline % BSpline instance, determines the basis functions for interpolation
        degree % spline degree
        cp % dim-by-ncp matrix of control points
        dim % dimension of control points
        ncp % num control points
        msi % highest segment index (starting count from 0)
        t_max % Parameter space [0,t_max]
    end
    
    properties(Constant)
        tess_default = 20
    end
    
    methods(Static)
        % Reads a spline from a file
        % Format:
        % First line contains two integers separated by a blankspace: the
        % degree of the curve, and the number of control points
        % All remaining lines are two floating-point numbers separated by a
        % blankspace, containing the coordinates of a control point
        function obj = import(filename)
            file = fopen(filename);
            deg_n = textscan(file, '%u', 2);
            deg = double(deg_n{1}(1));
            n = double(deg_n{1}(2));
            cp = textscan(file, '%f', n*2);
            cp = reshape(cp{1},2,n);
            
            obj = SplineCurve(deg, cp);
        end
    end
    
    methods
        % Generate a new plane spline curve with degree deg, and control
        % points cp, a 2-by-n matrix.
        function obj = SplineCurve(deg, cp)
            obj.spline = BSpline(deg);
            obj.cp = cp;
            obj.dim = size(obj.cp, 1);
            obj.ncp = size(obj.cp, 2);
            obj.degree = deg;
            obj.t_max = size(cp, 2) - deg;
            obj.msi = obj.t_max - 1;
        end
        
        function ret = evaluateAbstract(obj, t, func)
            t = t(:)';
            segi = min(obj.msi, max(0, floor(t)));
            s = t - segi;
            ret = zeros(obj.dim, length(t));
            for i=1:obj.degree+1
                ret = ret + obj.cp(:,segi+i) .* func(obj.spline.pieces(obj.degree+2-i), s);
            end
        end
        
        function weights = evaluateWeightsAbstract(obj, t, func)
            t = t(:)';
            segi = min(obj.msi, max(0, floor(t)));
            s = t - segi;
            weights = zeros(obj.degree+1, length(t));
            for i=1:obj.degree+1
                weights(i,:) = func(obj.spline.pieces(obj.degree+2-i), s);
            end
        end
        
        function first_cpi = firstControlPointIndexAt(obj, t)
            t = t(:)';
            first_cpi = min(obj.msi, max(0, floor(t))) + 1;
        end
        
        function weights = evaluateWeights(obj, t)
            weights = obj.evaluateWeightsAbstract(t, @(piece,s) piece.evaluate(s));
        end
        
        function weights = evaluateWeightsD(obj, t)
            weights = obj.evaluateWeightsAbstract(t, @(piece,s) piece.evaluateD(s));
        end
        
        function weights = evaluateWeightsD2(obj, t)
            weights = obj.evaluateWeightsAbstract(t, @(piece,s) piece.evaluateD2(s));
        end
        
        % Evaluate the spline curve at parameter t, which can be a row
        % vector
        function ret = evaluate(obj, t)
            ret = obj.evaluateAbstract(t, @(piece,s) piece.evaluate(s));
        end
        
        % Evaluate the first derivative of the spline curve at parameter t,
        % which can be a row vector
        function ret = evaluateD(obj, t)
            ret = obj.evaluateAbstract(t, @(piece,s) piece.evaluateD(s));
        end
        
        % Evaluate the second derivative of the spline curve at parameter t,
        % which can be a row vector
        function ret = evaluateD2(obj, t)
            ret = obj.evaluateAbstract(t, @(piece,s) piece.evaluateD2(s));
        end
        
        % Evaluate the third derivative of the spline curve at parameter t,
        % which can be a row vector
        function ret = evaluateD3(obj, t)
            ret = obj.evaluateAbstract(t, @(piece,s) piece.evaluateD3(s));
        end
        
        % Find all inflection points on the spline curve. t_infl is a
        % vector of parameters at which inflections occur. p_infl is a 
        % 2-by-n matrix containing the inflection points.
        function [t_infl, p_infl] = findInflectionPoints(obj)
            t_samp = linspace(0, obj.t_max, obj.t_max * 10 + 1);
            infl = obj.inflectionIndicator(t_samp);
            curv_pos = infl > 0;
            infl_segs = find(curv_pos(1:end-1) ~= curv_pos(2:end));
            t_infl = zeros(length(infl_segs), 1);
            for i=1:length(infl_segs)
                si = infl_segs(i);
                t_seg = infl(si) / (infl(si) - infl(si+1));
                t_init = t_samp(si) + t_seg * (t_samp(si+1) - t_samp(si));
                t_final = NumOpt.newtonsMethod(@(t) obj.inflectionIndicator(t), @(t) obj.inflectionIndicatorD(t), t_init, 1.e-12, 10);
                t_infl(i) = t_final;
            end
            if nargout > 1
                p_infl = obj.evaluate(t_infl);
            end
        end
        
        function ret = inflectionIndicator(obj, t)
            dr = obj.evaluateD(t);
            ddr = obj.evaluateD2(t);
            ret = dr(1,:) .* ddr(2,:) - dr(2,:) .* ddr(1,:);
        end
        
        function ret = inflectionIndicatorD(obj, t)
            dr = obj.evaluateD(t);
            dddr = obj.evaluateD3(t);
            ret = dr(1,:) .* dddr(2,:) - dr(2,:) .* dddr(1,:);
        end
        
        % Evaluate the curvature at t, a row vector
        function ret = curvature(obj, t)
            dr = obj.evaluateD(t);
            ddr = obj.evaluateD2(t);
            dr_norm_squared = sum(dr.*dr, 1);
            ret = (dr(1,:) .* ddr(2,:) - dr(2,:) .* ddr(1,:)) ./ (dr_norm_squared .^ 1.5);
        end
        
        % Evaluate the curvature derivative at t, a row vector
        function ret = curvatureD(obj, t)
            dr = obj.evaluateD(t);
            ddr = obj.evaluateD2(t);
            dddr = obj.evaluateD3(t);
            dr_norm_squared = sum(dr.*dr, 1);
            ret = (dr(1,:) .* dddr(2,:) - dr(2,:) .* dddr(1,:)) ./ (dr_norm_squared .^ 1.5) ...
                - 3*(dr(1,:) .* ddr(2,:) - dr(2,:) .* ddr(1,:)).*sum(dr.*ddr,1) ./ (dr_norm_squared .^ 2.5);
        end
        
        function ret = curvatureSD(obj, t)
            dr = obj.evaluateD(t);
            ddr = obj.evaluateD2(t);
            dw = obj.evaluateWeightsD(t);
            ddw = obj.evaluateWeightsD2(t);
            dr_norm = sqrt(sum(dr.*dr, 1));
            
            k_num = dr(1,:) .* ddr(2,:) - dr(2,:) .* ddr(1,:);
            k_denom = dr_norm .^ 3;
            
            dkx1 = (dw .* ddr(2,:) - dr(2,:) .* ddw) ./ k_denom;
            dkx2 = (-3. * k_num .* dr(1,:) .* dw) ./ (dr_norm.^5);
            dky1 = (dr(1,:) .* ddw - dw .* ddr(1,:)) ./ k_denom;
            dky2 = (-3. * k_num .* dr(2,:) .* dw) ./ (dr_norm.^5);
            
            ret = zeros(2 * (obj.degree + 1), length(t));
            ret(1:2:end) = dkx1 + dkx2;
            ret(2:2:end) = dky1 + dky2;
        end
        
        function [ds_dq, dp_dq] = inflectionPointSD(obj, t_infl)
            dr = obj.evaluateD(t_infl);
            ddr = obj.evaluateD2(t_infl);
            dddr = obj.evaluateD3(t_infl);
            dk_ds = dr(1,:) * dddr(2,:) - dr(2,:) * dddr(1,:);
            dw = obj.evaluateWeightsD(t_infl)';
            ddw = obj.evaluateWeightsD2(t_infl)';
            dk_dq = [dw * ddr(2); dw * (-ddr(1))] + [-dr(2) * ddw; dr(1) * ddw];
            ds_dq = -dk_dq / dk_ds;
            if nargout > 1
                w = obj.evaluateWeights(t_infl)';
                dp_dq = dr * ds_dq(:)';
                dp_dq(1,1:2:end) = dp_dq(1,1:2:end) + w;
                dp_dq(2,2:2:end) = dp_dq(2,2:2:end) + w;
            end
        end
        
        % Approximate the arc length of the spline curve
        function arc_length = arcLength(obj, tess)
            if nargin < 2
                tess = SplineCurve.tess_default;
            end
            
            t = linspace(0, obj.t_max, tess*obj.t_max+1);
            p = obj.evaluate(t);
            arc_length = sum(sqrt(sum((p(:,2:end) - p(:,1:end-1)).^2, 1)));
        end
        
        function dal = arcLengthSD(obj, tess)
            if nargin < 2
                tess = SplineCurve.tess_default;
            end
            
            t = linspace(0, obj.t_max, tess*obj.t_max+1);
            p = obj.evaluate(t);
            cpi = obj.firstControlPointIndexAt(t);
            cpi_all = cpi + (0:obj.degree)';
            weights = obj.evaluateWeights(t);
            
            to = p(:,2:end) - p(:,1:end-1);
            len = sqrt(sum(to.^2, 1));
            dir = to ./ len;
            
            dal = zeros(2, obj.ncp);
            for i=1:2
                contr_1 = dir(i,:).*weights(:,2:end);
                contr_2 = (-dir(i,:)).*weights(:,1:end-1);
                dal(i,:) = accumarray(...
                    reshape([cpi_all(:,2:end), cpi_all(:,1:end-1)],[],1), ...
                    reshape([contr_1, contr_2],[],1), ...
                    [obj.ncp 1])';
            end
        end
        
        function plotCurve(obj, color, tess, seg_points, line_style)
            if nargin < 2 || isempty(color)
                color = [0 0 0];
            end
            if nargin < 3 || isempty(tess)
                tess = SplineCurve.tess_default;
            end
            if nargin < 4 || isempty(seg_points)
                seg_points = true;
            end
            if nargin < 5 || isempty(line_style)
                line_style = '-';
            end
            
            t = linspace(0,obj.t_max, tess * obj.t_max + 1);
            p = obj.evaluate(t);
            plot(p(1,:),p(2,:),line_style,'Color',color);
            if seg_points
                scatter(p(1,1:tess:end), p(2,1:tess:end),[],color,'.');
            end
        end
        
        function plotInflectionPoints(obj, color)
            [~,p] = obj.findInflectionPoints();
            scatter(p(1,:),p(2,:),16,color,'filled');
        end
        
        % Export the spline curve to a text file, in the same format as the
        % 'SplineCurve.import(...)' function expects
        function writeToFile(obj, path)
            file_out = fopen(path,"w");
            fprintf(file_out, "%u %u\n", obj.degree, obj.ncp);
            fprintf(file_out, "%.15e %.15e\n", obj.cp);
            fclose(file_out);
        end
    end
end