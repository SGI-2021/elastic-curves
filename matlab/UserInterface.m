classdef UserInterface < handle
    % User interface which allows the user to edit a spline curve in real
    % time and get the strip needed to form that curve
    
    properties(Access = private)
        cp  % control points
    end

    properties 
        fig % figure opened
        figstrip % second figure
        ax % axes
        axstrip % axes for strip
        
        lines %lines connecting control points
        splineplot % drawing of spline
        strip % optimised strip
        strip_patch % patch graphics object for drawing the strip outline

        spline % spline curve object
        num_samples % number of samples taken

        pointpatch % patch graphics object for drawing the control points
        mouseflag % set to true if the mouse is down and near a point, set to false otherwise
        selectedInd % current index
    end
    
    methods
         function obj = UserInterface(cp)
            % Set up the figure for plotting the spline such that points
            % can be moved around with the mouse
            obj.fig = figure('WindowButtonDownFcn', @(source, event) obj.mouseDown(source),...
                             'WindowButtonUpFcn', @(source, event) obj.mouseUp(source, event), ...
                             'WindowButtonMotionFcn', @(source, event) obj.mouseMove(source, event));

            obj.ax = axes(obj.fig);
            title('Spline curve');
            xlim([-1,1]);
            ylim([-1,1]);
            obj.cp = cp;
            
            % Initialise mouse flag to false
            obj.mouseflag = false;
            
            hold on;
            % Plot of the control points
            obj.pointpatch = patch('Vertices', cp', ...
                                   'Faces', [1:size(cp,2)]', ...
                                   'Parent', obj.ax, ...
                                   'Marker', '.');
            % Plot of lines between adjacent control points
            obj.lines = plot(obj.cp(1,:), obj.cp(2,:),'k--');

            % Plot of the spline curve
            obj.spline = SplineCurve(3, obj.cp);
            
            % Initialise selectedInd and numSamples
            obj.selectedInd = 0;
            obj.num_samples = 400;

            % Sample the curve, and approximate the arc-length parameter s at every
            % sample
            t_samples = linspace(0,obj.spline.t_max, obj.num_samples);
            gamma = obj.spline.evaluate(t_samples);
            obj.splineplot = plot(gamma(1,:), gamma(2,:),'LineWidth',2,'Color',[1 0 0]);

            % Plotting a second figure with the elastic strip
            obj.figstrip = figure();
            obj.axstrip = axes(obj.figstrip);
            obj.strip_patch = patch(obj.axstrip);
            axis tight equal;
            title('Elastic strip');

            obj.drawStrip();

        end

        function mouseDown(obj, source, event)
            % Called when a mouse button is pushed 
            % LEFT CLICK DRAG: MOVE POINT
            % SHIFT CLICK: ADD POINT
            % CONTROL CLICK : DELETE POINT

            if(strcmp(source.SelectionType, 'normal')) 
                % Checks to see if a left-click is made near a point and updates
                % the mouseflag to be used by mouseMove              
                
                [d, obj.selectedInd] = obj.closestPoint();

                if(d < 0.1) 
                    obj.mouseflag = true;
                else
                    obj.mouseflag = false;
                end

            % ADDING POINTS (SHIFT CLICK):
            elseif(strcmp(source.SelectionType, 'extend')) 
                % Checks if a shift-click is made near a control point
                % Adds another control point in between  the control point
                % that was shift-clicked and the next one. If the last 
                % control point is shift-clicked, adds a control point 
                % after the the previous end point

                if(isempty(obj.cp))
                    % If there are no control points add a control point
                    % where the cursor is
                    currpoint = obj.ax.CurrentPoint;
                    coord = [currpoint(1,1); currpoint(1,2)];
                    
                    % Draw control point
                    obj.updateSpline(coord);
                    obj.drawPoints();
                    
                    return
                end

                if(size(obj.cp,2) == 1)
                    % If there is only one control point, add another 
                    % control point where the cursor is
                    currpoint = obj.ax.CurrentPoint;
                    coord = [currpoint(1,1); currpoint(1,2)];
                    c = [obj.cp, coord];

                    % Draw control point
                    obj.updateSpline(c);
                    obj.drawPoints();

                end

                % Computes nearest control point
                [d, closestInd] = obj.closestPoint();

                if(d < 0.1)
                    % Adds control point
                    if(closestInd ~= size(obj.cp,2)) 
                        % Find midpoint between selected point and next
                        % point
                        midpoint = ((obj.cp(:,closestInd) + obj.cp(:,closestInd+1))'./2)';
                        
                        % Add the midpoint to the list of control points at
                        % the index after the selected point
                        c1 = obj.cp(:,1:closestInd);
                        c2 = midpoint;
                        c3 = obj.cp(:,closestInd+1:end);
                        c = [c1 c2 c3];

                        obj.updateSpline(c);

                    elseif(closestInd == size(obj.cp,2))
                        % Find midpoint between last point and penultimate
                        % point
                        midpoint = ((obj.cp(:,end) + obj.cp(:,end-1))'./2)'; 
                        direction = midpoint - obj.cp(:,end-1);

                        % Create a new last point that is the same distance
                        % away from the last point as the midpoint is from
                        % the last point
                        newlastpoint = obj.cp(:,end) + direction;
                        
                        %Add the new last point to the control points
                        c = [obj.cp, newlastpoint];
                        obj.updateSpline(c);
                    end

                    % Draws contol point
                    obj.drawPoints();
                    obj.drawSpline();
                    obj.drawStrip();

                else
                    obj.mouseflag = false;
                end 
                               
            % DELETING POINTS (CONTROL-CLICK)
            elseif(strcmp(source.SelectionType, 'alt'))
                % Checks for a control-click
                % Deletes the control point that is control-clicked
                
                [d, closestInd] = obj.closestPoint();

                if(d < 0.1) 
                    c = obj.cp;
                    c(:,closestInd) = [];

                    obj.updateSpline(c);

                    obj.drawPoints();
                    obj.drawSpline();
                    obj.drawStrip();

                else
                    obj.mouseflag = false;
                end
                
            end
            
        end

        % Updates mouseflag when mouse is not being pressed
        function mouseUp(obj, source, event)
            % Set mouseflag to false since mouse is up
            obj.mouseflag = false;
            
            % Reset selectedInd
            obj.selectedInd = 0;
        end

        % Moves a control point 
        function mouseMove(obj, source, event)
            
            if(obj.mouseflag == true)
                 % Activates when the mouse is pressed down and moving
                 
                 currpoint = obj.ax.CurrentPoint;
                 coords = [currpoint(1,1); currpoint(1,2)];
                 c = obj.cp;
                 c(:, obj.selectedInd) = coords;
                 obj.updateSpline(c);

                 % Draws points, spline and strip in real time
                 obj.drawPoints();
                 obj.drawSpline();
                 obj.drawStrip();
                 
            end
        
        end
        
        % Finds the closest control point when mouse is clicked
        function [M, Ind] = closestPoint(obj)
            
            % 2x3 matrix where mouse click is
            currpoint = obj.ax.CurrentPoint;
            % x and y coordinates of the mouse click
            coord = [currpoint(1,1); currpoint(1,2)];

            % Distance between mouse click and each of the points
            d = obj.cp - coord;
            dnorm = sum(d.^2,1);

            % Returns the index of the point closest to the mouse as well
            % as the distance of this point from the mouse
            [M, Ind] = min(dnorm);

        end
        
        % Updates and draws control points
        function drawPoints(obj)
            obj.pointpatch.Vertices = obj.cp';
            obj.pointpatch.Faces = [1:size(obj.cp,2)]';
        end

        % Draws spline from the control poins
        function drawSpline(obj)
            obj.lines.XData =  obj.cp(1,:);
            obj.lines.YData = obj.cp(2,:);
            
            t_samples = linspace(0,obj.spline.t_max, obj.num_samples);
            gamma = obj.spline.evaluate(t_samples);
            obj.splineplot.XData = gamma(1,:);
            obj.splineplot.YData = gamma(2,:);
        end

        % Updates spline information
        function updateSpline(obj, cp)
            obj.spline.setControlPoints(cp);
            obj.cp = cp;
        end
        
        % Draws the corresponding elastic strip
        function drawStrip(obj)
            % Sample the curve, and approximate the arc-length parameter s at every
            % sample
            t_samples = linspace(0,obj.spline.t_max,obj.num_samples);
            gamma = obj.spline.evaluate(t_samples);
            kappa = obj.spline.curvature(t_samples);
            to = gamma(:,2:end)-gamma(:,1:end-1);
            seg_lens = sqrt(sum(to.^2,1));
            s = [0 cumsum(seg_lens)];
            
            % Initialize stiffness optimizer
            [~, gamma_infl] = obj.spline.findInflectionPoints();
            lpopt = LPStiffnessOptimizer(gamma, kappa, gamma_infl);
            
            % Solve the linear programme to find K
            K = lpopt.optimizeWithInflections();

            if (lpopt.err == 0)
                % Generate simple strip outline and plot
                max_width = 0.05;
                half_width = 0.5*max_width/max(K)*K;
                p_simple = GeometryGenerator.generateProfileOutlineLoop(s,half_width,0.01);
                SvgTools.exportCurves('spline-curve/strip_simple.svg', {p_simple}, 1e3/0.352778);
    
                obj.strip_patch.Vertices = p_simple';
                obj.strip_patch.Faces = [1:size(p_simple,2)-1; 2:size(p_simple,2)]';
                
                figure(obj.figstrip); % make the strip figure current, so the next line affects that figure, and not the spline figure
                title('Elastic strip');
                axis tight equal;

            elseif(lpopt.err == 1)
                % Do not plot anything if the curve is infeasible
                obj.strip_patch.Vertices = [];
                obj.strip_patch.Faces = [];
                
                figure(obj.figstrip);
                title('\color{red}Curve infeasible');
                return;  
            end

        end

    end

end

