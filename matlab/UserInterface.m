classdef UserInterface < handle
    %USERINTERFACE Summary of this class goes here
    %   Detailed explanation goes here
    % Method for each mouse event (click, drag, up)
    % How to distiguish between adding new points and dragging around
    % doc figure -> see also -> figure properties for documentations
    
    properties%(Access = private)
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

        pointpatch % 
        mouseflag
        selectedInd % current index
    end
    
    methods
         function obj = UserInterface(cp)
            obj.fig = figure('WindowButtonDownFcn', @(source, event) obj.mouseDown(source),...
                             'WindowButtonUpFcn', @(source, event) obj.mouseUp(source, event), ...
                             'WindowButtonMotionFcn', @(source, event) obj.mouseMove(source, event));
            

            obj.ax = axes(obj.fig);
            title('Spline curve');
            xlim([-1,1]);
            ylim([-1,1]);
            obj.cp = cp;

            obj.mouseflag = false;

            hold on;
            obj.pointpatch = patch('Vertices', cp', ...
                                   'Faces', [1:size(cp,2)]', ...
                                   'Parent', obj.ax, ...
                                   'Marker', '.');
            obj.lines = plot(obj.cp(1,:), obj.cp(2,:),'k--');

            obj.spline = SplineCurve(3, obj.cp);

            obj.selectedInd = 0;
            obj.num_samples = 400;
            % Sample the curve, and approximate the arc-length parameter s at every
            % sample
            t_samples = linspace(0,obj.spline.t_max, obj.num_samples);
            gamma = obj.spline.evaluate(t_samples);
            obj.splineplot = plot(gamma(1,:), gamma(2,:),'LineWidth',2,'Color',[1 0 0]);

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
            
            % MOVE POINTS (LEFT CLICK):
            

            

            if(strcmp(source.SelectionType, 'normal')) 
                % Checks to see if a left-click is made and updates
                % the mouseflag to be used by mouseMove              
                obj.mouseflag = true;
                obj.selectedInd = obj.closestPoint();

                currpoint = obj.ax.CurrentPoint;
                coord = [currpoint(1,1); currpoint(1,2)];

 
            % ADDING POINTS (SHIFT CLICK):
            elseif(strcmp(source.SelectionType, 'extend')) 
                % Checks if a shift-click is made near a control point
                % Adds another control point in between  the control point
                % that was shift-clicked and the next one. If the last 
                % control point is shift-clicked, adds a control point in
                % between the last and second from last control points             

                if(isempty(obj.cp))
                    % If there are no control points add a control point
                    % where the cursors is
                    currpoint = obj.ax.CurrentPoint;
                    coord = [currpoint(1,1); currpoint(1,2)];

                    obj.updateSpline(coord);
                    obj.drawPoints();
                    
                    return
                end

                % Computes nearest control point
                closestInd = obj.closestPoint();
                
                % Adds control point
                if(closestInd ~= size(obj.cp,2)) 
                    midpoint = ((obj.cp(:,closestInd) + obj.cp(:,closestInd+1))'./2)';
                    c1 = obj.cp(:,1:closestInd);
                    c2 = midpoint;
                    c3 = obj.cp(:,closestInd+1:end);
                    c = [c1 c2 c3];
                    obj.updateSpline(c);
                else
                    midpoint = ((obj.cp(:,end) + obj.cp(:,end-1))'./2)'; 
                    direction = midpoint - obj.cp(:,end-1);
                    newlastpoint = obj.cp(:,end) + direction;
                    obj.updateSpline([obj.cp, newlastpoint]);
                end
                
                % Draws contol point
                obj.drawPoints();
                obj.drawSpline();
                               
            % DELETING POINTS (CONTROL-CLICK)
            elseif(strcmp(source.SelectionType, 'alt'))
                % Checks for a control-click
                % Deletes the control point that is control-clicked
                
                closestInd = obj.closestPoint();
                c = obj.cp;
                c(:,closestInd) = [];

                obj.updateSpline(c);

                obj.drawPoints();
                obj.drawSpline();
                
            end
        end

     
        % Updates mouseflag when mouse is not being pressed
        function mouseUp(obj, source, event)
            obj.mouseflag = false;
            obj.drawStrip();

            obj.selectedInd = 0;
        end

        % Moves a Control Point 
        function mouseMove(obj, source, event)
            
            if(obj.mouseflag == true)
                 % Activates when the mouse is pressed down and moving
                 
                 currpoint = obj.ax.CurrentPoint;
                 coords = [currpoint(1,1); currpoint(1,2)];
                 c = obj.cp;
                 c(:, obj.selectedInd) = coords;
                 obj.updateSpline(c);

                 
                 obj.drawPoints();
                 obj.drawSpline();
                 
            end
        
        end
        
        % Finds the closest control point when mouse is clicked
        function Ind = closestPoint(obj)
            % 2x3 matlab matrix where mouse click is
            currpoint = obj.ax.CurrentPoint;
            % x and y coordinates of the mouse click
            coord = [currpoint(1,1); currpoint(1,2)];

            d = obj.cp - coord;

            
            [~, Ind] = min(sum(d.^2,1));
            
        end
        

        % Updates and draws control points
        function drawPoints(obj)
            obj.pointpatch.Vertices = obj.cp';
            obj.pointpatch.Faces = [1:size(obj.cp,2)]';
            %axis tight equal;
        end

        function drawSpline(obj)
            obj.lines.XData =  obj.cp(1,:);
            obj.lines.YData = obj.cp(2,:);
            
            t_samples = linspace(0,obj.spline.t_max, obj.num_samples);
            gamma = obj.spline.evaluate(t_samples);
            obj.splineplot.XData = gamma(1,:);
            obj.splineplot.YData = gamma(2,:);
            %axis tight equal;

        end

        function updateSpline(obj, cp)
            obj.spline.setControlPoints(cp);
            obj.cp = cp;

        end

        function drawStrip(obj)

            t_samples = linspace(0,obj.spline.t_max,obj.num_samples);
            gamma = obj.spline.evaluate(t_samples);
            kappa = obj.spline.curvature(t_samples);
            to = gamma(:,2:end)-gamma(:,1:end-1);
            seg_lens = sqrt(sum(to.^2,1));
            s = [0 cumsum(seg_lens)];

            [~, gamma_infl] = obj.spline.findInflectionPoints();
            lpopt = LPStiffnessOptimizer(gamma, kappa, gamma_infl);

            K = lpopt.optimizeWithInflections();

            if isempty(K)
                fprintf('Curve is infeasible!');
                return;
            end

            max_width = 0.05;
            half_width = 0.5*max_width/max(K)*K;
            p_simple = GeometryGenerator.generateProfileOutlineLoop(s,half_width,0.01);
            SvgTools.exportCurves('spline-curve/strip_simple.svg', {p_simple}, 1e3/0.352778);

            %obj.figstrip = patch(obj.axstrip);
            obj.strip_patch.Vertices = p_simple';
            obj.strip_patch.Faces = [1:size(p_simple,2)-1; 2:size(p_simple,2)]';
            
            figure(obj.figstrip); % make the strip figure current, so the next line affects that figure, and not the spline figure
            axis tight equal;

        end

       

    end
end

