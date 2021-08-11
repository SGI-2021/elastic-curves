classdef UserInterface < handle
    %USERINTERFACE Summary of this class goes here
    %   Detailed explanation goes here
    % Method for each mouse event (click, drag, up)
    % How to distiguish between adding new points and dragging around
    % doc figure -> see also -> figure properties for documentations
    
    properties
        fig % figure opened
        ax % axes

        cp  % control points
        lines %lines connecting control points
        splineplot % drawing of spline

        spline % spline curve object
        num_samples % number of samples taken

        gamma
        kappa % curvature
        s

        pointpatch % 
        mouseflag
    end
    
    methods
         function obj = UserInterface(cp)
            obj.fig = figure('WindowButtonDownFcn', @(source, event) obj.mouseDown(source),...
                             'WindowButtonUpFcn', @(source, event) obj.mouseUp(source, event), ...
                             'WindowButtonMotionFcn', @(source, event) obj.mouseMove(source, event));

            %obj.fig = figure('WindowButtonDownFcn', @(source, event) obj.mouseDown(source, event))
                             %'WindowButtonMotionFcn', @(source, event) obj.mouseMove(source, event));
                             %'WindowButtonUpFcn', @(source, event) obj.mouseUp(source, event));
            obj.ax = axes(obj.fig);
            xlim([0,1]);
            ylim([0,1]);
            obj.cp = cp;

            obj.mouseflag = false;

            hold on;
            obj.pointpatch = patch('Vertices', cp', ...
                                   'Faces', [1:size(cp,2)]', ...
                                   'Parent', obj.ax, ...
                                   'Marker', '.');
            obj.lines = plot(obj.cp(1,:), obj.cp(2,:),'k--');

            obj.spline = SplineCurve.import('rect_spline1.txt');
            scale = 0.08;
            obj.spline.cp = obj.spline.cp * scale;
            obj.num_samples = 400;

            % Sample the curve, and approximate the arc-length parameter s at every
            % sample
            t_samples = linspace(0,obj.spline.t_max, obj.num_samples);
            obj.gamma = obj.spline.evaluate(t_samples);
            obj.kappa = obj.spline.curvature(t_samples);
            to = obj.gamma(:,2:end)-obj.gamma(:,1:end-1);
            seg_lens = sqrt(sum(to.^2,1));
            obj.s = [0 cumsum(seg_lens)];

            obj.splineplot = plot(obj.gamma(1,:), obj.gamma(2,:),'LineWidth',2,'Color',[1 0 0]);

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
 
            % ADDING POINTS (SHIFT CLICK):
            elseif(strcmp(source.SelectionType, 'extend')) 
                % Checks if a shift-click is made near a control point
                % Adds another control point in between  the control point
                % that was shift-clicked and the next one. If the last 
                % control point is shift-clicked, adds a control point in
                % between the last and second from last control points             

                if(isempty(obj.cp) || size(obj.cp,2) == 1)
                    % If there are no control points add a control point
                    % where the cursors is
                    currpoint = obj.ax.CurrentPoint;
                    coord = [currpoint(1,1); currpoint(1,2)];

                    obj.cp = coord;
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
                    obj.cp = [c1 c2 c3];
                else
                    midpoint = ((obj.cp(:,end) + obj.cp(:,end-1))'./2)';
                    obj.cp = [obj.cp(:,1:end-1), midpoint, obj.cp(:,end)];
                end
                
                % Draws contol point
                obj.drawPoints();
                %obj.drawSpline();
                               
            % DELETING POINTS (CONTROL-CLICK)
            elseif(strcmp(source.SelectionType, 'alt'))
                % Checks for a control-click
                % Deletes the control point that is control-clicked
                
                closestInd = obj.closestPoint();

                obj.cp(:,closestInd) = [];

                obj.drawPoints();
                %obj.drawSpline();
                
            end
        end

     
        % Updates mouseflag when mouse is not being pressed
        function mouseUp(obj, source, event)
            obj.mouseflag = false;
            obj.spline.setControlPoints(obj.cp);
            obj.drawSpline();
        end

        % Moves a Control Point 
        function mouseMove(obj, source, event)
            closestInd = obj.closestPoint();
            if(obj.mouseflag == true)
                 % Activates when the mouse is pressed down and moving
                 currpoint = obj.ax.CurrentPoint;
                 coords = [currpoint(1,1); currpoint(1,2)];
                 obj.cp(:, closestInd) = coords;
                 obj.drawPoints();
                 
            end
        
        end
        
        % Finds the closest control point when mouse is clicked
        function Ind = closestPoint(obj)
            % 2x3 matlab matrix where mouse click is
            currpoint = obj.ax.CurrentPoint;
            % x and y coordinates of the mouse click
            coord = [currpoint(1,1); currpoint(1,2)];

            d1 = (obj.cp(1,:)' - repmat(coord(1), size(obj.cp,2), 1));
            d2 = (obj.cp(2,:)' - repmat(coord(2), size(obj.cp,2), 1));
            d = [d1,d2];

            [~, Ind] = min(normrow(d), [], 'linear');
            
        end
        
        % Updates and draws control points
        function drawPoints(obj)
            obj.pointpatch.Vertices = obj.cp';
            obj.pointpatch.Faces = [1:size(obj.cp,2)]';
        end

        function drawSpline(obj)
            obj.lines = plot(obj.cp(1,:), obj.cp(2,:),'k--');
            %obj.splineplot = plot(obj.gamma(1,:), obj.gamma(2,:),'LineWidth',2,'Color',[1 0 0]);

        end

       

    end
end

