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
        spline % spline curve object

        pointpatch % 
        mouseflag
    end
    
    methods
         function obj = UserInterface(cp)
            obj.fig = figure('WindowButtonDownFcn', @(source, event) obj.mouseDown(source, event));
            obj.ax = axes(obj.fig);
            xlim([0,1]);
            ylim([0,1]);
            obj.cp = cp;
            obj.pointpatch = patch('Vertices', cp', 'Faces', [1:size(cp,2)]', 'Parent', obj.ax, 'Marker', 'o');
            obj.mouseflag = false;
        end

        function mouseDown(obj, source, event)
            % called when the mouse is pushed down
            % left click for move and add, maybe hold down shift for add
            % right click for deletion
            % figure property : selection type 
            if(strcmp(obj.fig.SelectionType, 'normal')) % moving points
                obj.mouseflag = true;

            elseif(strcmp(obj.fig.SelectionType, 'extend')) % adding points
                % 2x3 matlab matrix where mouse click is
                currpoint = obj.ax.CurrentPoint;
                % x and y coordinates of the mouse click
                coord = [currpoint(1,1); currpoint(1,2)];

                % appending the coordinates to the end of the cp matrix
                obj.cp = [obj.cp, coord];
            elseif(strcmp(obj.fig.SelectionType, 'alt')) % deleting points


            end

            currpoint = obj.ax.CurrentPoint;
            fprintf('Current point: %.5f %.5f \n', currpoint(1,1), currpoint(1,2));
            obj.cp(:,1) = currpoint(1, 1:2)';
            obj.drawPoints();
            
        end
        
        function mouseUp(obj, source, event)
            % called when the mouse is released
            obj.mouseflag = false;
        end

        function mouseMove(obj, source, event)
            % called when the mouse is pressed down and moving
            
        end

        function closestPoint(obj)
            % 2x3 matlab matrix where mouse click is
            currpoint = obj.ax.CurrentPoint;
            % x and y coordinates of the mouse click
            coord = [currpoint(1,1); currpoint(1,2)];
            
            [~, Ind] = min()
            
        end

        function drawPoints(obj)
            obj.pointpatch.Vertices = obj.cp';
            obj.pointpatch.Faces = [1:size(obj.cp,2)]';
        end

        function drawSpline(obj)
        
        end

       

    end
end

