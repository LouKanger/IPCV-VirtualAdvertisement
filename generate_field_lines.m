function [field_lines, field_points] = generate_field_lines(width, height)
%GENERATE_FIELD_LINES Generates the field lines and points of a given soccer field
% dimension in meters
% 
%   INPUT
%   width   :   width of the soccer field in meters
%   height  :   height of the soccer field in meters
%
%   OUTPUT
%   field_lines : Nx3 array of N lines in the form ax+by+c=0 --> [a,b,c]
%   field_points : all real intersection points (no virtual intersections)
%
% Author: L.W.J. Kanger, University of Twente

% general dimensions of lines (ref. https://the18.com/soccer-learning/soccer-field-dimensions-markings-explained-maximum-minimium-size-length-width)
center_circ_pos = [width/2, height/2];
center_circ_rad = 9.15;
goal_width = 7.32;  % distance from one pole to the other pole
goal_heigh = 2.44;  % 3d height (not used in this function)
goal_area = 5.5;
pen_area = 16.5;


% define all the straight lines
side_line_top = [0, 0; width, 0];
side_line_bot = [0, height; width, height];
goal_line_left = [0, 0; 0, height];
goal_line_right = [width, 0; width, height];
mid_line = [width/2, 0; width/2, height];

% NOTE: following lines have wrong order [x1, x2; y1, y2]. added transpose as fix
goal_area_1_top = [0, goal_area; 
    height/2 - goal_width/2 - goal_area, height/2 - goal_width/2 - goal_area]';
goal_area_1_bot = [0, goal_area;
    height/2 + goal_width/2 + goal_area, height/2 + goal_width/2 + goal_area]';
goal_area_1_side = [goal_area, goal_area; 
    height/2 - goal_width/2 - goal_area, height/2 + goal_width/2 + goal_area]';

pen_area_1_top = [0, pen_area; 
    height/2 - goal_width/2 - pen_area, height/2 - goal_width/2 - pen_area]';
pen_area_1_bot = [0, pen_area;
    height/2 + goal_width/2 + pen_area, height/2 + goal_width/2 + pen_area]';
pen_area_1_side = [pen_area, pen_area; 
    height/2 - goal_width/2 - pen_area, height/2 + goal_width/2 + pen_area]';

goal_area_2_top = [width, width - goal_area;
    height/2 - goal_width/2 - goal_area, height/2 - goal_width/2 - goal_area]';
goal_area_2_bot = [width, width - goal_area;
    height/2 + goal_width/2 + goal_area, height/2 + goal_width/2 + goal_area]';
goal_area_2_side = [width - goal_area, width - goal_area;
    height/2 - goal_width/2 - goal_area, height/2 + goal_width/2 + goal_area]';

pen_area_2_top = [width, width - pen_area;
    height/2 - goal_width/2 - pen_area, height/2 - goal_width/2 - pen_area]';
pen_area_2_bot = [width, width - pen_area;
    height/2 + goal_width/2 + pen_area, height/2 + goal_width/2 + pen_area]';
pen_area_2_side = [width - pen_area, width - pen_area;
    height/2 - goal_width/2 - pen_area, height/2 + goal_width/2 + pen_area]';

% store the points of the lines in a specific struct format (vertical then hor)
% vertical is sorted left to right, horizontal is sorted top to bottom (as seen
% in image)
point1 = [
    % vertical lines (sorted left to right left goal and mid)
    goal_line_left(1,:); goal_area_1_side(1,:); pen_area_1_side(1,:);
    mid_line(1,:); 
    
    % vertical lines (sorted left to right right goal only)
    pen_area_2_side(1,:); goal_area_2_side(1,:);
    goal_line_right(1,:);
    
    % horizontal lines (sorted top to bot left goal side)
    side_line_top(1,:); pen_area_1_top(1,:); goal_area_1_top(1,:);
    goal_area_1_bot(1,:); pen_area_1_bot(1,:); side_line_bot(1,:);
    
    % horizontal lines (sorted top to bot right goal side)
    goal_area_2_top(1,:); goal_area_2_bot(1,:); pen_area_2_top(1,:);
    pen_area_2_bot(1,:)];

point2 = [
    % vertical lines (sorted left to right left goal and mid)
    goal_line_left(2,:); goal_area_1_side(2,:); pen_area_1_side(2,:);
    mid_line(2,:); 
    
    % vertical lines (sorted left to right right goal only)
    pen_area_2_side(2,:); goal_area_2_side(2,:);
    goal_line_right(2,:);
    
    % horizontal lines (sorted top to bot left goal side)
    side_line_top(2,:); pen_area_1_top(2,:); goal_area_1_top(2,:);
    goal_area_1_bot(2,:); pen_area_1_bot(2,:); side_line_bot(2,:);
    
    % horizontal lines (sorted top to bot right goal side)
    goal_area_2_top(2,:); goal_area_2_bot(2,:); pen_area_2_top(2,:);
    pen_area_2_bot(2,:)];

% Convert the points to lines in homogeneous coordinate system
field_lines = zeros(length(point1), 3);
for i = 1:length(point1)
    field_lines(i,:) = cross([point1(i,:), 1.0], [point2(i,:), 1.0]);
end

field_points = table2struct(table(point1, point2));
end

% random
% point1 = [side_line_top(1,:); side_line_bot(1,:); goal_line_left(1,:); 
%     goal_line_right(1,:); mid_line(1,:); goal_area_1_top(1,:);
%     goal_area_1_bot(1,:); goal_area_1_side(1,:); pen_area_1_top(1,:);
%     pen_area_1_bot(1,:); pen_area_1_side(1,:); goal_area_2_top(1,:);
%     goal_area_2_bot(1,:); goal_area_2_side(1,:); pen_area_2_top(1,:);
%     pen_area_2_bot(1,:); pen_area_2_side(1,:);];
% 
% point2 = [side_line_top(2,:); side_line_bot(2,:); goal_line_left(2,:); 
%     goal_line_right(2,:); mid_line(2,:); goal_area_1_top(2,:);
%     goal_area_1_bot(2,:); goal_area_1_side(2,:); pen_area_1_top(2,:);
%     pen_area_1_bot(2,:); pen_area_1_side(2,:); goal_area_2_top(2,:);
%     goal_area_2_bot(2,:); goal_area_2_side(2,:); pen_area_2_top(2,:);
%     pen_area_2_bot(2,:); pen_area_2_side(2,:);];


% sorted vertical and horizontal
% point1 = [goal_line_left(1,:); goal_line_right(1,:); mid_line(1,:); 
%     goal_area_1_side(1,:); pen_area_1_side(1,:); goal_area_2_side(1,:);
%     pen_area_2_side(1,:);
%     side_line_top(1,:); side_line_bot(1,:); goal_area_1_top(1,:);
%     goal_area_1_bot(1,:);  pen_area_1_top(1,:); pen_area_1_bot(1,:); 
%     goal_area_2_top(1,:); goal_area_2_bot(1,:); pen_area_2_top(1,:);
%     pen_area_2_bot(1,:)];
% 
% point2 = [goal_line_left(2,:); goal_line_right(2,:); mid_line(2,:); 
%     goal_area_1_side(2,:); pen_area_1_side(2,:); goal_area_2_side(2,:);
%     pen_area_2_side(2,:);
%     side_line_top(2,:); side_line_bot(2,:); goal_area_1_top(2,:);
%     goal_area_1_bot(2,:);  pen_area_1_top(2,:); pen_area_1_bot(2,:); 
%     goal_area_2_top(2,:); goal_area_2_bot(2,:); pen_area_2_top(2,:);
%     pen_area_2_bot(2,:)];










