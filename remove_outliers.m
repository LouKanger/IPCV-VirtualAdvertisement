function [lines_vert, lines_hor] = remove_outliers(lines_vert, lines_hor)
%REMOVE_OUTLIERS removes the outliers in the vertical and horizontal line sets.
%This is done by combing the fact that the vertical lines are sorted such that
%their rho value from the Hough transform is ascending and descending for the
%horizontal lines. Inspecting the result shows that the theta values are then
%sorted in the opposite direction except for the outliers. Hence, these are
%filtered by looking at the increase and decrease in theta value of two
%consecutive detected lines for horizontal and vertical respectively.
%
% Author: L.W.J. Kanger, University of Twente
%
%   Parameters
%   ----------
%   lines_vert : struct
%       A Nx1 struct where N is the number of vertical lines. The struct has 4
%       fields: point1, point2, theta and rho. Point1 and point2 indicate the
%       start and end pixel coordinates of the found line [row,col]. Theta and
%       rho are the corresponding Hough transform coordinates, theta is the
%       angle and rho the distance from the origin.
%       Note: this should be the output of the GROUP_LINES function.
%   lines_hor : struct
%       A Nx1 struct where N is the number of horizontal lines. The struct has 4
%       fields: point1, point2, theta and rho. Point1 and point2 indicate the
%       start and end pixel coordinates of the found line [row,col]. Theta and
%       rho are the corresponding Hough transform coordinates, theta is the
%       angle and rho the distance from the origin.
%       Note: this should be the output of the GROUP_LINES function.
% 
%   Returns
%   -------
%   lines_vert : struct
%       A Mx1 struct where M is the remaining number of lines. The struct 
%       is built up in the same way as the original lines struct, but without
%       the outliers.
%   lines_hor : struct
%       A Mx1 struct where M is the remaining number of lines. The struct 
%       is built up in the same way as the original lines struct, but without
%       the outliers.
%

% Find outliers untill none are found
found_outlier = false;
while true
    % Assume no outliers
    found_outlier = false;
    
    % Loop over all vertical lines but the first
    for i = 2:size(lines_vert, 1)
        % If next line has larger angle then prev, it is a outlier. Remove it
        % reset the counter and examine the array again
        if lines_vert(i).theta > lines_vert(i-1).theta
            lines_vert(i) = [];
            found_outlier = true;
            break
        end
    end
    
    % Break loop if no outliers are found
    if found_outlier == false
        break
    end
end

% Find outliers untill none are found
found_outlier = false;
while true
    % Assume no outliers
    found_outlier = false;
    
    % Loop over all horizontal lines but the first
    for i = 2:size(lines_hor, 1)
        % If next line has smaller angle then prev, it is a outlier. Remove it
        % reset the counter and examine the array again
        if lines_hor(i).theta < lines_hor(i-1).theta
            lines_hor(i) = [];
            found_outlier = true;
            break
        end
    end
    
    % Break loop if no outliers are found
    if found_outlier == false
        break
    end
end
end

