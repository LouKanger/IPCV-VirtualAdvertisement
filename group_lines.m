function [lines_vert, lines_hor] = group_lines(lines)
%GROUP_LINES groups the lines to vertical or horizontal based on their angle in
%the image. The resulting groups are sorted in such a way that the vertical
%lines are sorted to appear from left to right in the image. For the horizontal
%lines this sorting is done as well but then from top to bottom.
%
% Author: L.W.J. Kanger, University of Twente
%
%   Parameters
%   ----------
%   lines : struct
%       A Nx1 struct where N is the number of lines found. The struct has 4
%       fields: point1, point2, theta and rho. Point1 and point2 indicate the
%       start and end pixel coordinates of the found line [row,col]. Theta and
%       rho are the corresponding Hough transform coordinates, theta is the
%       angle and rho the distance from the origin. The lines are sorted based
%       on the theta value in increasing order. 
%       Note: this should be the output of the HOUGH_LINE_DETECTION function.
% 
%   Returns
%   -------
%   lines_vert : struct
%       A Nx1 struct where N is the number of vertical lines found. The struct 
%       is built up in the same way as the original lines struct.
%   lines_hor : struct
%       A Nx1 struct where N is the number of horizontal lines found. The struct 
%       is built up in the same way as the original lines struct.
%

% Group the lines based on their angle in the image
lines_vert = [];
lines_hor = [];
for i = 1:length(lines)
    % Horizontal line if angle is negative
    if lines(i).theta < 0
        lines_hor = [lines_hor, lines(i)];
    else
        lines_vert = [lines_vert, lines(i)];
    end
end

% Order the vertical lines left to right and horizontal lines top to bottom
lines_vert = table2struct(sortrows(struct2table(lines_vert), 'rho'));
lines_hor = table2struct(sortrows(struct2table(lines_hor), 'rho', 'descend'));
end



