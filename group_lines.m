function [lines_vert, lines_hor] = group_lines(lines, tmax)
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
%   tmax : double
%       Value that indicates when a line is considered to be horizontal in the
%       image. All lines that have an angle in the image less than tmax are
%       considered to be horizontal lines.
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
    % horizontal line if angle in image is less than x degrees
    if abs(90 - abs(lines(i).theta)) < tmax
        lines_hor = [lines_hor, lines(i)];
    else
        lines_vert = [lines_vert, lines(i)];
    end
end

% Order the vertical lines left to right and horizontal lines top to bottom
lines_vert = table2struct(sortrows(struct2table(lines_vert), 'rho'));
lines_hor = table2struct(sortrows(struct2table(lines_hor), 'rho', 'descend'));
end



