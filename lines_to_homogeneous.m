function [hlines] = lines_to_homogeneous(lines)
%LINES_TO_HOMOGENEOUS transforms a set of lines that are in 2D Euclidean space
%to the homogeneous space.
% Author: L.W.J. Kanger, University of Twente
%
%   Parameters
%   ----------
%   lines : struct
%       A Nx1 struct where N is the number of lines. The struct has 4
%       fields: point1, point2, theta and rho. Point1 and point2 indicate the
%       start and end pixel coordinates of the found line [row,col]. Theta and
%       rho are the corresponding Hough transform coordinates, theta is the
%       angle and rho the distance from the origin. The lines are sorted based
%       on the theta value in increasing order. 
%       Note: this input should be of the structure as the output of the 
%       HOUGH_LINE_DETECTION function to ensure correct functionality.
% 
%   Returns
%   -------
%   lines_h : double
%       A Nx3 matrix where N is the number of lines and 3 are the three
%       dimensions of the homogeneous space.
%

% Convert to homogeneous coordinates [x,y] -> [cx, cy, c]
hlines = zeros(length(lines),3);
for k = 1:length(lines)
    hlines(k,:) = [cosd(lines(k).theta), sind(lines(k).theta), -lines(k).rho];
end

end

