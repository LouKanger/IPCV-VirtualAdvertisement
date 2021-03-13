function [lines] = hough_line_detection(BW, num_peaks, fill_gap, min_length)
%HOUGH_LINE_DETECTION Detects and contructs lines in the given binary image using a
%Hough transformation. The lines detected are dependent on the number of Hough
%peaks, the maximum gap between to line segments and the minimum length of the
%found lines.
% Author: L.W.J. Kanger, University of Twente
%
%   Parameters
%   ----------
%   BW : uint8
%       An NxM image on which a Hough transform is applied to detect lines.
%   num_peaks : int
%       An integer that indicates how many Hough peaks should be used, which
%       affects how many lines are found.
%   fill_gap : int
%       An integer that indicates what the maximum gap between to line segements
%       can be to still be considered the same line.
%   min_length : int
%       An integer that indicates what the minimum length of the found lines
%       should be. Any line shorter than this value will not be returned.
% 
%   Returns
%   -------
%   lines : struct
%       A Nx1 struct where N is the number of lines found. The struct has 4
%       fields: point1, point2, theta and rho. Point1 and point2 indicate the
%       start and end pixel coordinates of the found line [row,col]. Theta and
%       rho are the corresponding Hough transform coordinates, theta is the
%       angle and rho the distance from the origin. The lines are sorted based
%       on the theta value in increasing order.
%

% Perform a Hough transform on the binary lines/edges
[H, theta, rho] = hough(BW);

% Find the peaks of the Hough transform
hough_peaks  = houghpeaks(H, num_peaks, 'threshold', 0);

% Construct the lines from the found peaks
lines = houghlines(BW, theta, rho, hough_peaks,...
    'FillGap', fill_gap, 'MinLength' ,min_length);

% sort the lines struct based on the theta value (increasing order)
lines = table2struct(sortrows(struct2table(lines), 'theta'));
end

