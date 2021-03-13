function [lines] = remove_duplicate_lines(lines, dtheta, drho)
%REMOVE_DUPLICATE_LINES Removes duplicate lines based on the angle difference
%and distance to the origin difference of two lines. NOTE: the input should be
%the output of the DETEC_LINES function to ensure the correct format.
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
%   dtheta : double
%       Value that indicates how much the angle between two lines can be to
%       still be considered two different lines.
%   drho : double
%       Value that indicates how much the distance between two lines can be to
%       still be considered two different lines.
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

% Loop over all the lines and check if duplicate; remove the shortest one
remove_idx = [];
for i = 1:length(lines)
    line1 = lines(i);
    for j = 1:length(lines)
        line2 = lines(j);
        % Consider only the lines that are nearly parallel
        if j ~= i && abs(line1.theta - line2.theta) < dtheta
            % if rho1 is close to rho2, line1 is on same field line as line2
            if abs(line1.rho - line2.rho) < drho
                % calculate length of both lines
                d1 = norm(line1.point1 - line1.point2);
                d2 = norm(line2.point1 - line2.point2);
                
                % store longest line (=more accurate)
                if d1 < d2
                    remove_idx = [remove_idx, i];
                else
                    remove_idx = [remove_idx, j];
                end
            end
        end
    end
end

% Remove the duplicate lines
lines(unique(remove_idx)) = [];

end

