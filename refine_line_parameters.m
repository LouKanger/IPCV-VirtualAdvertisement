function [hlines1, hlines2] = refine_line_parameters(BW, hlines1, hlines2, sr, num_iters)
%REFINE_LINE_PARAMETERS refines the line parameters of two sets of lines of a
%binary image. The set of lines should be in homogeneous coordinates, but they
%do not have to have the same number of lines.
% Author: L.W.J. Kanger, University of Twente
%
%   Parameters
%   ----------
%   BW : uint8
%       A RxC matrix that contains the edges on which the lines are detected.
%   hlines1 : double
%       A Nx3 matrix that contains the lines in the homogeneous coordinate
%       system.
%   hlines2 : double
%       A Nx3 matrix that contains the lines in the homogeneous coordinate
%       system.
%   sr : double
%       Specifies the maximum distance bewteen the edge pixels and the lines
%   num_iters : int
%       Specifies the number of iterations that are performed to refine the
%       parameters of the lines
% 
%   Returns
%   -------
%   hlines1 : double
%       A Nx3 matrix that contains the refined lines in the homogeneous 
%       coordinate system.
%   hlines2 : double
%       A Nx3 matrix that contains the refined lines in the homogeneous 
%       coordinate system.
%

% perform the line refinement
for iter = 1 : num_iters
    for k = 1 : max( [size(hlines1,1), size(hlines2,1)] )
        L1 = [];
        L2 = [];
        for i = 1:size(BW, 1)
            for j = 1:size(BW, 2)
                if BW(i,j) == 1
                    y = i; x = j;
                    % find white pixels for first line set refinement
                    if k <= size(hlines1,1)
                        if abs( dot(hlines1(k,:), [x, y, 1]) ) <= sr
                            L1 = [L1; [x, y]];
                        end
                    end
                    % find white pixels for the second line set refinement
                    if k <= size(hlines2,1)
                        if abs( dot(hlines2(k,:), [x, y, 1]) ) <= sr
                            L2 = [L2; [x, y]];
                        end
                    end
                end
            end
        end
        
        % refinement of the first set of lines
        if k <= size(hlines1,1)
            fit1 = fit(L1(:,1), L1(:,2), 'poly1'); %,'Robust','LAR');
            
            % Construct the new line in homogeneous coordinates
            my = 1 / fit1.p2;
            mx = -fit1.p1 * my;
            d = 1 / sqrt(mx^2 + my^2);
            hlines1(k,:) = [mx * d, my * d, -d];
        end
        
        % refinement of the second set of lines
        if k <= size(hlines2,1)
            fit2 = fit(L2(:,1), L2(:,2), 'poly1');   %,'Robust','LAR');
            
            % Construct the new line in homogeneous coordinates
            my = 1 / fit2.p2;
            mx = -fit2.p1 * my;
            d = 1 / sqrt(mx^2 + my^2);
            hlines2(k,:) = [mx * d, my * d, -d];
        end
    end
end
end

