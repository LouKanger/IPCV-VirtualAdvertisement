function plot_hlines(lines_h, im, color)
%UNTITLED2 Summary of this function goes here
%   INPUT
%   width   :   width of the soccer field in meters
%   height  :   height of the soccer field in meters
%
%   OUTPUT
%   field_lines : Nx3 array of N lines in the form ax+by+c=0 --> [a,b,c]
%   field_points : all real intersection points (no virtual intersections)
%
% Author: L.W.J. Kanger, University of Twente

for k = 1:size(lines_h, 1)
    % lines_h(k) = [nx, ny, -d] but also [a, b, c] --> ax + by + c = 0
    x1 = 0;
    y1 = (-lines_h(k,3) - lines_h(k,1) * x1) / lines_h(k, 2);
    if y1 < 0
        y1 = 0;
        x1 = (-lines_h(k,3) - lines_h(k,2) * y1) / lines_h(k, 1);
    elseif y1 > size(im, 1)
        y1 = size(im, 1);
        x1 = (-lines_h(k,3) - lines_h(k,2) * y1) / lines_h(k, 1);
    end
    
    x2 = size(im, 2);
    y2 = (-lines_h(k,3) - lines_h(k,1) * x2) / lines_h(k, 2);
    if y2 < 0
        y2 = 0;
        x2 = (-lines_h(k,3) - lines_h(k,2) * y2) / lines_h(k, 1);
    elseif y2 > size(im, 1)
        y2 = size(im, 1);
        x2 = (-lines_h(k,3) - lines_h(k,2) * y2) / lines_h(k, 1);
    end
    
    plot([x1, x2], [y1, y2],'LineWidth',2,'Color',color);
end
end

