function [white_mask] = detect_white_pixels(I, tl, td, tau)
%DETECT_WHITE_PIXELS Detects white line pixels in an image based on the given
%threshold values for the luminance and the width of the white lines.
% Author: L.W.J. Kanger, University of Twente
%
%Implementation of the algorithm that is formulated in the article by Farin, 
%Dirk & Krabbe, Susanne & With, Peter & Effelsberg, Wolfgang. (2004). Robust 
%camera calibration for sport videos using court models. 5307. 80-91. 
%DOI: 10.1117/12.526813. 
%
%   Parameters
%   ----------
%   I : uint8
%       An NxMx3 image used to find white line pixels. It should be in RGB
%       format otherwise the conversion to YCbCr space does not go correctly.
%   tl : int
%       An integer in the range 0 to 255 that acts as the luminance threshold.
%       Only pixels that are bright enough are white line pixel candidates.
%   td : int
%       A positive integer (can be zero) that acts as the difference threshold
%       between a white line pixel candidate and the pixel tau places above,
%       below, left or right of it.
%   tau : int
%       A positive integer that specifies half the line width. 
% 
%   Returns
%   -------
%   white_mask : uint8
%       An NxM matrix of ones and zeros indicating which parts of the image have
%       white line pixels (1) and which not (0) such that the non white line
%       pixels can be filtered out.
%

% Convert RGB image to YCbCr image and only take the Y (luminace) part
im_masked_YCbCr = rgb2ycbcr(I);
g = im_masked_YCbCr(:,:,1);

% loop over the image and find the white pixels
white_mask = zeros(size(g));
for i = (1+tau):(size(g, 1)-tau)
    for j = (1+tau):(size(g, 2)-tau)
        % check if bright enough
        if (g(i,j) >= tl)
            % check vertically for dark enough pixel values
            if (g(i,j)-g(i-tau,j) > td) && (g(i,j) - g(i+tau,j) > td)
                white_mask(i,j) = 1;
            % check horizontally for dark enough pixel values
            elseif (g(i,j)-g(i,j-tau) > td) && (g(i,j) - g(i,j+tau) > td)
                white_mask(i,j) = 1;
            end
        end
    end
end

white_mask = uint8(white_mask);
end

