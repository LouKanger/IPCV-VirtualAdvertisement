function [field_mask] = construct_field_mask(I, num)
%CONSTRUCT_FIELD_MASK Constructs a binary field mask to only consider the 
% playing field and not the surroundings such as the boardings and part of the 
% stadium.
% Author: L.W.J. Kanger, University of Twente
%
%   Parameters
%   ----------
%   I : uint8
%       An NxMx3 image to construct a binary field mask of. The image must be in
%       RGB format such that the conversion to HSV goes correctly.
%   num : int
%       An integer that indicates how many color points around the dominant
%       color peak should be conserved.
% 
%   Returns
%   -------
%   field_mask : NxM matrix of type uint8
%       Matrix of ones and zeros indicating which part of the image to keep (1)
%       and which part of the image to discard (0) such that only the playing
%       field is left.
%

% Convert to HSV space and create a historgram on the H values
im_hsv = rgb2hsv(I);
[counts, x] = imhist(im_hsv(:,:,1));

% Based on colour histogram, we want only those pixels at dominant peak
peak_idx = find(counts == max(counts));

% Take the points +- num from the peak (so total is 2*num + 1)
lambdas = x((peak_idx-num):(peak_idx+num));

% Apply thresholding on the HSV image to get the field mask
field_mask = (im_hsv(:,:,1) > min(lambdas)) .* (im_hsv(:,:,1) < max(lambdas));
end

