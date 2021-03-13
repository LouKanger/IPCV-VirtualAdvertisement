function [field_mask] = enhance_field_mask(field_mask)
%ENHANCE_FIELD_MASK enhances the field mask produced by CONSTRUCT_FIELD_MASK
% Author: L.W.J. Kanger, University of Twente
%
%   Parameters
%   ----------
%   field_mask : NxM matrix of type uint8
%       Matrix of ones and zeros indicating which part of the image to keep (1)
%       and which part of the image to discard (0) such that only the playing
%       field is left. This should be the output of CONSTRUCT_FIELD_MASK.
% 
%   Returns
%   -------
%   field_mask : NxM matrix of type uint8
%       Matrix of ones and zeros indicating which part of the image to keep (1)
%       and which part of the image to discard (0) such that only the playing
%       field is left.

% Enhance field mask using morphological operations
field_mask = bwareafilt(logical(field_mask),1); % keep largest area
field_mask = imfill(field_mask, 'holes');       % fill the remaining holes
field_mask = uint8(field_mask);                 % convert to binary image mask
end

