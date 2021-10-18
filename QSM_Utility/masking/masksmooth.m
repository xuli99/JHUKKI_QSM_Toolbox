function maskOut = masksmooth(maskIn, r)
% function newMask = masksmooth(oldMask)

if nargin < 2
    r = 1;
end

% do some mask smoothing
se = strel('sphere', r);
maskOut = imclose(maskIn, se);