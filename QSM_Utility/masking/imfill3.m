function imgfill = imfill3(img)
%% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
%
% use imfill to fill holes in 3D volume
% imfill in 3 directions and only fill in those holes for all directions

N = size(img);

imgXfill = zeros(N);
imgYfill = zeros(N);
imgZfill = zeros(N);

for ii = 1:N(1)
    imgXfill(ii,:,:) = imfill(squeeze(img(ii,:,:)), 'holes');
end

for jj = 1:N(2)
    imgYfill(:,jj,:) = imfill(squeeze(img(:,jj,:)), 'holes');
end

for kk = 1:N(3)
    imgZfill(:,:,kk) = imfill(squeeze(img(:,:,kk)), 'holes');
end

imgfill = (imgXfill + imgYfill + imgZfill) >= 3;