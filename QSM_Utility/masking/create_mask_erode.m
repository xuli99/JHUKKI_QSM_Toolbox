function [maskErode,maskLargest,maskThreshold] = create_mask_erode(dataMagnitude,radiusDisk,percentileThresh)
% [maskErode,maskLargest,maskThreshold] = create_mask_erode(dataMagnitude,radiusDisk,percentileThresh)
%
%=================================
% PURPOSE
%=================================
%
% creates a new mask from magnitude volume by eroding outside of image with
% a structural disk
%
%=================================
% INPUTS
%=================================
%
% dataMagnitude     volume [nRows,nCols,nSlices]
%
% radiusDisk        radius of disk for structural element erosion; default = 6;
%
% percentileThresh  threshold for percentile of initial mask; default = 70;
%                   usually, values of 60 to 70 will suffice
%
%=================================
% OUTPUT
%=================================
%
% maskErode         binary image of eroded mask, based on threshold and disk
%                   erosion; [nRows,nCols,nSlices]
%
% maskLargest       binary image of the largest "region" after thresholding and
%                   cumulative summing to pick out the largest continuous region
%                   [nRows,nCols,nSlices]
%
% maskThreshold     original mask from thresholding [nRows,nCols,nSlices]
%
%=================================
% INFO
%=================================
% Author:          Jonathan Farrell, PhD; Issel Anne Lim
% Created:         July 13, 2010
% Last Modified:   2010-09-03
%
%=================================
% EXAMPLE
%=================================
%
% dataMagnitude = readraw('20100826az_Magnitude.img','float',[nRows,nCols,nSlices]);
% radiusDisk = 6; percentileThresh = 70;
% [maskErode,maskLargest,maskThreshold] = create_mask_erode(dataMagnitude,radiusDisk,percentileThresh);
%
%=================================
% BODY OF CODE
%=================================

if (nargin == 1)
    percentileThresh = 80;
    radiusDisk = 6;
end

disp('Creating mask...')
maskThreshold = zeros(size(dataMagnitude));
maskThreshold(dataMagnitude > prctile(dataMagnitude(:),percentileThresh)) = 1;

maskThreshold = imfill(maskThreshold,8,'holes');
% note: 4 = 2D four-connected neighborhood; 8 = 2D 8-connected
% neighborhood; 6 = 3D six-connected neighborhood;

% Choose the largest region
maskLargest = maskThreshold;
maskAll = cjlabel(maskThreshold,3);
for ii = 1:size(maskThreshold,3)
    viewMap = maskAll(:,:,ii).*dataMagnitude(:,:,ii);
    m0Slice = dataMagnitude(:,:,ii);
    nMean = mean(viewMap(:))./max(m0Slice(:));
    if nMean < 0.01
        maskLargest(:,:,ii) = maskAll(:,:,ii,2);
    end
end

% STREL creates morphological structuring element
structureElement = strel('disk',radiusDisk);
% Note: default "N" of strel is 4; erode based on this structureElement
maskLargest = imerode(maskLargest,structureElement);
 
% Choose largest region
maskLargest = cjlabel(maskLargest,1);

maskErode = maskLargest;

h = waitbar(0,'Eroding mask for each slice');
for ii = 1:size(maskLargest,3)   
    maskSlice = maskLargest(:,:,ii);
    maskErode(:,:,ii) = imopen(maskSlice,structureElement);
    waitbar(ii/size(maskLargest,3),h);
end
close(h);

maskErode = cjlabel(maskErode,1);