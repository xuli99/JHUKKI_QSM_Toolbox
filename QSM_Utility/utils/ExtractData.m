function [xMean, xStd, Vol, xMedian]  = ExtractData(x, ROImap, roiIndex, erodesize, maskBrain)
% [xMean, xStd, Vol, xMedian]  = ExtractData(x, ROImap, roiIndex, erodesize, maskBrain)
%% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
% 
% Extract data statistics using ROImap and roiIndex
% Assumed ROImap and x are matched
% roiIndex is nxm, n: different ROIs, m index number corresponding to
% indexes in ROImap
% 
% Example use
% [chiMean(iiSub, :, iidatasetTotal), chiStd(iiSub, :, iidatasetTotal)] = ExtractData(DispData, parcelAllDat, roiIndex);
% [chiMean(iiSub, :, iidatasetTotal), chiStd(iiSub, :, iidatasetTotal)] = ExtractData(DispData, parcelAllDat, roiIndex, 1, maskErode);
% 
% Updated 2024-07-12, X.L.,added xMedian output

%%
% Basic check
if nargin < 3
    error('not enough input ...')
elseif nargin < 4
    erodesize = 0;      % default no erode
    maskBrain = ones(size(x));
elseif nargin < 5
    maskBrain = ones(size(x));
end

if size(x) ~= size(ROImap) 
    error('ROImap does not match the data itself ...')
end

se = strel('disk', erodesize);

%%
[nROI, nIndex] = size(roiIndex);       % roiIndex tells in which part of the brain to extract data   

xMean = zeros(nROI, 1);
xStd = zeros(nROI, 1);
Vol = zeros(nROI, 1);

xMedian = zeros(nROI, 1);

for iiROI = 1:nROI

    maskROI = zeros(size(ROImap));

    % making the maskROI
    for roiIndexii = 1:nIndex
        maskROI = maskROI | (ROImap == roiIndex(iiROI, roiIndexii));
    end       
    maskROI = maskROI > 0;
    
    % mask erosion and correction
    if erodesize > 0
        maskROI = imerode3dslice(maskROI, se);
    end            
    maskROI = (maskROI.*maskBrain) > 0;         % exclude regions outside of brain     
    
    % doing statistics
    temp = x(maskROI);
    xMean(iiROI) = mean(temp(:));
    xStd(iiROI) = std(temp(:));
    Vol(iiROI) = length(temp(:));

    xMedian(iiROI) = median(temp(:));      % adding xMedian
 
end
