function [ CSFmask] =CSFmaskThresh(dataMap, thresh, maskBrain, voxSize)
%[ CSFmask ] =CSFmaskThresh(dataMap, thresh, maskBrain)
%   Get the central CSF mask based on threshold and maskBrain if available
%   Updated 2018-04-09, X.L. make it to handle bad R2* data
%  Ref: Liu et al., 2018 79(5):2795-2803

if nargin < 3
    maskBrain = ones(size(dataMap));
    voxSize = [1,1,1];
elseif nargin < 3
    voxSize = [1,1,1];
end

% maskflag = 1;       % default using maskflag

% Get initial data Mask by thresholding
maskBrain = maskBrain > 0;
dataMask = (dataMap < thresh).*maskBrain;     % inital mask
dataMask = dataMask > 0;
CSFmask = zeros(size(dataMask));

% Set size of the Central region
CenterRegionRadius = [30, 30, 30];            % in mm
CenterRegionRadiusVoxel = floor(CenterRegionRadius./voxSize);

% further trim based on connectivity and centroid
% get Centroid of Brain region
BrainCenter = regionprops(maskBrain, 'Centroid', 'Area');       % x, y, z

% in case maskBrain is discontinous and there are more than one center
if size(BrainCenter, 1) > 1
    totalArea = 0;
    temp = zeros(size(BrainCenter(1).Centroid));
    
    for centerii = 1:size(BrainCenter, 1)
        temp = temp + (BrainCenter(centerii).Centroid).*BrainCenter(centerii).Area;      % weighted sum
        totalArea = totalArea + BrainCenter(centerii).Area;
    end
    BrainCenter = [];
    BrainCenter.Centroid = temp./totalArea;
end

% get Central region mask
CenterRegion = zeros(size(dataMap));
CenterRegion(floor(BrainCenter.Centroid(2)-CenterRegionRadiusVoxel(2)) : floor(BrainCenter.Centroid(2)+CenterRegionRadiusVoxel(2)), ...
                floor(BrainCenter.Centroid(1)-CenterRegionRadiusVoxel(1)) : floor(BrainCenter.Centroid(1)+CenterRegionRadiusVoxel(1)),...
                floor(BrainCenter.Centroid(3)-CenterRegionRadiusVoxel(3)) : floor(BrainCenter.Centroid(3)+CenterRegionRadiusVoxel(3))) = 1;

% get inital CSF mask in the Center (overlap between dataMask and CenterRegion)
CSFmask0 = CenterRegion & dataMask;            
cc6mask = bwconncomp(CSFmask0, 6);
numPixels = cellfun(@numel, cc6mask.PixelIdxList);
[~, I] = sort(numPixels, 'descend');

% Pick the largest regions and combine
CSFmask0 = zeros(size(dataMask));
numRegionComb = 3;      % 2 or 3
for iiRegion = 1:numRegionComb
    CSFmask0(cc6mask.PixelIdxList{I(iiRegion)}) = 1;
end

% Extend to connected regions (use overlap of 10)
cc6mask = bwconncomp(dataMask, 6);
numPixels = cellfun(@numel, cc6mask.PixelIdxList);
[~, I] = sort(numPixels, 'descend');

for iiRegion = 1:cc6mask.NumObjects
    if sum(CSFmask0(cc6mask.PixelIdxList{I(iiRegion)})) > 0
        CSFmask(cc6mask.PixelIdxList{I(iiRegion)}) = 1;
    end
end

end

