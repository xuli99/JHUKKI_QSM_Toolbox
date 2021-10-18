function wG = gradient_mask_all(rawImage, Mask, percentage, kspaceFlag, wG_thresh)
%  wG = gradient_mask_all(rawImage, Mask, percentage, wG_thresh)
% Generate weighting using gradient of rawImage
% with percentage of voxels in Mask to be edges
%   
% Authors: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
% Updated 2020-06

if nargin < 4
    kspaceFlag = 0; % default in image sapce
    wG_thresh = [];
elseif nargin < 5
    wG_thresh = [];
end

rawImage = rawImage.*Mask;

if kspaceFlag == 1      % calculated gradient in k-space
    N = size(rawImage);
    datatype = class(rawImage);
    
    % gradient operator in k-space
    [k2,k1,k3] = meshgrid(0:N(2)-1, 0:N(1)-1, 0:N(3)-1);
    k1 = cast(k1, datatype);
    k2 = cast(k2, datatype);
    k3 = cast(k3, datatype);
    
    fd1 = -1 + exp(2*pi*1i*k1/N(1));    % dimention 1
    fd2 = -1 + exp(2*pi*1i*k2/N(2));    % dimention 2
    fd3 = -1 + exp(2*pi*1i*k3/N(3));    % dimention 3

    rawImagek = fftn(rawImage);

    % gradient of initial estimation: 4D with gradx, grady, gradz
    rawImageGrad = real(cat(4, ifftn(rawImagek.*fd1), ifftn(rawImagek.*fd2), ifftn(rawImagek.*fd3)));  
else
    TV = TVOP;                  % total varition class (3D finite difference operator)
    rawImageGrad = TV*rawImage;
end

wG = abs(rawImageGrad);             % gradient norm

% initial guess of the threshold
wGtemp = sqrt(sum(wG.^2, 4));       % average over dimentions
if isempty(wG_thresh)
    wG_thresh = prctile(wGtemp(Mask>0), percentage*100)+1e-6;  % threshold based on gradient norm, all positive, dealt with 0 initial threshold
end

denominator = sum(Mask(:)>0);          % voxel number in brain mask
numerator = sum(wG(:)>wG_thresh);

if  (numerator/denominator)>3*percentage       % in toal 3*percentage
    while (numerator/denominator)>3*percentage
        wG_thresh = wG_thresh*1.05;
        numerator = sum(wG(:)>wG_thresh);
    end
else
    while (numerator/denominator)<3*percentage
        wG_thresh = wG_thresh*.95;
        numerator = sum(wG(:)>wG_thresh);
    end
end

% TV weighting
wG = (wG<=wG_thresh);     % edge is 0, non-edge is 1 for regularization, note wG.^2 = wG
