function mask_unrelyPhase = create_mask_unrelyPhase(GREPhase, GradMag_thresh, maskErode)
%% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
%
% This function creates unreliable phase mask for QSM
% updated 2019-06-25, X.L.
% updated 2021-01-01, X.L.

if nargin < 2
    GradMag_thresh = pi/2;  % practically determined
end

Nkernel = [5, 5, 5];    % may need to change according to image resolution

TV = TVOP;

num_echoes = size(GREPhase, 4);
mask_unrelyPhase = zeros(size(GREPhase(:,:,:,1)));  % 3D mask
Adapt_thresh = GradMag_thresh;

for echoii = 1:num_echoes
    phaseGrad = TV*GREPhase(:,:,:,echoii);    
    
    if exist('medfilt3', 'file')
       phaseGradMag = medfilt3(sqrt(sum(phaseGrad.^2, 4)./3), 'zeros', Nkernel); % matlab internal one
    else
       disp('matlab version before R2016, use mymedfilt3.')
       phaseGradMag = mymedfilt3(sqrt(sum(phaseGrad.^2, 4)./3), Nkernel);        
    end
    
    if ~isempty(maskErode)
        Adapt_thresh = prctile(phaseGradMag(maskErode>0), 99);
    end    
    mask_unrelyPhase = mask_unrelyPhase | (phaseGradMag > min([GradMag_thresh, Adapt_thresh]));
end

