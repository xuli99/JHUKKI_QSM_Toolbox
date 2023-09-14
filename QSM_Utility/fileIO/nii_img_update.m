function nii = nii_img_update(nii, Params)
% function nii = nii_img_update(nii, Params)
%
% Author: Xu Li
% Affiliation: Radiology @ JHU - Kirby Center
% Contact via xuli@mri.jhu.edu
% 
% flip nii.img if needed to fit the nifti_hdr
% check Params.nifti_flp_sli and Params.nifti_flp


% if Params has nifti_flp_sli, flip slices if needed
if isfield(Params, 'nifti_flp_sli')
    if Params.nifti_flp_sli     % flip slices to fit nifti_hdr
        img_data = nii.img; img_data = flip(img_data, 3); nii.img = img_data;
    end
end

if isfield(Params, 'nifti_flp')
    img_data = nii.img;
    for k = 1:3, if Params.nifti_flp(k), img_data = flip(img_data, k); end; end
    nii.img = img_data;
end

