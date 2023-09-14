%% function save_nii_img_only(headerfilename,savefilename,images)
%
% Description: Wrapper script for saving with target image header
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 24 January 2017
% Date modified: 28 June 2020 (v0.8.0)
%
% cleaned up X.L.
function save_nii_img_only(headerfilename,savefilename,img,datatype)

    if nargin < 4
        datatype = 16;
    end
    
    % load target image header
    nii = load_untouch_nii(headerfilename);
    
    % fill nii.img with img and update dim with multi-echo
    nii.img = single(img);
    nii.hdr.dime.datatype = datatype;
    nii.hdr.dime.dim(5) = size(img,4);
    nii.hdr.dime.dim(1) = ndims(img);

    % assume the input image contains the true values
    nii.hdr.dime.scl_inter = 0;
    nii.hdr.dime.scl_slope = 1;

    save_untouch_nii(nii,savefilename);
end