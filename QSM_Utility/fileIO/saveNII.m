function []=saveNII(saveVar, fileName, Params, permuteflag, ext, prec)
    % Wrapper function for NIFTI save operation
    % updated for v2.5, X.L., 2014-09-08
    % updated for v2.8, X.L., 2016-12-08
    % updated to include .img/.hdr format, X.L., 2017-05-31
    % 2019-02-18, to include dynamics data
    % 2019-06-05, included .nii.gz and prec, X.L.
    % 2020-04-07, updated sliceOri option for saving saggital images
    % 2021-01-01, added SliceOriSave
    
    if nargin < 4
        permuteflag = 1;        % default setting, compatible with old version
        ext = '.nii.gz';
        prec = 16;
    elseif nargin < 5
        ext = '.nii.gz';           % default using nii
        prec = 16;
    elseif nargin < 6
        prec = 16;
    end
    
    if permuteflag == 1
        saveVar = permute(saveVar, [2, 1, 3:length(size(saveVar))]);
        Params.voxSize = Params.voxSize([2, 1, 3]);
    end
    
    if isfield(Params, 'SliceOriSave')
        Params.sliceOri = Params.SliceOriSave;
    end
    
    if ~isfield(Params, 'sliceOri')
        Params.sliceOri = 1;    % 1:axial, 2:saggital, 3:coronal, default to axial
    end
    
    if strcmpi(ext, '.nii') || strcmpi(ext, '.nii.gz')
        nii = make_nii(saveVar, Params.voxSize, [], prec);       % single, float 32
        % save as nifti
        switch Params.sliceOri
            case 1
                nii.hdr.hist.qform_code = 1;
                nii.hdr.hist.sform_code = 0;
                nii.hdr.hist.quatern_b = 0;
                nii.hdr.hist.quatern_c = 0;
                nii.hdr.hist.quatern_d = 1;        
            case 2
                nii.hdr.hist.qform_code = 1;
                nii.hdr.hist.sform_code = 0;
                nii.hdr.hist.quatern_b = -0.5;
                nii.hdr.hist.quatern_c = 0.5;
                nii.hdr.hist.quatern_d = -0.5;  
            case 3
                nii.hdr.hist.qform_code = 1;
                nii.hdr.hist.sform_code = 0;
                nii.hdr.hist.quatern_b =  0;
                nii.hdr.hist.quatern_c = 0.7071;
                nii.hdr.hist.quatern_d = -0.7071;  
                % nii.hdr.hist.qoffest_x = 0;   % may need to change
            otherwise
                disp('unknow slice orientation')
        end
        save_nii(nii, [fileName, ext]);
    else
        % save as analyze
        nii = make_nii(flip(saveVar, 2), Params.voxSize, [], prec);
        save_nii(nii, [fileName, '']);
    end
    
    clear nii
end

