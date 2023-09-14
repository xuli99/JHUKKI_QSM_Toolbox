function []=saveNII(saveVar, fileName, Params, permuteflag, ext, prec)
    % Wrapper function for NIFTI save operation
    % updated for v2.5, X.L., 2014-09-08
    % updated for v2.8, X.L., 2016-12-08
    % updated to include .img/.hdr format, X.L., 2017-05-31
    % 2019-02-18, to include dynamics data
    % 2019-06-05, included .nii.gz and prec, X.L.
    % 2020-04-07, updated sliceOri option for saving saggital images
    % 2021-01-01, added SliceOriSave
    % 2022-03-22, added option to use loaded NIFTI hdr from DICOM (dicm2nii)
    % 2023-05-08, added option for NIFTI hdr from nifti
    % 2023-09-11, bug fix for saving NIFTI RAS with flipped slice order
    
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

        % if Params has saved nifti_hdr, which should be RAS
        if isfield(Params, 'nifti_hdr')
            % check dim
            % nii = nii_tool('init', saveVar);
            % if nii.hdr.dim == Params.nifti_hdr.dim
            if isfield(Params.nifti_hdr, 'dime')
                dim = Params.nifti_hdr.dime.dim;
            else
                dim = Params.nifti_hdr.dim;
            end

            if nii.hdr.dime.dim == dim
                % if Params has nifti_flp_sli, flip slices if needed
                if isfield(Params, 'nifti_flp_sli')
                    if Params.nifti_flp_sli     % flip slices to fit nifti_hdr
                        img_data = nii.img;
                        img_data = flip(img_data, 3);
                        nii.img = img_data;
                    end
                end

                % flip if needed.
                if isfield(Params, 'nifti_flp')
                    img_data = nii.img;
                    for k = 1:3
                        if Params.nifti_flp(k), img_data = flip(img_data, k); end
                    end
                    nii.img = img_data;
                end
                nii.hdr = hdr_update(nii.hdr, Params.nifti_hdr);

                % nii_tool('save', nii, [fileName, ext])
            else
                error('data dimension does not match saved NIFTI header.')
            end

        else
            % save as nifti without PatientPosition
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
        end

        save_nii(nii, [fileName, ext]);

    else
        % save as analyze
        nii = make_nii(flip(saveVar, 2), Params.voxSize, [], prec);
        save_nii(nii, [fileName, '']);
    end
    
    clear nii

    function hdr = hdr_update(hdr, new_hdr)
        if isfield(new_hdr, 'dime')
            new_hdr.pixdim = new_hdr.dime.pixdim;
            new_hdr.vox_offset = new_hdr.dime.vox_offset;

            new_hdr.qform_code = new_hdr.hist.qform_code;
            new_hdr.sform_code = new_hdr.hist.sform_code;
            new_hdr.quatern_b = new_hdr.hist.quatern_b;
            new_hdr.quatern_c = new_hdr.hist.quatern_c;
            new_hdr.quatern_d = new_hdr.hist.quatern_d;
            new_hdr.qoffset_x = new_hdr.hist.qoffset_x;
            new_hdr.qoffset_y = new_hdr.hist.qoffset_y;
            new_hdr.qoffset_z = new_hdr.hist.qoffset_z;
            new_hdr.srow_x = new_hdr.hist.srow_x;
            new_hdr.srow_y = new_hdr.hist.srow_y;
            new_hdr.srow_z = new_hdr.hist.srow_z;
            new_hdr.magic = new_hdr.hist.magic;

        end

        if isfield(hdr, 'dime')
            % multi-layer way
            hdr.dime.pixdim = new_hdr.pixdim; hdr.dime.vox_offset = new_hdr.vox_offset;
            hdr.hist.qform_code = new_hdr.qform_code;
            hdr.hist.sform_code = new_hdr.sform_code;
            hdr.hist.quatern_b = new_hdr.quatern_b; hdr.hist.quatern_c = new_hdr.quatern_c; hdr.hist.quatern_d = new_hdr.quatern_d;
            hdr.hist.qoffset_x = new_hdr.qoffset_x; hdr.hist.qoffset_y = new_hdr.qoffset_y; hdr.hist.qoffset_z = new_hdr.qoffset_z;
            hdr.hist.srow_x = new_hdr.srow_x; hdr.hist.srow_y = new_hdr.srow_y; hdr.hist.srow_z = new_hdr.srow_z;
            hdr.hist.magic = new_hdr.magic;            
        else
            % single-layer way
            hdr.pixdim = new_hdr.pixdim; hdr.vox_offset = new_hdr.vox_offset;
            hdr.qform_code = new_hdr.qform_code;
            hdr.sform_code = new_hdr.sform_code;
            hdr.quatern_b = new_hdr.quatern_b; hdr.quatern_c = new_hdr.quatern_c; hdr.quatern_d = new_hdr.quatern_d;
            hdr.qoffset_x = new_hdr.qoffset_x; hdr.qoffset_y = new_hdr.qoffset_y; hdr.qoffset_z = new_hdr.qoffset_z;
            hdr.srow_x = new_hdr.srow_x; hdr.srow_y = new_hdr.srow_y; hdr.srow_z = new_hdr.srow_z;
            hdr.magic = new_hdr.magic;
        end
    end

end

