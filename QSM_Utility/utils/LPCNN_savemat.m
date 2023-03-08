function status = LPCNN_savemat(freqMap_file_mat)

[data_path, base_name, ~] = fileparts(freqMap_file_mat);
freqMap_file = fullfile(data_path, base_name);
mask_file = fullfile(data_path, [base_name, '_mask']);
dipole_file = fullfile(data_path, [base_name, '_kernel.mat']);

if ~exist(dipole_file, 'file')
    S = load(freqMap_file_mat);

    % save freqMap to NII
    saveNII(S.freqMap, freqMap_file, S.Params, 1);  

    % save mask to NII
    saveNII(squeeze(S.maskErode.*1), mask_file, S.Params, 1);  

    % and save dipole kernel to .mat
    if prod(S.Params.sizeVol) > 1e7 % in case of ultra high resolution
        datatype = 'single';
    else
        datatype = 'double';
    end
    D = conv_kernel_rot_c0(S.Params, S.Params.TAng, datatype);
    D = ifftshift(D);
    save(dipole_file, 'D')
end

status = 0;
