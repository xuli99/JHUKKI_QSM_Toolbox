function read_nifti_combine(nifti_dirs, cleanup, parrecflag)
% function read_nifti_combine(nifti_dirs, cleanup, parrecflag)
% 
% combine multi-echo GRE data (NIFTI from dcm2niix) and read Params from json file and save to header.mat 
% 
% nifti_dirs: folders with nifti output from dcm2niix
% cleanup   : cleanup the original NIFTI and json files from dcm2niix
% 
%% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu

% 2024-03-21, X.L., for SIEMENS data with diff series number for mag vs. phase

if nargin < 2
    cleanup = 0;        % default no cleanup
    parrecflag = 0;
elseif nargin < 3
    parrecflag = 0;
end

% if multiple folders
if ~iscell(nifti_dirs)
    nifti_dirs = {nifti_dirs};
end

%% echo combination
for nifti_ii = 1:length(nifti_dirs)

    nifti_dir = nifti_dirs{nifti_ii};

    disp(nifti_dir);

    % combine nifti file from each echo
    json_test = dir(fullfile(nifti_dir, '*_ph.json'));

    if length(json_test) > 1
        % if with multi-echoes
        json_list = dir(fullfile(nifti_dir, '*_e1.json'));            % magnitude echo 1
        filename_prefix_mag = extractBefore(json_list(1).name, '_e1.json');     % prefix without "_" now
    
        json_list_ph = dir(fullfile(nifti_dir, '*_e*_ph.json')); 
        num_echo = length(json_list_ph);
        filename_prefix_phase = extractBefore(json_list_ph(1).name, '_e1_ph.json');     % prefix without "_" now

    elseif length(json_test) == 1
        % if with single-echo
        num_echo = 1;
        json_list_ph = json_test;
        filename_prefix_phase = extractBefore(json_list_ph(1).name, '_ph.json');  % prefix without "_"

        json_list = dir(fullfile(nifti_dir, '*.json'));
        for json_ii = 1:lenght(json_list)
            if ~contains(json_list(json_ii).name, json_list_ph(1).name)
                filename_prefix_mag = extractBefore(json_list(json_ii).name, '.json');
            end
        end

    else
        error('There is no phase data.')
    end

    % magnitude & phase
    img_mag     = [];
    img_phase   = [];

    for kecho = 1:num_echo
        if num_echo > 1
            curr_echo_mag_filename      = strcat(filename_prefix_mag, '_e', num2str(kecho),'.nii.gz');
            curr_echo_phase_filename    = strcat(filename_prefix_phase, '_e', num2str(kecho),'_ph.nii.gz');
        else
            curr_echo_mag_filename      = strcat(filename_prefix_mag, '.nii.gz'); % prefix without "_"
            curr_echo_phase_filename    = strcat(filename_prefix_phase, '_ph.nii.gz');
        end

        img_mag     = cat(4,img_mag, load_nii_img_only(fullfile(nifti_dir,curr_echo_mag_filename)));
        img_phase   = cat(4,img_phase, load_nii_img_only(fullfile(nifti_dir,curr_echo_phase_filename)));
        
    end

    % check phase scaling 
    phase_range = max(img_phase(:)) - min(img_phase(:));
    if phase_range > 2*pi + 100*eps
        disp('correct phase scaling ...')
        img_phase = img_phase./phase_range*(2*pi);
    end

    % find common prefix & save
    filename_prefix = intersect(filename_prefix_mag, filename_prefix_phase, 'stable');

    output_mag_filename     = strcat(filename_prefix, '_GRE_mag.nii.gz');
    output_phase_filename   = strcat(filename_prefix, '_GRE_phase.nii.gz');

    save_nii_img_only(fullfile(nifti_dir,curr_echo_mag_filename),fullfile(nifti_dir,output_mag_filename),img_mag);
    save_nii_img_only(fullfile(nifti_dir,curr_echo_phase_filename),fullfile(nifti_dir,output_phase_filename),img_phase);

    disp(['GRE mag & phase data saved at ', nifti_dir, '.'])

    % Extract Params from json files
    nii_phase = load_untouch_nii(fullfile(nifti_dir,output_phase_filename));

    % Params only for NIFTI output, thus should be single volume
    Params.nifti_hdr = nii_phase.hdr;   % nifti head for output, with multi-layer
    Params.nifti_hdr.dime.dim(5) = 1;   % should be echo combined
    Params.nifti_hdr.dime.pixdim(5) = 0;
    Params.nifti_hdr.dime.dim(1) = 3;

    % read TE from json file
    json_list = dir(fullfile(nifti_dir, [filename_prefix_phase, '*_ph.json']));
    json_filenames = cell(length(json_list), 1);
    for json_ii = 1:length(json_list)
        json_filenames{json_ii} = fullfile(json_list(json_ii).folder, json_list(json_ii).name);
    end
    Params.TEs = readTE_dcm2niix_JSON(json_filenames);
    Params.nEchoes = length(Params.TEs);

    % read other params from json file, B0, TR etc.
    Params = readParams_dcm2niix_JSON(json_filenames{1}, Params);

    Params.sizeVol = nii_phase.hdr.dime.dim(2:4);  % 2-4, pixdim, 5:echoes, 6:dynamics?
    Params.voxSize = nii_phase.hdr.dime.pixdim(2:4);
    Params.fov = Params.sizeVol.*Params.voxSize;
    Params.nDynamics = nii_phase.hdr.dime.dim(6);  % need check 

    [Params.b0dir, Params.TAng] = get_B0_dir_from_nifti(nii_phase);

    % par/rec to NIFTI test, convert to LAS, b0 needs correction
    if parrecflag == 1
       disp('usning NIFTI file converted from par/rec, check with caution...')
       negz = diag([1, 1, -1]);   %
       Params.b0dir = negz*Params.b0dir; Params.TAng = Params.TAng*negz;
    end

    % save header .mat file
    output_header_filename     = strcat(filename_prefix, '_header.mat');
    save(fullfile(nifti_dir, output_header_filename), "Params", '-v7.3');

    disp('GRE header .mat file saved.')

    % clean up if selected
    if cleanup
        if num_echo > 1
            % if multi-echo
            delete(fullfile(nifti_dir, [filename_prefix_mag, '_e*.json']));
            delete(fullfile(nifti_dir, [filename_prefix_mag, '_e*.nii.gz']));
            delete(fullfile(nifti_dir, [filename_prefix_phase, '_e*.json']));
            delete(fullfile(nifti_dir, [filename_prefix_phase, '_e*.nii.gz']));
        else
            % if single-echo
            delete(fullfile(nifti_dir, [filename_prefix_mag, '.json']));
            delete(fullfile(nifti_dir, [filename_prefix_mag, '.nii.gz']));
            delete(fullfile(nifti_dir, [filename_prefix_phase, '_ph.json']));
            delete(fullfile(nifti_dir, [filename_prefix_phase, '_ph.nii.gz']));
        end

    end

end