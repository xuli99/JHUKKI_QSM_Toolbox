function GRE_preparation_2nifti(dcm, output_dir, cleanup)
% function GRE_preparation_2nifti(dcm, output_dir, cleanup)
% 
% preparation for multi-echo GRE data, convert from DICOM file/folder to NIFTI using dcm2niix  
% then combine multi-echo GRE data and read Params from json file and save to .mat 
% 
% input_dir: input DICOM file or folder
% ouput_dir: optional, default to nifti
% cleanup  : cleanup the original NIFTI and json files from dcm2niix
% 
%% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu

if ismac
    dcm2niix_path = '/Users/lixu/opt/anaconda3/bin/';
else
    dcm2niix_path = '/usr/local/bin/';
end

if nargin < 2
    output_dir = [];
    cleanup = 1;
elseif nargin < 3
    cleanup = 1;
end

% if multiple files
if ~iscell(dcm)
    dcm = {dcm};
end

%% dcm2niix conversion & echo combination
for dcm_ii = 1:length(dcm)

    dcm_name = dcm{dcm_ii};

    if isfile(dcm_name)
        dcm_pdir = dir(dcm_name).folder;                    % for enhanced DICOM
    elseif isfolder(dcm_name)
        if dcm_name(end) == filesep, dcm_name = dcm_name(1:end-1); end
        dcm_pdir_parts = strsplit(dcm_name, filesep);       % for classical DICOM dir
        dcm_pdir = strjoin(dcm_pdir_parts(1:end-1), filesep);
    end

    disp('Transfer: ')
    disp(dcm_name);

    if isempty(output_dir)
        nii_output = fullfile(dcm_pdir, 'nifti'); % default inter output
        if ~exist(nii_output, "dir")
            mkdir(nii_output)
        end
        output_dir_internal = nii_output;
    else
        output_dir_internal = output_dir;
    end

    cmd = [dcm2niix_path, 'dcm2niix -z y -v y -w 1 -f %p_%s -p y -o ', output_dir_internal, ' ', dcm_name];
    system(cmd);

    % combine nifti file from each echo    
    json_list = dir(fullfile(output_dir_internal, '*_e1.json'));    % magnitude echo 1
    filename_prefix = extractBefore(json_list(1).name, 'e1.json');

    json_list_ph = dir(fullfile(output_dir_internal, '*_e*_ph.json')); 
    num_echo = length(json_list_ph);

    % magnitude & phase
    img_mag     = [];
    img_phase   = [];

    for kecho = 1:num_echo
        curr_echo_mag_filename      = strcat(filename_prefix, 'e', num2str(kecho),'.nii.gz');
        curr_echo_phase_filename    = strcat(filename_prefix, 'e', num2str(kecho),'_ph.nii.gz');
        
        img_mag     = cat(4,img_mag, load_nii_img_only(fullfile(output_dir_internal,curr_echo_mag_filename)));
        img_phase   = cat(4,img_phase, load_nii_img_only(fullfile(output_dir_internal,curr_echo_phase_filename)));
    end

    output_mag_filename     = strcat(filename_prefix, 'GRE_mag.nii.gz');
    output_phase_filename   = strcat(filename_prefix, 'GRE_phase.nii.gz');

    save_nii_img_only(fullfile(output_dir_internal,curr_echo_mag_filename),fullfile(output_dir_internal,output_mag_filename),img_mag);
    save_nii_img_only(fullfile(output_dir_internal,curr_echo_phase_filename),fullfile(output_dir_internal,output_phase_filename),img_phase);

    disp(['GRE mag & phase data saved at ', output_dir_internal, '.'])

    % Extract Params from json files
    nii_phase = load_untouch_nii(fullfile(output_dir_internal,output_phase_filename));

    % Params only for NIFTI output, thus should be single volume
    Params.nifti_hdr = nii_phase.hdr;   % nifti head for output, with multi-layer
    Params.nifti_hdr.dime.dim(5) = 1;   % should be echo combined
    Params.nifti_hdr.dime.pixdim(5) = 0;
    Params.nifti_hdr.dime.dim(1) = 3;

    % read TE from json file
    json_list = dir(fullfile(output_dir_internal, [filename_prefix, '*_ph.json']));
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

    % save header .mat file
    output_header_filename     = strcat(filename_prefix, 'header.mat');
    save(fullfile(output_dir_internal, output_header_filename), "Params", '-v7.3');

    disp('GRE header .mat file saved.')

    % clean up if selected
    if cleanup
        delete(fullfile(output_dir_internal, [filename_prefix, 'e*.json']));
        delete(fullfile(output_dir_internal, [filename_prefix, 'e*.nii.gz']));
    end

end