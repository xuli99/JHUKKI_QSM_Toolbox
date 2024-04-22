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

% 2024-03-21, X.L., for SIEMENS data with diff series number for mag vs. phase
%                   dcm2niix may have phase scaling problem, need to check

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
        dcm_name_dir = dir(dcm_name);                    % for enhanced DICOM
        dcm_pdir = dcm_name_dir.folder;
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

    % for debuging add -v y
    cmd = [dcm2niix_path, 'dcm2niix -z y -v n -w 1 -f %p_%s -p y -o ', output_dir_internal, ' ', dcm_name];
    system(cmd);

    % combine nifti files, call read_nifti_combine.m
    parrec_flag = 0;
    if isfile(dcm_name) && contains(dcm_name(end-3:end), ["par","rec"], 'IgnoreCase', true)
        disp('par/rec files detected.')
        parrec_flag = 1;
    end
    
    read_nifti_combine(output_dir_internal, cleanup, parrec_flag);

end