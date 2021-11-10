function [GREPhase, GREMag,Params] = read_DICOM(DICOMdir, Params, verbose)
% function [GREPhase, GREMag,Params] = read_DICOM(DICOMdir, Params, verbose)
%
%% Author: Xu Li
%% EDITED BY JIRI VAN BERGEN - 2014
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
%
% read in DICOM format SWI images (Mag or Phase)
%
% input:
%        DICOMdir: DICOM directory
%        Params: The basic params from Constants.m
% output:
%        data: image data
%        Params.voxSize
%        Params.sizeVol
%        Params.B0: main field strength
%        Params.TE echo time (array if multi-echo, TE1 otherwise)
%        Params.TAng
%        Params.TR
%
% created on 2013-11-15, modified from DICOM reader from MEDI toolbox
% 
% Updated on 2014-10-10 for old Philips DICOM format, X.L.
% Updated on 2014-10-22 for Philips DICOM format, X.L.
% Updated on 2015-04-07 for SIEMENS uncombined data with ASCII code, X.L.
% Updated 2017-04-19, X.L. no round TEs
% Updated 2019-02-15, changed to original 5D/6D data format
% Updated 2019-04-22, added in support for GE data, X.L.
% Updated 2019-05-20, added in support for TOSHIBA/Canon data, X.L.
% Updated 2019-11-08, fixed a bug for info.ImageType
% Updated 2019-12-18, updated according to new DICOM tags (no info.EchoNumbers)
% Updated 2020-10-21, fixed a bug for info.EchoNumbers
% Updated 2020-10-31, for EchoNumber(s), space in filenames and
%                       SliceOriSave
% Updated 2021-04-03, used fullfile
% Updated 2021-04-26, fixed in case no tag on EchoNumber(s)
% Updated 2021-06-26, add update for cluster version
% Updated 2021-10-28, for reverse slice stack condition, making TAng always R.H.S.

if nargin < 1
    % uigetdir get the DICOMdir
    DICOMdir = uigetdir(pwd, 'Select DICOM directory');
    if DICOMdir == 0
        error('DICOM folder not selected.')
    end
    verbose = 1;
elseif nargin < 3
    verbose = 1;
end

% Fancy waitbar
textWaitbar = 'Reading in DICOM-files';
if ~isfield(Params, 'cluster')
    multiWaitbar(textWaitbar, 0);
else
    disp(textWaitbar);
end

% Basics
cd(DICOMdir);
[Params.PathName, Params.FileBaseName,~] = fileparts(pwd);
Params.PathName = fullfile(Params.PathName, Params.FileBaseName);
Params.FileBaseName = strrep(Params.FileBaseName, ' ', '_');

% Get files
filelist = dir(DICOMdir);

ii=1;
while ii<=length(filelist)
    if filelist(ii).isdir==1 || ~isDICOM(filelist(ii).name)
        filelist = filelist([1:ii-1 ii+1:end]);   % skip folders
    else
        ii=ii+1;
    end
end

filenameBase = fullfile(DICOMdir, filelist(1).name);

info = dicominfo(filenameBase);
taglist = fieldnames(info);
if contains(taglist, 'EchoNumber')
    tagEchoNumber = taglist{contains(taglist, 'EchoNumber')};
else
    tagEchoNumber = [];
end

try
    % check manufacture
    if contains(info.Manufacturer, 'SIEMENS', 'IgnoreCase', true)
        Manufacturer = 1;
    elseif contains(info.Manufacturer, 'Philips', 'IgnoreCase', true)
        Manufacturer = 2;
    elseif contains(info.Manufacturer, 'GE', 'IgnoreCase', true)
        Manufacturer = 3;
    elseif contains(info.Manufacturer, 'TOSHIBA', 'IgnoreCase', true)
        Manufacturer = 4;
    else
        error('unknow Manufacturer ... ')
    end
catch
    % Other cases
    if contains(info.ManufacturerModelName, 'Skyra', 'IgnoreCase', true)
        Manufacturer = 1;
    else
        error('Unknown Manufacturer ... ')
    end
end

% basic information
Params.sizeVol(1) = single(info.Width);
Params.sizeVol(2) = single(info.Height);

Params.voxSize(1) = single(info.PixelSpacing(1));       % mm
Params.voxSize(2) = single(info.PixelSpacing(2));
Params.voxSize(3) = single(info.SliceThickness);        % mm

% CF = info.ImagingFrequency*1e6;                    % central frequency, Hz
Params.B0 = info.MagneticFieldStrength;                 % in T
Params.TR = info.RepetitionTime*1e-3;                   % TR in sec

% Angulation matrix
Affine2D = reshape(info.ImageOrientationPatient,[3 2]);
Params.TAng = [Affine2D, cross(Affine2D(:,1), Affine2D(:,2))];

% Dynamics imaging
if isfield(info, 'NumberOfTemporalPositions')        % Philips/TOSHIBA
    Params.DynamicNum = info.NumberOfTemporalPositions;
else
    Params.DynamicNum = 1;
end

% coil element combination
if ~isfield(Params, 'coilNum')                              % Params not set yet
    Params.coilNum = 1;                                     % default, single coil or combined
    if isfield(info, 'Private_0051_100f') && ~isfield(Params, 'coilName')         % SIEMENS    
        Params.coilName(1) = cellstr(char((info.Private_0051_100f(:))'));        % coil name, could be ASCII code
    end
end

% Number of Echoes
minSlice = 1e10;
maxSlice = -1e10;

if isfield(info, tagEchoNumber)
    NumEcho = info.(tagEchoNumber);          % first check EchoNumber
elseif isfield(info, 'EchoTime')    
    NumEcho = 1;                             % at least one
elseif isfield(info, 'EchoTrainLength')
    NumEcho = info.EchoTrainLength;
else
    NumEcho = 0;
end

TE = [];
for ii = 1:length(filelist)
    info = dicominfo(fullfile(DICOMdir, filelist(ii).name));
    if info.SliceLocation<minSlice                  % find lowest slice
        minSlice = info.SliceLocation;              % in mm
        minLoc = info.ImagePositionPatient;
    end
    if info.SliceLocation>maxSlice                  % find highest slice
        maxSlice = info.SliceLocation;
        maxLoc = info.ImagePositionPatient;
    end
    
    if isfield(info, tagEchoNumber) && (info.(tagEchoNumber)>NumEcho) % in case multi-echo
        NumEcho = info.(tagEchoNumber);
    end
    
    if isfield(info, 'EchoTime')
        if isempty(find(TE==info.EchoTime, 1))
            TE(end+1) = info.EchoTime;
            TE = sort(TE);
        end
    end
    
    if isfield(info, 'Private_0051_100f')              % get all the coil names
        if isempty(find(strcmpi(Params.coilName, cellstr(char((info.Private_0051_100f(:))'))) > 0, 1, 'first'))  % new coil name
            Params.coilName(length(Params.coilName)+1) = cellstr(char((info.Private_0051_100f(:))'));
        end   
        Params.coilNum = length(Params.coilName);
    end
    
    if verbose        
        disp(['searching DICOM ', num2str(ii), ' ...']);
    end
end

if length(TE) ~= NumEcho && ~isempty(tagEchoNumber)
    error('length of TEs does not equal NumEcho, need to check ...')
else
    % NumEcho was guessed
    NumEcho = length(TE);
end

% Number of Slices
if isfield(info, 'SpacingBetweenSlices')
    Params.sizeVol(3) = round(norm(maxLoc - minLoc)/info.SpacingBetweenSlices) + 1;
else
    Params.sizeVol(3) = round(norm(maxLoc - minLoc)/Params.voxSize(3)) + 1;
end

slicenorm = (maxLoc - minLoc)/(Params.sizeVol(3)-1);
if dot(Params.TAng(:,3), slicenorm) < 0  % needed to reverse slice stack
    reverse_flag = 1;
end

%% read in imaging data
if Manufacturer ~= 3   % for SIEMENS or PHILIPS or TOSHIBA, can read in Mag & Phase directly
    GREPhase = single(zeros([Params.sizeVol NumEcho, Params.DynamicNum, Params.coilNum]));
    GREMag = single(zeros([Params.sizeVol NumEcho, Params.DynamicNum, Params.coilNum]));
else
    % For GE data, read in Real/Imaginary then convert to Mag/Phase
    GREReal = single(zeros([Params.sizeVol NumEcho, Params.DynamicNum, Params.coilNum]));
    GREImag = single(zeros([Params.sizeVol NumEcho, Params.DynamicNum, Params.coilNum]));
    
    GREPhase = single(zeros([Params.sizeVol NumEcho, Params.DynamicNum, Params.coilNum]));
    GREMag = single(zeros([Params.sizeVol NumEcho, Params.DynamicNum, Params.coilNum]));
end

for i = 1:length(filelist)
    info = dicominfo(fullfile(DICOMdir, filelist(i).name));

    if isfield(info, 'SpacingBetweenSlices') && Manufacturer ~= 3   % slice
        slice = int32(round(norm(info.ImagePositionPatient-minLoc)/info.SpacingBetweenSlices) +1);      % slice number
    else
        slice = int32(round(norm(info.ImagePositionPatient-minLoc)/Params.voxSize(3)) +1);      % slice number
    end
    
    if isfield(info, 'TemporalPositionIdentifier')          % dynamics
        dynamic = info.TemporalPositionIdentifier;
    else
        dynamic = 1;
    end
    
    if isfield(info, 'Private_0051_100f')                   % SIEMENS coil
        coil = find(strcmpi(Params.coilName, cellstr(char((info.Private_0051_100f(:))')))>0, 1, 'first');        % coil name
    else
        coil = 1;
    end        

    if ~isfield(info, 'EchoNumber')
        info.EchoNumber = find(TE==info.EchoTime);
    end
    
    if ~isfield(info, 'LargestImagePixelValue') && isfield(info, 'BitsStored')
        info.LargestImagePixelValue = 2^(info.BitsStored - 1);
    end
    
    % read in image data
    if Manufacturer < 3         % for SIEMENS and PHILIPS data
        if strcmpi(info.ImageType(18), 'P') || contains(info.ImageType, '\P\')             % phase
            switch Manufacturer
                case 1
                    GREPhase(:,:,slice,info.EchoNumber, dynamic, coil)  = (single(dicomread(fullfile(DICOMdir, filelist(i).name))')*info.RescaleSlope+info.RescaleIntercept)/single(info.LargestImagePixelValue)*pi;    %phase
                case 2
                    if isfield(info, 'RealWorldValueMappingSequence')
                        GREPhase(:,:,slice,info.EchoNumber, dynamic, coil)  = 1e-3*(single(dicomread(fullfile(DICOMdir, filelist(i).name))')*info.RealWorldValueMappingSequence.Item_1.RealWorldValueSlope+info.RealWorldValueMappingSequence.Item_1.RealWorldValueIntercept); %phase
                    elseif isfield(info, 'RescaleSlope') && matches(info.RescaleType, ["milliradials", "mrad"])
                        GREPhase(:,:,slice,info.EchoNumber, dynamic, coil)  = 1e-3*(single(dicomread(fullfile(DICOMdir, filelist(i).name))')*info.RescaleSlope+info.RescaleIntercept); %phase
                    else
                        error('unknown scaling ...')
                    end
            end

        elseif strcmpi(info.ImageType(18), 'M') || contains(info.ImageType, '\M\')         % magnitude
            GREMag(:,:,slice,info.EchoNumber, dynamic, coil)  = single(dicomread(fullfile(DICOMdir, filelist(i).name))');     %magnitude
        end
    elseif Manufacturer == 3
        % for GE data
        if mod(info.InstanceNumber,4)==0
            GREImag(:,:,slice,info.EchoNumber, dynamic, coil)  = single(dicomread(fullfile(DICOMdir, filelist(i).name))');   %imaginary
        elseif mod(info.InstanceNumber,4)==3
            GREReal(:,:,slice,info.EchoNumber, dynamic, coil)  = single(dicomread(fullfile(DICOMdir, filelist(i).name))');   %real
        end     
    elseif Manufacturer == 4
        % for TOSHIBA
        if info.WindowCenter == 0   % phase
            GREPhase(:,:,slice,info.EchoNumber, dynamic, coil)  = (single(dicomread(fullfile(DICOMdir, filelist(i).name))'))/single(info.WindowWidth)*2*pi;
        else
            GREMag(:,:,slice,info.EchoNumber, dynamic, coil)  = single(dicomread(fullfile(DICOMdir, filelist(i).name))');     %magnitude
        end
    end
    
    % Only show waitbar every 10 iterations
    if(mod(i, 50) == 0)
        if ~isfield(Params, 'cluster')
            multiWaitbar(textWaitbar, (i/length(filelist)));
        else
            disp([num2str(100*(i/length(filelist))), '% Done.']);
        end
    end
end

% Close waitbar
if ~isfield(Params, 'cluster')
    multiWaitbar( 'CloseAll' );
else
    disp('Done.')
end

if Manufacturer == 3    % for GE data with phase shift artifacts
    temp = ifft(fftshift(fft(GREReal + 1i.*GREImag, [], 3), 3), [], 3);    
    temp = conj(temp);  % for GE data, use the conj (-phase)
    GREMag = abs(temp);
    GREPhase = angle(temp);
    clear temp
end

if Manufacturer == 4        % TOSHIBA data is flipped?
    GREMag = flip(GREMag, 3);
    GREPhase = flip(GREPhase, 3);
end

% reverse slice stack
if reverse_flag == 1
    GREMag = GREMag(:,:,end:-1:1,:,:,:);
    GREPhase = GREPhase(:,:,end:-1:1,:,:,:);
end

% Post-process
Params.nEchoes  = NumEcho;
Params.voxSize  = (Params.voxSize);
Params.sizeVol  = double(Params.sizeVol);
Params.fov      = (round(Params.voxSize.*Params.sizeVol));
Params.TEs      = TE./1000;
Params.nDynamics= Params.DynamicNum;

[vcos, LPS] = max(abs(Params.TAng(:,3)));
if vcos > 0.7
    switch LPS
        case 1
            Params.SliceOriSave = 2;    % Sagittal
        case 2
            Params.SliceOriSave = 3;    % Coronal
        case 3
            Params.SliceOriSave = 1;    % Axial
        otherwise
            error('unexpected orientation.')
    end
end
% data in 6D array, ncol, nrow, nslice, necho, ndyanmics, ncoil
