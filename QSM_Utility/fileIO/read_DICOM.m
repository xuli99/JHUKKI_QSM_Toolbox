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
% Updated 2021-10-28, for reverse slice stack condition, making TAng always R.H.S. (LPS)
% Updated 2021-11-18, bug fix
% Updated 2022-03-23, X.L. added option to use RAS NIFTI
% Updated 2022-11-28, X.L. deal with non-image DICOM files
% Updated 2023-09-11, X.L. bug fix for saving NIFTI RAS from DICOM LPS

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
Params.sizeVol(1) = single(info.Width);                 % columns
Params.sizeVol(2) = single(info.Height);                % rows

Params.voxSize(1) = single(info.PixelSpacing(2));       % mm, column spacing
Params.voxSize(2) = single(info.PixelSpacing(1));       % row spacing
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
valid_index = ones(length(filelist),1, "logical");
for ii = 1:length(filelist)
    info = dicominfo(fullfile(DICOMdir, filelist(ii).name));

    % in case there are non-image DICOM files
    if ~isfield(info, 'SliceLocation') || ~isfield(info, 'ImageType')
        valid_index(ii) = 0;
        continue;
    end

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

% in case there are non-image DICOM files, cleanup the filelist and update info
filelist = filelist(valid_index);
info = dicominfo(fullfile(DICOMdir, filelist(end).name));

if length(TE) ~= NumEcho && ~isempty(tagEchoNumber)
    error('length of TEs does not equal NumEcho, need to check ...')
else
    % NumEcho was guessed
    NumEcho = length(TE);
end

% Number of Slices
if isfield(info, 'SpacingBetweenSlices')
    Params.sizeVol(3) = round(norm(maxLoc - minLoc)/info.SpacingBetweenSlices) + 1;
    % in case of overcontiguous slices, not recommended
    if info.SpacingBetweenSlices < Params.voxSize(3)
        Params.voxSize(3) = info.SpacingBetweenSlices;
    end
else
    Params.sizeVol(3) = round(norm(maxLoc - minLoc)/Params.voxSize(3)) + 1;
end

slicenorm = (maxLoc - minLoc)/(Params.sizeVol(3)-1);

reverse_flag = 0;
if dot(Params.TAng(:,3), slicenorm) < 0  % needed to reverse slice stack for LPS
    reverse_flag = 1;
end

% % ----- NIFTI affine and hdr, tested for SIEMENS/PHILIPS data only
if Manufacturer <= 2
    testmat = abs(Params.TAng);
    [~, ixyz] = max(testmat);
    if ixyz(2) == ixyz(1), testmat(ixyz(2), 2) = 0; [~, ixyz(2)] = max(testmat(:,2)); end
    if any(ixyz(3) == ixyz(1:2)), ixyz(3) = setdiff(1:3, ixyz(1:2)); end
    
    pixdim = Params.voxSize;
    dim = Params.sizeVol;
    R = [Params.TAng * diag(pixdim) minLoc];
    R(:,3) = slicenorm;
    R(1:2,:) = -R(1:2, :);  % LPS to RAS
    
    flp = R(ixyz+[0 3 6])<0; % flip an axis if true
    d = det(R(:,1:3)) * prod(1-flp*2); % det after all 3 axis positive
    if d<0
       flp(1) = ~flp(1);    % right storage
    end
    rotM = diag([1-flp*2 1]); % 1 or -1 on diagnal
    rotM(1:3, 4) = (dim-1) .* flp; % 0 or dim-1
    R = R / rotM; % xform matrix after flip
    img_affine = eye(4); 
    img_affine(1:3,:) = R;

    img_data = zeros(dim,'single');
    nii = nii_tool('init', img_data);
    nii = nii_tool('update', nii, img_affine); % update only sform
    nii.hdr.sform_code = 1;
    nii.hdr.pixdim(:) = 1;
    nii.hdr.pixdim(2:4) = pixdim; 

    % from dicm2nii
    R0 = normc(R(:, 1:3)); iSL=3;
    sNorm = null(R0(:, setdiff(1:3, iSL))');
    if sign(sNorm(ixyz(iSL))) ~= sign(R(ixyz(iSL),iSL)), sNorm = -sNorm; end
    R0(:,iSL) = sNorm;

    % qform
    nii.hdr.qform_code = 1;
    nii.hdr.qoffset_x = R(1,4);
    nii.hdr.qoffset_y = R(2,4);
    nii.hdr.qoffset_z = R(3,4);
    [q, nii.hdr.pixdim(1)] = dcm2quat(R0); % 3x3 dir cos matrix to quaternion
    nii.hdr.quatern_b = q(2);
    nii.hdr.quatern_c = q(3);
    nii.hdr.quatern_d = q(4);

    Params.nifti_affine  = img_affine;
    Params.nifti_flp     = flp;
    Params.nifti_hdr     = nii.hdr;
    Params.nifti_flp_sli = reverse_flag;
end

%% read in imaging data
if Manufacturer ~= 3   % for SIEMENS or PHILIPS or TOSHIBA, can read in Mag & Phase directly
    GREPhase = single(zeros([Params.sizeVol NumEcho, Params.DynamicNum, Params.coilNum]));
    GREMag = single(zeros([Params.sizeVol NumEcho, Params.DynamicNum, Params.coilNum]));

%     % ****** In case there is no Phase, but with only Real/Imaginary for testing ******     
%     GREReal = single(zeros([Params.sizeVol NumEcho, Params.DynamicNum, Params.coilNum]));
%     GREImag = single(zeros([Params.sizeVol NumEcho, Params.DynamicNum, Params.coilNum]));
    
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
                    elseif isfield(info, 'RescaleSlope') && ( matches(info.RescaleType, ["milliradials", "mrad"]) || contains(info.RescaleType, 'normalized') )
                        GREPhase(:,:,slice,info.EchoNumber, dynamic, coil)  = 1e-3*(single(dicomread(fullfile(DICOMdir, filelist(i).name))')*info.RescaleSlope+info.RescaleIntercept); %phase
                    else
                        error('unknown scaling ...')
                    end
            end

        elseif strcmpi(info.ImageType(18), 'M') || contains(info.ImageType, '\M\')         % magnitude
            GREMag(:,:,slice,info.EchoNumber, dynamic, coil)  = single(dicomread(fullfile(DICOMdir, filelist(i).name))');     %magnitude

%        % ****** for testing  
%         elseif strcmpi(info.ImageType(18), 'R') || contains(info.ImageType, '\R\')         % real
%             GREReal(:,:,slice,info.EchoNumber, dynamic, coil)  = (single(dicomread(fullfile(DICOMdir, filelist(i).name))')*info.RescaleSlope+info.RescaleIntercept);    %real
% 
%         elseif strcmpi(info.ImageType(18), 'I') || contains(info.ImageType, '\I\')         % imaginary
%             GREImag(:,:,slice,info.EchoNumber, dynamic, coil)  = (single(dicomread(fullfile(DICOMdir, filelist(i).name))')*info.RescaleSlope+info.RescaleIntercept);    %imaginary
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

% % ****** for testing  
% temp = GREReal + 1i.*GREImag;
% GREPhase = angle(temp);

% reverse slice stack to make LPS, saved to Params.nifti_flp_sli, flip back
% for NIFTI format with saved nifti_hdr 
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


%% Subfunctions from dicm2nii, Convert 3x3 direction cosine matrix to quaternion
% Simplied from Quaternions by Przemyslaw Baranski 
function [q, proper] = dcm2quat(R)
% [q, proper] = dcm2quat(R)
% Retrun quaternion abcd from normalized matrix R (3x3)
proper = sign(det(R));
if proper<0, R(:,3) = -R(:,3); end

q = sqrt([1 1 1; 1 -1 -1; -1 1 -1; -1 -1 1] * diag(R) + 1) / 2;
if ~isreal(q(1)), q(1) = 0; end % if trace(R)+1<0, zero it
[mx, ind] = max(q);
mx = mx * 4;

if ind == 1
    q(2) = (R(3,2) - R(2,3)) /mx;
    q(3) = (R(1,3) - R(3,1)) /mx;
    q(4) = (R(2,1) - R(1,2)) /mx;
elseif ind ==  2
    q(1) = (R(3,2) - R(2,3)) /mx;
    q(3) = (R(1,2) + R(2,1)) /mx;
    q(4) = (R(3,1) + R(1,3)) /mx;
elseif ind == 3
    q(1) = (R(1,3) - R(3,1)) /mx;
    q(2) = (R(1,2) + R(2,1)) /mx;
    q(4) = (R(2,3) + R(3,2)) /mx;
elseif ind == 4
    q(1) = (R(2,1) - R(1,2)) /mx;
    q(2) = (R(3,1) + R(1,3)) /mx;
    q(3) = (R(2,3) + R(3,2)) /mx;
end
if q(1)<0, q = -q; end % as MRICron

%% 
function v = normc(M)
vn = sqrt(sum(M .^ 2)); % vn = vecnorm(M);
vn(vn==0) = 1;
v = bsxfun(@rdivide, M, vn);