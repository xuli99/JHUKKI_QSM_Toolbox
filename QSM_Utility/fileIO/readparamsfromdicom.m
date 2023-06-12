function Params = readparamsfromdicom(dicomheader, Params)
% function Params = readparamsfromdicom(dicomheader, Params)
%
%% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
%
% extracted scan information from dicom header and save to Params structure
%
% input:
%        dicomheader: dicomheader extracted from enhanced DICOM format
%        Params: The basic params structure
% output:
%       Params.FileBaseName
%       Params.nEchoes  = NumEcho;
%       Params.voxSize
%       Params.sizeVol
%       Params.B0: main field strength
%       Params.TE echo time (array if multi-echo, TE1 otherwise)
%       Params.TAng
%       Params.TR
%
%       Other optional info
%       Params.DynamicNum
%       Params.coilNum % added default 2019-10-07

% created on 2019-07-09
% bug fix on 2023-01-30

pfgs = dicomheader.PerFrameFunctionalGroupsSequence;
items = fieldnames(pfgs);
nimag = length(items);    

Params.B0 = dicomheader.MagneticFieldStrength;

% change data format to single, otherwise it is uint16
Params.sizeVol(1)   = cast(dicomheader.Width, 'double');
Params.sizeVol(2)   = cast(dicomheader.Height, 'double');
Params.sizeVol(3)   = cast(dicomheader.MRSeriesNrOfSlices, 'double');

Params.voxSize(1)   = cast(pfgs.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing(1), 'double');       % mm
Params.voxSize(2)   = cast(pfgs.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing(2), 'double');
Params.voxSize(3)   = cast(pfgs.Item_1.PixelMeasuresSequence.Item_1.SliceThickness, 'double');        % mm

Params.fov      = (round(Params.voxSize.*double(Params.sizeVol)));
Params.nEchoes     = dicomheader.MRSeriesNrOfEchoes;
Params.TR          = dicomheader.MRSeriesRepetitionTime(1);
Params.nDynamics   = dicomheader.MRSeriesNrOfDynamicScans;

if ~isfield(dicomheader, 'coilNum')
    Params.coilNum     = 1;     % default
end

% Angulation matrix (same for all the images)
Affine2D = reshape(dicomheader.PerFrameFunctionalGroupsSequence.Item_1.PlaneOrientationSequence.Item_1.ImageOrientationPatient,[3 2]);
Params.TAng = [Affine2D, cross(Affine2D(:,1), Affine2D(:,2))];

% SliceOriSave
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

% TEs
TE = zeros([Params.nEchoes 1]);
for ii = 1:nimag    
    imagecho = pfgs.(items{ii}).PrivatePerFrameSq.Item_1.EchoNumber;
    if TE(imagecho) == 0  % needs update
        TE(imagecho) = pfgs.(items{ii}).PrivatePerFrameSq.Item_1.EchoTime;
    end
end
Params.TEs          = TE./1000;     % in sec

