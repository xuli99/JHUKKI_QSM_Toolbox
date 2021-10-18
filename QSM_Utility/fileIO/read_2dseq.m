function [ img , hdr ] = read_2dseq(path2dseq)

% function [ img , hdr ] = read_2dseq(path2dseq)
%
% This function reads image data in the Bruker format and return it into
% an array where each pixel is in double format. Also, it output the
% spatial resolution (nom_resol) in mm.
%
% Chern-Chyi (Cecil) Yen @ CMS/LFMI/NINDS/NIH
% CentOS 6.9 & Matlab 2018a

% History:
% 07/21/'06 First version by Fernando Paiva
% 11/28/'06 +IMND support
% 12/07/'06 +3D sequences support
% 12/14/'11 rewrite by Cecil
%           -IMND support
% 12/15/'11 *support multi-echoes data
% 12/16/'11 *incorporate 2dseq reading scripts from Kai-Hsiang Chuang and Chien-Yuan Lin
% 12/19/'11 *support dASL data for Afonso Silva
% 12/20/'11 +convert Unix time to local time and recorded the loading time
% 12/21/'11 +read offset in three directions
% 12/22/'11 +read orientation and gap
% 12/23/'11 +output parameters to screen
% 12/27/'11 *read matrix from RECO_size, because the matrix size will changed by zero-filling during reconstruction
% 12/28/'11 *use switch/case instead of if/elseif to simplify the script
% 12/29/'11 +output more orientation information
% 12/30/'11 *optimize the script by reducing the usage of strfind
% 2/29/'12 *use SpatResol(3) instead of SliceThick for 3D mode and fix error in caculating the fov with gap
% 8/12/'13 +work for RARE
% 8/13/'13 +read shuffled images
% 10/16/'14 +read complex images
% 5/30/'18 *rewrite and tidy the code for marmoset connectome project
% 5/31/'18 *adopt some ideas from Bruker2nifti by Cristina Chavarrias
% 6/1/'18 +read AFI and FieldMap
% 6/4/'18 *use cell structure for output parameters
% 6/5/'18 *use visu_pars
% 6/6''18 +read gradient orientation
% 6/7/'18 *use memmapfile to speed up 2 times
% 6/8/'18 *read parameters by fread to speed up 10 times
% 11/16/'18 *change newline to char(10)

% updated by X.L. to add other relavant header info -- 2019
% updated by X.L. add 2D slice info saved in method, in case it is not saved in visu_par
% updated by X.L., fixed offset/slope parsing for multi-echoes

if nargin < 1
    path2dseq  = uigetdir('Select the directory containing the targeted 2dseq');
end

% tic
fprintf('Reading from %s',path2dseq)

if exist(fullfile(path2dseq,'2dseq'),'file')
    imagefile = fullfile(path2dseq,'2dseq');
else
    error('2dseq does not exist!')
end

if exist([path2dseq,'/../../method'],'file')
    methodfile = [path2dseq,'/../../method'];
else
    error('method does not exist!')
end

if exist([path2dseq,'/../../acqp'],'file')
    acqpfile = [path2dseq,'/../../acqp'];
else
    error('acqp does not exist!')
end

if exist([path2dseq,'/reco'],'file')
    recofile = [path2dseq,'/reco'];
else
    error('reco does not exist!')
end

if exist([path2dseq,'/visu_pars'],'file')
    visu_parsfile = [path2dseq,'/visu_pars'];
else
    error('visu_pars does not exist!')
end

%%%%%%%%%%%%%% Parameter Reading %%%%%%%%%%%%%%
% Default JACMP-DX Header
% ##TITLE=Parameter List
% ##JCAMPDX=4.24
% ##DATATYPE=Parameter Values
% ##ORIGIN=Bruker BioSpin MRI GmbH
% ##OWNER=

fprintf('.')
fp = fopen(methodfile,'r');
method = fread(fp,[1 inf],'*char');
[idx_sta,idx_end] = regexp(method,'##\$\w*='); 
fclose(fp);

for idx = 1:length(idx_sta)
    switch method(idx_sta(idx)+3:idx_end(idx)-1)
        case 'Method' %string: Measuring method
            Method = method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1);
        case 'PVM_NAverages' %integer: Number of times the signal is accumulated prior to the storage on disk and the reconstruction.
            hdr.NAverages = str2double(method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1));
        case 'PVM_NRepetitions' %integer: Number of repetitions (executions) of the experiment.
            NRepetitions = str2double(method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1));
        case 'PVM_RepetitionTime' %(ms)
            RepetitionTime = str2double(method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1));
        case 'PVM_SPackArrReadOrient' %array of strings: Read gradient orientation in the slice package ( L_R / A_P / H_F)
            SPackArrReadOrient = method(idx_end(idx)+find_newline(method,idx_end(idx))+1:idx_end(idx)+find_newline(method,idx_end(idx))+3);
        case 'PVM_SPackArrSliceOrient' %array of strings: General orientation of each slice package (axial, sagittal, coronal)
            SPackArrSliceOrient = sscanf(method(idx_end(idx)+find_newline(method,idx_end(idx))+1:idx_end(idx)+find_newline(method,idx_end(idx))+8),'%s');
        case 'PVM_ObjOrderScheme' %string: Selects the order in which the slices are excited in a multislice experiment. (Interlaced, Reverse_sequential, Sequential)
            ObjOrderScheme = method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1);
        case 'FairMode' %string: Type of FAIR experiment (SELECTIVE, NONSELECTIVE, INTERLEAVED, INTERLEAVED2)
            FairMode = method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1);
        case 'FairTIR_NExp' %integer: Number of different TIR values repeated with each inversion mode (selective and non-selective).
            FairTIR_NExp = str2double(method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1));
        case 'CASL_AcqOrder' %string: (Interleaved, Dynamic)
            CASL_AcqOrder = method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1);
        case 'CASL_LabelImages'
            CASL_LabelImages = str2double(method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1));
        case 'CASL_ControlImages'
            CASL_ControlImages = str2double(method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1));
        case 'MPRAGE_Selection' %string: (MP1RAGE, MP2RAGE, MP3RAGE)
            MPRAGE_Selection = method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1);
        case 'PVM_DwNDiffDir' %integer
            DwNDiffDir = str2double(method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1));
        case 'PVM_ScanTime' %integer
            ScanTime = str2double(method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1));
        case 'PVM_ScanTimeStr' %integer
            ScanTimeHMS = sscanf(method(idx_end(idx)+find_newline(method,idx_end(idx))+1:...
                idx_end(idx)+find_newline(method,idx_end(idx))+find_newline(method, idx_end(idx)+find_newline(method,idx_end(idx))+1)),'<%dh%dm%ds%dms>',4);
            ScanTime = duration(ScanTimeHMS(1:3)');
        case {'PVM_EchoTime'; 'FirstEchoTime'} % in ms
            TE1 = str2double(method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1));
        case 'EchoSpacing'
            deltaTE = str2double(method(idx_end(idx)+1:idx_end(idx)+find_newline(method,idx_end(idx))-1));
        case 'PVM_SPackArrNSlices'
            SlicePacksSlices = sscanf(method(idx_end(idx)+find_newline(method,idx_end(idx))+1:idx_end(idx)+find_newline(method,idx_end(idx))+13*3),'%d', 1);
        case 'PVM_SPackArrSliceDistance'
            SPackArrSliceDistance = sscanf(method(idx_end(idx)+find_newline(method,idx_end(idx))+1:idx_end(idx)+find_newline(method,idx_end(idx))+19*1),'%f',1);
    end
end
fprintf('.')
fp = fopen( acqpfile, 'r' );
acqp = fread(fp,[1 inf],'*char');
[idx_sta,idx_end] = regexp(acqp,'##\$\w*=');
fclose(fp);
for idx = 1:length(idx_sta)
    switch acqp(idx_sta(idx)+3:idx_end(idx)-1)
        case 'ACQ_n_echo_images' %integer: It contains the number of images per slice within the slice loop in 2dseq, typically the number of echo-images.
            NEchoes = str2double(acqp(idx_end(idx)+1:idx_end(idx)+find_newline(acqp,idx_end(idx))-1));
        case 'ACQ_CalibratedRG' %array of number: Linear reciever gain 0.25-203
            hdr.CalibratedRG = sscanf(acqp(idx_end(idx)+find_newline(acqp,idx_end(idx))+1:idx_end(idx)+find_newline(acqp,idx_end(idx))+6*16),'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f',[1 16]);
    end
end
fprintf('.')
fp = fopen(recofile,'r');
reco = fread(fp,[1 inf],'*char');
[idx_sta,idx_end] = regexp(reco,'##\$\w*=');
fclose(fp);

CombineMode = 'SumOfSquares'; % as default

for idx = 1:length(idx_sta)
    switch reco(idx_sta(idx)+3:idx_end(idx)-1)
        case 'RecoCombineMode' %string. either SumOfSquares, AddImages, or ShuffleImages
            CombineMode = reco(idx_end(idx)+1:idx_end(idx)+find_newline(reco,idx_end(idx))-1);
        case 'RecoNumInputChan' %the parameter describes the structure of the input data. When RecoNumInputChan > 1, reconstruction assumes, that the raw data file consists of RecoNumInputChan blocks of size RECO_inp_size[0] forming the first dimension of the data file.
            NumInputChan = str2double(reco(idx_end(idx)+1:idx_end(idx)+find_newline(reco,idx_end(idx))-1));
        case 'RECO_image_type' %string
            image_type = reco(idx_end(idx)+1:idx_end(idx)+find_newline(reco,idx_end(idx))-1);
        case 'RECO_map_offset'
            % Reco_offset = sscanf(reco(idx_end(idx)+find_newline(reco,idx_end(idx))+1:idx_end(idx)+find_newline(reco,idx_end(idx))+16),'%f',1);
            tmp_stra=reco(idx_end(idx)+find_newline(reco,idx_end(idx))+1:idx_sta(idx+1)-1);
            Reco_offset = parse_pattern1(tmp_stra);
        case 'RECO_map_slope'
            % Reco_slope = sscanf(reco(idx_end(idx)+find_newline(reco,idx_end(idx))+1:idx_end(idx)+find_newline(reco,idx_end(idx))+16),'%f',1);
            tmp_stra=reco(idx_end(idx)+find_newline(reco,idx_end(idx))+1:idx_sta(idx+1)-1);
            Reco_slope = parse_pattern1(tmp_stra);            
        case 'RecoScaleChan' % array of numbers
            hdr.RecoScaleChan = sscanf(reco(idx_end(idx)+find_newline(reco,idx_end(idx))+1:idx_end(idx)+find_newline(reco,idx_end(idx))+9*16),'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f',[1 16]);
        case 'RecoPhaseChan' % array of numbers
            hdr.RecoPhaseChan = sscanf(reco(idx_end(idx)+find_newline(reco,idx_end(idx))+1:idx_end(idx)+find_newline(reco,idx_end(idx))+8*16),'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f',[1 16]);
    end
end
fprintf('.')
fp = fopen( visu_parsfile, 'r' );
visu_pars = fread(fp,[1 inf],'*char');
[idx_sta,idx_end] = regexp(visu_pars,'##\$\w*=');
fclose(fp);
for idx = 1:length(idx_sta)
    switch visu_pars(idx_sta(idx)+3:idx_end(idx)-1)
        case 'VisuCoreDim'
            VisuCoreDim = str2double(visu_pars(idx_end(idx)+1:idx_end(idx)+find_newline(visu_pars,idx_end(idx))-1));
        case 'VisuCoreSize'
            VisuCoreSize = sscanf(visu_pars(idx_end(idx)+find_newline(visu_pars,idx_end(idx))+1:idx_end(idx)+find_newline(visu_pars,idx_end(idx))+5*3),'%d %d %d',3);
        case 'VisuCoreExtent'
            VisuCoreExtent = sscanf(visu_pars(idx_end(idx)+find_newline(visu_pars,idx_end(idx))+1:idx_end(idx)+find_newline(visu_pars,idx_end(idx))+19*3),'%f %f %f',3);
        case 'VisuCoreOrientation'
            VisuCoreOrientation = (sscanf(visu_pars(idx_end(idx)+find_newline(visu_pars,idx_end(idx))+1:idx_end(idx)+find_newline(visu_pars,idx_end(idx))+19*9),'%f %f %f %f %f %f %f %f %f',[3 3]))';
        case 'VisuCorePosition'
            VisuCorePosition = sscanf(visu_pars(idx_end(idx)+find_newline(visu_pars,idx_end(idx))+1:idx_end(idx)+find_newline(visu_pars,idx_end(idx))+19*3),'%f %f %f',3);
        case 'VisuCoreDataOffs'
            hdr.scl_inter = sscanf(visu_pars(idx_end(idx)+find_newline(visu_pars,idx_end(idx))+1:idx_end(idx)+find_newline(visu_pars,idx_end(idx))+19*1),'%f',1);
        case 'VisuCoreDataSlope'
            hdr.scl_slope = (sscanf(visu_pars(idx_end(idx)+find_newline(visu_pars,idx_end(idx))+1:idx_end(idx)+find_newline(visu_pars,idx_end(idx))+19*1),'%f',1));
        case 'VisuCoreWordType'
            VisuCoreWordType = visu_pars(idx_end(idx)+1:idx_end(idx)+find_newline(visu_pars,idx_end(idx))-1);
        case 'VisuCoreSlicePacksSlices'
            VisuCoreSlicePacksSlices = sscanf(visu_pars(idx_end(idx)+find_newline(visu_pars,idx_end(idx))+1:idx_end(idx)+find_newline(visu_pars,idx_end(idx))+13*3),'(%d,%d) ',[2 3]);
            VisuCoreSlicePacksSlices = sum(VisuCoreSlicePacksSlices(2,:));
        case 'VisuCoreSlicePacksSliceDist'
            VisuCoreSlicePacksSliceDist = sscanf(visu_pars(idx_end(idx)+find_newline(visu_pars,idx_end(idx))+1:idx_end(idx)+find_newline(visu_pars,idx_end(idx))+19*1),'%f',1);
        case 'VisuSubjectPosition'
            hdr.VisuSubjectPosition = visu_pars(idx_end(idx)+1:idx_end(idx)+find_newline(visu_pars,idx_end(idx))-1);
    end
end

%%%%%%%%%%%%%% Raw Data Reading %%%%%%%%%%%%%%

switch VisuCoreWordType
    case '_32BIT_SGN_INT'
        precision = 'int32';
        hdr.datatype = 8;
        hdr.bitpix = 32;
    case '_16BIT_SGN_INT'
        precision = 'int16';
        hdr.datatype = 4;
        hdr.bitpix = 16;
    case '_8BIT_UNSGN_INT'
        precision = 'uint8';
        hdr.datatype = 2;
        hdr.bitpix = 8;
    case '_32BIT_FLOAT'
        precision = 'single';
        hdr.datatype = 16;
        hdr.bitpix = 32;
end

%PV loops
%   dim1 dim2 dim3 Echo Slice Cycle Repetition Channel Complex

switch VisuCoreDim
    case 2
        hdr.size(1:2) = VisuCoreSize(1:2);
        if ~exist('VisuCoreSlicePacksSlices', 'var')
            VisuCoreSlicePacksSlices = SlicePacksSlices;
        end
        hdr.size(3) = VisuCoreSlicePacksSlices;
        hdr.dim(1:2) = VisuCoreExtent(1:2)./VisuCoreSize(1:2);
        if ~exist('VisuCoreSlicePacksSliceDist', 'var')
            VisuCoreSlicePacksSliceDist = SPackArrSliceDistance;
        end
        hdr.dim(3) = VisuCoreSlicePacksSliceDist;
        VisuCoreSize(3) = 1;
    case 3
        hdr.size = VisuCoreSize;
        hdr.dim(1:3) = VisuCoreExtent./VisuCoreSize;
        VisuCoreSlicePacksSlices = 1;
    otherwise
        error('Data is not 2D or 3D format!')
end

if ~exist('NEchoes','var') || strcmp(Method,'<Bruker:FieldMap>')
    NEchoes = 1;
end

if exist('MPRAGE_Selection','var')
    Cycle = sscanf(MPRAGE_Selection,'MP%dRAGE');
elseif exist('DwNDiffDir','var')
    Cycle = DwNDiffDir;
elseif  exist('FairMode','var') && strcmp(FairMode,'INTERLEAVED')
    Cycle = FairTIR_NExp+1;
elseif (exist('CASL_AcqOrder','var') && strcmp(CASL_AcqOrder,'INTERLEAVED'))
    Cycle = CASL_LabelImages+CASL_ControlImages;
elseif (strcmp(Method,'<Bruker:FieldMap>') && strcmp(CombineMode,'ShuffleImages')) || strcmp(Method,'<User:rpAFI>')
    Cycle = 2;
else
    Cycle = 1;
end

if ~strcmp(CombineMode,'ShuffleImages')
    NumInputChan = 1;
end

if ~exist('NRepetitions','var')
    NRepetitions = 1;
end
fprintf('.')
mem_img = memmapfile(imagefile,'Format',{precision,[VisuCoreSize(1),VisuCoreSize(2),VisuCoreSize(3),NEchoes,VisuCoreSlicePacksSlices,Cycle,NRepetitions,NumInputChan],'raw'});
if strcmp(image_type,'COMPLEX_IMAGE')
    img = complex(mem_img.Data(1).raw,mem_img.Data(2).raw);
    hdr.datatype = 32; %Requiring data to convert to single (float32)
    hdr.bitpix = 64;
else
    img = mem_img.Data.raw;
end

% change to float and do scaling
img = cast(img, 'single');
img = img./Reco_slope + Reco_offset; 
hdr.datatype = 32; %Requiring data to convert to single (float32)
hdr.bitpix = 64;

% Writing headers
hdr.method = Method;
hdr.ReadOrient = SPackArrReadOrient;
hdr.SliceOrient = SPackArrSliceOrient;
hdr.ObjOrderScheme = ObjOrderScheme;
hdr.ScanTime = ScanTime;
hdr.RepetitionTime = RepetitionTime;
hdr.VisuCoreOrientation = VisuCoreOrientation;
hdr.VisuCorePosition = VisuCorePosition;
hdr.VisuCoreDim = VisuCoreDim;
hdr.NEchoes = NEchoes;
hdr.image_type = image_type;
hdr.NumInputChan = NumInputChan;
hdr.NRepetitions = NRepetitions;

% Get TEs
if (NEchoes > 1) && ~isempty('deltaTE')
    TEn = (NEchoes-1)*deltaTE + TE1;
    hdr.TEs = (TE1:deltaTE:TEn)./1000;        % in sec
    hdr.echoNums = 1:NEchoes;
else
    hdr.TEs = TE1./1000;
    hdr.echoNums = 1;
end

% fprintf('read_2dseq finished in %5.2f sec\n',toc);
fprintf('done\n');

function idx_out = find_newline(str,idx_in)
    % Find newline after certain index
    idx_out = 1;
    while str(idx_in+idx_out)~=char(10)
        idx_out = idx_out+1;
    end
end

function data_out = parse_pattern1(str_in)
    str_list_new = splitlines(str_in);
    if contains(str_list_new{1},'@')
        str_list_new_split=split(str_list_new{1},{'(',')'});
        data_out = str2double(str_list_new_split{2});  % inside ()
    else
        str_list_new_split=split(str_list_new{1},' ');
        data_out= str2double(str_list_new_split{1});   % the first one
    end
end

end