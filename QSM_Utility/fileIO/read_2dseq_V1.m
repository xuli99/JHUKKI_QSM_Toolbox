function img = read_2dseq_V1(folderName, procno)
% img = read_2dseq_V1(fname, expno, procno)
% 
%% Author: Jiadi Xu
%% Author: Xu Li (xuli@mri.jhu.edu)
% Affiliation: Radiology @ JHU
% 
% read image data in file "2dseq" reconstructed from Bruker ParaVision 5
% the file is in the directory <fname>/<expno>/pdata/<procno>/2dseq
% 
% certain acquision and reconstruction parameters are stored in the
% following files: 
%  <fname>/<expno>/acqp:    acquisition parameter
%  <fname>/<expno>/fid:     raw data
%  <fname>/<expno>/method:  PVM method paramter
%  <fname>/<expno>/pdata/<procno>/reco:  reconstruction parameter
% 
%   Image Size = sizeof(word)*NR*NI*RECO_size
%       where
%   NR: repetition time (dynamics)
%   NI: slice number X echo number
%   RECO_size: Nx*Ny*Nz
%
%   parameters extracted from 'method': READDIRECTION
%   parameters extracted from 'acqp':   NI, NR, NSLICES, NECHOES
%   parameters extracted from 'reco':   RECOSIZE, WordType, ByteOrder, 
%                                       MapMode, MapOffset, MapSlope

% updated X.L, 2019-09

method_file = fullfile(folderName, 'method');
acqp_file   = fullfile(folderName, 'acqp');
reco_file   = fullfile(folderName, 'pdata', num2str(procno), 'reco');
img_name    = fullfile(folderName, 'pdata', num2str(procno), '2dseq');

%% read paramters from method file

para_list = fopen(method_file);
if para_list == -1
    error('Could not open method file');
end
tline = fgetl(para_list);
while ischar(tline)
    
   matches = strfind(tline, '##$PVM_SPackArrReadOrient');  % Read direction
   num = length(matches);
   if num > 0    
     READDIRECTION = strtrim(fgetl(para_list));
   end
   
   tline = fgetl(para_list);
end
fclose(para_list);

%% read paramters from acqp file

para_list = fopen(acqp_file);
if para_list == -1
    error('Could not open acqp file');
end
tline = fgetl(para_list);
while ischar(tline)
    matches = strfind(tline, '##$NI=');                     % NI: Nslice*Nechoes
    num = length(matches);
    if num > 0
      NI = str2double(tline(strfind(tline, '=')+1:end));
    end

    matches = strfind(tline, '##$NR=');                     % NR: Number of Repetition
    num = length(matches);
    if num > 0
      NR = str2double(tline(strfind(tline, '=')+1:end));
    end
    
    matches = strfind(tline, '##$NSLICES=');                % Number of Slices
    num = length(matches);
    if num > 0
      NSLICES = str2double(tline(strfind(tline, '=')+1:end));
      Nslices = NSLICES;                                    % for 2D sequence
    end
    
    matches = strfind(tline, '##$NECHOES=');                % Number of Echoes
    num = length(matches);
    if num > 0
      NECHOES = str2double(tline(strfind(tline, '=')+1:end));
    end
        
    tline = fgetl(para_list);
end
fclose(para_list);

%% read paramters from reco file

para_list = fopen(reco_file);
if para_list == -1
   error('Could not open reco file');
end
tline = fgetl(para_list);
while ischar(tline)
    matches = strfind(tline, '##$RECO_size=');          % RECOSIZE
    num = length(matches);
    if num > 0       
       RECOSIZE = str2num(fgetl(para_list));
       Nrows = RECOSIZE(1);
       Ncolumns = RECOSIZE(2);
       if length(RECOSIZE) > 2       
            Nslices = RECOSIZE(3);                      % for 3D sequence
       end
    end
   
    matches = strfind(tline, '##$RECO_wordtype=_');     % RECO_wordtype
    num = length(matches);
    if num > 0
       WORDTYPE = tline(strfind(tline, '=')+2:end);
       BIT = ['int', WORDTYPE(1:2)];       
    end
    
    matches = strfind(tline, '##$RECO_byte_order=');     % RECO_byteorder, litlle or big endian
    num = length(matches);
    if num > 0
       BYTEORDER = tline(strfind(tline, '=')+1:end);  
    end
    
    matches = strfind(tline, '##$RECO_map_mode=');          % RECO_map_mode
    num = length(matches);
    if num > 0
      RECOMAPMODE = tline(strfind(tline, '=')+1:end);
    end
    
    matches = strfind(tline, '##$RECO_map_offset=');         % MapOffset
    num = length(matches);
    if num > 0       
       RECOMAPOFFSET = str2num(fgetl(para_list));
       RECOMAPOFFSET = RECOMAPOFFSET(1);
    end
    
    matches = strfind(tline, '##$RECO_map_slope=');          % MapSlope
    num = length(matches);
    if num > 0       
       RECOMAPSLOPE = str2num(fgetl(para_list));
       RECOMAPSLOPE = RECOMAPSLOPE(1);
    end    
    
    matches = strfind(tline, '##$RecoNumInputChan=');       % Number of Input Channel
    num = length(matches);
    if num > 0
        RecoNumInputChan = str2double(tline(strfind(tline, '=')+1:end));
    else
        RecoNumInputChan = 1;
    end    
    
    matches = strfind(tline, '##$RecoCombineMode=');        % Reconstruction Combination Mode
    num = length(matches);
    if num > 0    
        RecoCombineMode = tline(strfind(tline, '=')+1:end);
    else
        RecoCombineMode = 'Normal';
    end   
    
    tline = fgetl(para_list);
end
fclose(para_list);


%% read image data

Bruker_recon = fopen(img_name,'r', BYTEORDER(1));
if Bruker_recon==-1
    error('Could not open 2dseq file');
end

if strcmpi(RecoCombineMode, 'ShuffleImages')
    if (NR > 1) && (NECHOES > 1)
        sizeRecAll = [Nrows, Ncolumns, Nslices, NECHOES, NR, RecoNumInputChan];
        sizeD3 = Nslices;
        sizeD4 = NECHOES;
        sizeD5 = NR;
        sizeD6 = RecoNumInputChan;

    elseif (NECHOES > 1)
        sizeRecAll = [Nrows, Ncolumns, Nslices, NECHOES, RecoNumInputChan];
        sizeD3 = Nslices;
        sizeD4 = NECHOES;
        sizeD5 = RecoNumInputChan;
        sizeD6 = 1;
    else
        sizeRecAll = [Nrows, Ncolumns, Nslices, RecoNumInputChan];
        sizeD3 = Nslices;
        sizeD4 = RecoNumInputChan;
        sizeD5 = 1;
        sizeD6 = 1;          
    end
else
    if (NR > 1) && (NECHOES > 1)
        sizeRecAll = [Nrows, Ncolumns, Nslices, NECHOES, NR];
        sizeD3 = Nslices;
        sizeD4 = NECHOES;
        sizeD5 = NR;
        sizeD6 = 1;

    elseif (NECHOES > 1)
        sizeRecAll = [Nrows, Ncolumns, Nslices, NECHOES];
        sizeD3 = Nslices;
        sizeD4 = NECHOES;
        sizeD5 = 1;
        sizeD6 = 1;
    else
        sizeRecAll = [Nrows, Ncolumns, Nslices];
        sizeD3 = Nslices;
        sizeD4 = 1;
        sizeD5 = 1;
        sizeD6 = 1;          
    end
end


img = zeros(sizeRecAll);

h = waitbar(0, ['Reading ', num2str(img_name)]);

SliceInd = 0;
TotalSlice = prod(sizeRecAll(3:end));

for ii = 1:sizeD6
    for jj = 1:sizeD5
        for kk = 1:sizeD4
            for ll = 1:sizeD3
                img(:,:,ll,kk,jj,ii) = fread(Bruker_recon, [Nrows, Ncolumns], BIT);    % read in slice by slice
                SliceInd = SliceInd + 1;
                waitbar(SliceInd/TotalSlice);
            end
        end        
    end
end

close(h)

fclose(Bruker_recon);


%% Scaling and permutation
switch READDIRECTION
    case 'L_R'
        dopermute = 0;      % right? need to confirm.
    case 'A_P'
        dopermute = 1;
    case 'H_F'
        dopermute = 0;
    otherwise
        error('check Read direction, only L_R and A_P allowed');
end

switch RECOMAPMODE
    case 'ABSOLUTE_MAPPING'
        img = img./RECOMAPSLOPE + RECOMAPOFFSET;          % absolute mapping scaling
        
    otherwise
        error('check RECO map mode');
end
      
if dopermute
    img = permute(img, [2, 1, 3:length(size(img))]);        % permute
end

end