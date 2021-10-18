function Params = params_bruker(fname, expno, procno)
% Params = params_bruker(fname, expno, procno)
% 
%% Author: Xu Li (xuli@mri.jhu.edu)
% Affiliation: Radiology @ JHU
% 
% read in scan parameters from Bruker ParaVision 5
% the image file is in the directory <fname>/<expno>/pdata/<procno>/2dseq
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
%   parameters extracted from 'method': TE1, deltaTE, TR, fov
%   parameters extracted from 'acqp':   NI, NR, NSLICES, NECHOES
%   parameters extracted from 'reco':   RECOSIZE, RecoCombineMode, RecoNumInputChan 

method_file = [ fname '/' num2str(expno) '/method'];
acqp_file = [ fname '/' num2str(expno) '/acqp'];
reco_file = [ fname '/' num2str(expno) '/pdata/', num2str(procno), '/reco'];

%% read paramters from method file

para_list = fopen(method_file);
if para_list == -1;
    error('Could not open method file');
end
tline = fgetl(para_list);
while ischar(tline)
    
   matches = strfind(tline, '##$PVM_EchoTime=');       % TE for single echo
   num = length(matches);
   if num > 0    
     Params.TE1 = str2num(tline(strfind(tline, '=')+1:end));
   end

   matches = strfind(tline, '##$FirstEchoTime=');       % TE1
   num = length(matches);
   if num > 0    
     Params.TE1 = str2num(tline(strfind(tline, '=')+1:end));
   end

   matches = strfind(tline, '##$EchoSpacing=');         % deltaTE
   num = length(matches);
   if num > 0    
     Params.deltaTE = str2num(tline(strfind(tline, '=')+1:end));
   end   

   matches = strfind(tline, '##$PVM_RepetitionTime=');  % TR
   num = length(matches);
   if num > 0    
     Params.TR = str2num(tline(strfind(tline, '=')+1:end));
   end   

   matches = strfind(tline, '##$PVM_Fov=');             % fov
   num = length(matches);
   if num > 0    
     Params.fov = str2num(fgetl(para_list));
   end      
   
   tline = fgetl(para_list);
end
fclose(para_list);

%% read paramters from acqp file

para_list = fopen(acqp_file);
if para_list == -1;
    error('Could not open acqp file');
end
tline = fgetl(para_list);
while ischar(tline)
    matches = strfind(tline, '##$NI=');                     % NI: Nslice*Nechoes
    num = length(matches);
    if num > 0
      NI = str2num(tline(strfind(tline, '=')+1:end));
    end

    matches = strfind(tline, '##$NR=');                     % NR: Number of Repetition
    num = length(matches);
    if num > 0
      NR = str2num(tline(strfind(tline, '=')+1:end));
    end
    
    matches = strfind(tline, '##$NSLICES=');                % Number of Slices
    num = length(matches);
    if num > 0
      NSLICES = str2num(tline(strfind(tline, '=')+1:end));
      Nslices = NSLICES;                                    % for 2D sequence
    end
    
    matches = strfind(tline, '##$NECHOES=');                % Number of Echoes
    num = length(matches);
    if num > 0
      NECHOES = str2num(tline(strfind(tline, '=')+1:end));
    end
        
    tline = fgetl(para_list);
end
fclose(para_list);

%% read paramters from reco file

para_list = fopen(reco_file);
if para_list == -1;
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
   
    matches = strfind(tline, '##$RecoNumInputChan=');   % Number of Input Channel
    num = length(matches);
    if num > 0
      RecoNumInputChan = str2num(tline(strfind(tline, '=')+1:end));
    end    
    
    matches = strfind(tline, '##$RecoCombineMode=');    % Reconstruction Combination Mode
    num = length(matches);
    if num > 0    
      RecoCombineMode = tline(strfind(tline, '=')+1:end);
    end    
    
    tline = fgetl(para_list);
    
end
fclose(para_list);

Params.nchannel = RecoNumInputChan;
Params.nechos = NECHOES;
Params.sizeVol = [Nrows, Ncolumns, Nslices];

if strcmpi(RecoCombineMode, 'ShuffleImages')
    if (NR > 1) && (NECHOES > 1)
        Params.sizeRecAll = [Nrows, Ncolumns, Nslices, NECHOES, NR, RecoNumInputChan];
    elseif (NECHOES > 1)
        Params.sizeRecAll = [Nrows, Ncolumns, Nslices, NECHOES, RecoNumInputChan];
    else
        Params.sizeRecAll = [Nrows, Ncolumns, Nslices, RecoNumInputChan];       
    end
else
    if (NR > 1) && (NECHOES > 1)
        Params.sizeRecAll = [Nrows, Ncolumns, Nslices, NECHOES, NR];
    elseif (NECHOES > 1)
        Params.sizeRecAll = [Nrows, Ncolumns, Nslices, NECHOES];
    else
        Params.sizeRecAll = [Nrows, Ncolumns, Nslices];          
    end
end

end