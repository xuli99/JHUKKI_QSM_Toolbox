function Params = ReadBrukerParams(folderName, procno, Params)

% updated X.L, 2019-09

% Do complicated stuff to get filename.......
if(filesep=='\')
    % Windows
    delimiter = '\\';
else
    % Linux
    delimiter = '/';
end
% Split up
parentFolders = textscan(folderName, '%s', 'delimiter', delimiter);
parentFolders = parentFolders{1};

% Second to last folder is baseName
Params.PathName = folderName;
Params.FileBaseName = parentFolders{length(parentFolders)-1};

% Load other stuff
method_file = fullfile(folderName, 'method');
acqp_file   = fullfile(folderName, 'acqp');
reco_file   = fullfile(folderName, 'pdata', num2str(procno), 'reco');
%% read paramters from method file

para_list = fopen(method_file);
if para_list == -1
    error('Could not open method file');
end
tline = fgetl(para_list);
while ischar(tline)
    
   matches = strfind(tline, '##$PVM_EchoTime=');       % TE for single echo
   num = length(matches);
   if num > 0    
     Params.TE1 = str2double(tline(strfind(tline, '=')+1:end));
   end

   matches = strfind(tline, '##$FirstEchoTime=');       % TE1
   num = length(matches);
   if num > 0    
     Params.TE1 = str2double(tline(strfind(tline, '=')+1:end));
   end

   matches = strfind(tline, '##$EchoSpacing=');         % deltaTE
   num = length(matches);
   if num > 0    
     Params.deltaTE = str2double(tline(strfind(tline, '=')+1:end));
   end   

   matches = strfind(tline, '##$PVM_RepetitionTime=');  % TR
   num = length(matches);
   if num > 0    
     Params.TR = str2double(tline(strfind(tline, '=')+1:end));
   end   

   matches = strfind(tline, '##$PVM_Fov=');             % fov
   num = length(matches);
   if num > 0    
     Params.fov = str2num(fgetl(para_list));         % read next line   
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
   
    matches = strfind(tline, '##$RecoNumInputChan=');   % Number of Input Channel
    num = length(matches);
    if num > 0
      RecoNumInputChan = str2double(tline(strfind(tline, '=')+1:end));
    else
      RecoNumInputChan = 1;                             % default if not recorded
    end    
    
    matches = strfind(tline, '##$RecoCombineMode=');    % Reconstruction Combination Mode
    num = length(matches);
    if num > 0    
      RecoCombineMode = tline(strfind(tline, '=')+1:end);
    else
      RecoCombineMode = 'Normal';                          % default if not recorded  
    end    
    
    tline = fgetl(para_list);
    
end
fclose(para_list);

Params.nchannel = RecoNumInputChan;
Params.nEchoes = NECHOES;
Params.sizeVol = [Nrows, Ncolumns, Nslices];
Params.voxSize = Params.fov./Params.sizeVol;

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

% Figure out TE's
Params.startEcho = Params.TE1;     % in ms
if isfield(Params, 'deltaTE')
    Params.deltaTE = Params.deltaTE;   % in ms
    Params.endEcho = (Params.nEchoes-1)*Params.deltaTE + Params.startEcho;
    Params.TEs = (Params.startEcho:Params.deltaTE:Params.endEcho)./1000;        % in sec
    Params.echoNums = 1:length(Params.TEs);
else
    Params.deltaTE = 0;     % single echo
    Params.TEs = Params.startEcho./1000;
    Params.echoNums = 1;
end

end