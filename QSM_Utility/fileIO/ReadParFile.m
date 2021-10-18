function header = ReadParFile(filenamePar, Params)
%%
% Author: Xu Li
% Extended by Jiri van Bergen
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu

% updated for v3.1 due to .par header change, by X.L. 202001

%% extract Angulation Parameters from Philips par file
% for Philips system ver 4.2 and above

% Basics
header = Params;

% Basics
[header.PathName, header.FileBaseName,~] = fileparts(filenamePar);

% Open file
fid = fopen(filenamePar);

% Header/rest
ParamsFlag = 0;

% read in the parameters from .par file
while (~feof(fid))
    
    tline = fgetl(fid);
    if isempty(tline)
        continue;
    end
    
    if (tline(1) == '#')        % header part
        isheader = 1;
        ispara_block = 0;
        isslice_block = 0;
        
    elseif (tline(1) == '.')    % basic parameter part
        isheader = 0;
        ispara_block = 1;
        isslice_block = 0;
        
        %  Look for number of slices
        if(contains(tline, 'number of slices', 'IgnoreCase', true))
            aa = tline(strfind(tline, ':')+1:end);
            aa = str2num(aa);
            header.nslices = aa(1);
        end
        
        %  Look for the FOV (ap, fh, rl)
        if(contains(tline, 'FOV (', 'IgnoreCase', true))
            aa = tline(strfind(tline, ':')+1:end);
            header.fov = str2num(aa);
        end
        
        %  Look for the patient position
        if(contains(tline, 'Patient position', 'IgnoreCase', true))
            aa = tline(strfind(tline, ':')+1:end);
            aa = aa(4:end);  
            aaind = strfind(aa, ' ');           
            if isempty(aaind)
                header.patient_position = aa(1:2); % brief info, e.g. HFS
                header.patient_orientation = aa(3);
            else               
                header.patient_position = [aa(1), aa(aaind(1)+1)];
                header.patient_orientation = aa(aaind(2)+1:end);
            end
        end
        
        %  Look for the Angulation Parameters (ap, fh, rl)
        if(contains(tline, 'Angulation midslice(', 'IgnoreCase', true))
            aa = tline(strfind(tline, ':')+1:end);
            header.ang = str2num(aa);
        end
        
        %  Look for number of echos
        if(contains(tline, 'number of echoes', 'IgnoreCase', true))
            aa = tline(strfind(tline, ':')+1:end);
            aa = str2num(aa);
            header.nEchoes = aa(1);
        end
        
        %  Look for number of dynamics
        if(contains(tline, 'number of dynamics', 'IgnoreCase', true))
            aa = tline(strfind(tline, ':')+1:end);
            aa = str2num(aa);
            header.nDynamics = aa(1);
        end
        
        %  Look for TR
        if(contains(tline, 'Repetition time [', 'IgnoreCase', true))
            aa = tline(strfind(tline, ':')+1:end);
            aa = str2num(aa);
            header.TR = aa(1);
        end
    else                        % slice parameter part
        isheader = 0;
        ispara_block = 0;
        isslice_block = 1;
        if contains(tline, 'SRT')
            sind1 = strfind(tline, 'SRT');
            sind2 = strfind(tline, 'route');
            tline(sind1:sind2+4)=[];        
        end        
        aa = str2num(tline);
        
        if ParamsFlag == 0      % have not got the required Parameters
            % SliceOri           
            header.sliceOri = aa(26);       % TRA/SAG/COR
            
            % Size
            header.sizeVol(1) = aa(10);     
            header.sizeVol(2) = aa(11);
                        
            % Start echo
            header.startEcho = aa(31);
            header.TEs = aa(31);
            
            % resolution
            header.voxSize(3) = aa(23);     
            header.voxSize(1) = aa(29);     
            header.voxSize(2) = aa(30);
            
            header.slicegap = aa(24);       
            header.voxSize(3) = header.voxSize(3) + header.slicegap;   
            
            ParamsFlag = 1;
        else
            
            if(aa(31) ~= header.TEs(end) && length(header.TEs) < header.nEchoes)
                 header.TEs = [header.TEs aa(31)];
            end
        end 
    end
end
fclose(fid);

Tpo = []; Tpp = [];
if strcmpi(header.patient_orientation,'supine') || (strcmpi(header.patient_orientation,'s'))
    Tpo = [1,0,0;0,1,0;0,0,1];
    rev_Tpo = [1,0,0;0,1,0;0,0,1];
elseif strcmpi(header.patient_orientation,'prone') || (strcmpi(header.patient_orientation,'p'))
    Tpo = [-1,0,0;0,-1,0;0,0,1];
    rev_Tpo = [-1,0,0;0,-1,0;0,0,1];  
elseif strcmpi(header.patient_orientation,'rd')
    Tpo = [0,-1,0;1,0,0;0,0,1];
    rev_Tpo = [0,1,0;-1,0,0;0,0,1];
elseif strcmpi(header.patient_orientation,'ld')
    Tpo = [0,1,0;-1,0,0;0,0,1];
    rev_Tpo = [0,-1,0;1,0,0;0,0,1];
end

if strcmpi(header.patient_position,'ff')
    Tpp = [0,-1,0;-1,0,0;0,0,1];
    rev_Tpp = [0,-1,0;-1,0,0;0,0,-1];
elseif strcmpi(header.patient_position,'hf')
    Tpp = [0,1,0;-1,0,0;0,0,-1];
    rev_Tpp = [0,-1,0;1,0,0;0,0,-1];
end

if ~isempty(Tpo) && ~isempty(Tpp)
    header.Tpom = Tpo*Tpp;
    header.Tpominv = rev_Tpp*rev_Tpo;
end

% for slice orientation
switch header.sliceOri 
    case 1  % TRA
        header.fov = [header.fov(3), header.fov(1), header.fov(2)]; 
        header.Tsom = [0,-1,0; -1,0,0; 0,0,1];                      
        header.Tsominv = [0,-1,0;-1,0,0;0,0,1];
    
    case 2  % SAG
        header.fov = [header.fov(1), header.fov(2), header.fov(3)];
        header.Tsom = [0,0,-1;0,-1,0;1,0,0];
        header.Tsominv = [0,0,1;0,-1,0;-1,0,0];
        
    case 3  % COR
        header.fov = [header.fov(3), header.fov(2), header.fov(1)];
        header.Tsom = [0,-1,0;0,0,1;1,0,0];
        header.Tsominv = [0,0,1;-1,0,0;0,1,0];
    
    otherwise 
        error('Wrong slice orientation parameter ...')
end

header.AngAP = header.ang(1);   % degree
header.AngFH = header.ang(2);   % degree    
header.AngRL = header.ang(3);   % degree 

header.TAng = Rmatrix_ang(header.AngRL, header.AngAP, header.AngFH, 0);     % rotation matrix
header.TAnginv = Rmatrix_ang(header.AngRL, header.AngAP, header.AngFH, 1);

header.TEs          = (header.TEs)./1000;           

header.sizeVol(3)   = header.nslices;
header.sizeRecAll   = [header.sizeVol(1), header.sizeVol(2), header.sizeVol(3), header.nEchoes];

header.fov = header.voxSize.*header.sizeVol;