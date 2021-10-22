%% Constant Parameters for QSM_Toolbox_v3.0
% Authors: Jiri van Bergen and Xu Li
% Updated by Xu Li, 2019-06-20
% Updated by Xu Li, 2021-06-23 for cluster version

handles.Params.QSMdir           = fileparts(mfilename('fullpath'));
handles.Params.QSMSettingsFile  = fullfile(handles.Params.QSMdir, '/QSM_ConstantsSaved.mat');

% check FSL installation
fsldir = getenv('FSLDIR');
if isempty(fsldir)
    % Use structure "handles" instead of graphic object handles
    disptxt = 'FSLDIR set to default: /usr/local/fsl, check fsl installation... ';
    disp(disptxt)
    if isfield(handles.Params, 'cluster') % for cluster
        writelog(logfile, [logtxt, '\n']);
    end
    setenv('FSLDIR','/usr/local/fsl');   
    setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
    fsldir = getenv('FSLDIR');
end        

% Check if file exists
if(exist(handles.Params.QSMSettingsFile, 'file') == 2 )
    % Advice message
    disptxt1 = 'Loading previouly saved parameters ... ';
    disptxt2 = 'FOR DEFAULT PARAMETERS REMOVE THE  QSM_ConstantsSaved.mat  FILE ';
    disp(disptxt1); disp(disptxt2);
    if isfield(handles.Params, 'cluster') % for cluster
        writelog(logfile, [disptxt1, '\n']); writelog(logfile, [disptxt2, '\n']);
    end
    % Load Params, only load those defined by user
    load(handles.Params.QSMSettingsFile);            
    handles.Params.saveOutput       = Params.saveOutput;
    handles.Params.B0               = Params.B0;
    
    handles.Params.FSLFolder        = fullfile(fsldir, 'bin/'); % system dependent
    handles.Params.FSLThreshold     = Params.FSLThreshold;
    
    handles.Params.EchoAvg          = Params.EchoAvg;
    handles.Params.ErodeRadius      = Params.ErodeRadius;   % mm
    handles.Params.SHARPradius      = Params.SHARPradius;   % mm 
    handles.Params.MaskThreshold    = Params.MaskThreshold;
    clear Params
else
    % Basic params for GUI - default or User defined
    handles.Params.saveOutput       = 1;
    handles.Params.B0               = '3';  % Tesla    
    handles.Params.FSLFolder        = fullfile(fsldir, 'bin/'); 
    handles.Params.FSLThreshold     = 0.4;  % default BET threshold
    
    handles.Params.EchoAvg          = 1;    % defualt using EchoAvg
    handles.Params.ErodeRadius      = 1;    % mm, mask erosion
    handles.Params.SHARPradius      = 8;    % mm, SHARP kernal radius
    handles.Params.MaskThreshold    = 65;   % backup masking code
end

%% Parameters not defined by Users
% Echo limits for EchoAvg
handles.TELowerLimit_1p5T       = 0/1000; % s
handles.TEUpperLimit_1p5T       = 120/1000; % s
handles.TELowerLimit_3T         = 0/1000; % s 
handles.TEUpperLimit_3T         = 60/1000; % s
handles.TELowerLimit_7T         = 0/1000; % s
handles.TEUpperLimit_7T         = 40/1000; % s
handles.TELowerLimit_11p7T      = 0/1000; % s
handles.TEUpperLimit_11p7T      = 30/1000; % s   

% Other Params
handles.Params.thresh_tsvd      = 0.05;             % good trade-off
handles.Params.gamma            = 42.57747892e6;    % Hz/T, updated
handles.Params.fileTypes        = 'method;*.par;*.PAR;*.DIC;*.dic;*.IMA;*.ima;*.DICOM;*.dicom;*.dcm;*.mat;*.1';

% Strings
handles.textReadyLoad           = 'Ready to be reconstructed';
handles.textReconstructing      = 'Reconstructing...';

% Save hObject if called from GUI (when createing VarB0)
try 
    if strcmp(hObject.Tag, 'VarB0')
        guidata(hObject, handles);
    end
catch ME
    if ((strcmp(ME.identifier,'MATLAB:undefinedVarOrClass')) && (contains(ME.message,'hObject.Tag')))
        writelog(logfile, ['cluster version.', '\n']);
    else
        rethrow(ME)
    end
end
