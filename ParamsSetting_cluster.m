% Parameters setting template 

% -----------------------  Parameters set in GUI or script ------------------------
% parameters don't need to change
handles.Params.saveOutput           = 1;              % default
handles.Params.LapPhaseCorrection   = 0;              % default

% Method Dictionary
handles.Params.UnwrappingMethodsDict    = {'Path', 'Laplacian', 'NonlinearFit + Path'}';
handles.Params.BgRemovalMethodsDict     = {'VSHARP', 'PDF', 'LBV+VSHARP', 'iRSHARP'}';
handles.Params.QSMSolverDict            = {'iLSQR','TKD','iTKD','MEDI','SFCR','nSFCR','FANSI','NDI','TFI'}';

% -----------parameters may need changes
handles.Params.thresh_tsvd          = 0.05;           % default, good trade-off

% -----------parameter input from GUI
% setup B0
handles.Params.B0                   = 3;              % Tesla

% Method selection
handles.VarUnwrappingMethod.Value   = 1;           % default, select according to Dict
handles.VarBgRemoval.Value          = 3;           % default, select according to Dict
handles.VarQSMSolver.Value          = 5;           % default, select according to Dict

% Phase pre-processing
handles.Params.UnwrappingMethod     = handles.VarUnwrappingMethod.Value;
handles.Params.phase2DprocFlag      = 0;           % default, edit, 0-1 
handles.Params.TemplateEcho         = 1;           % default, edit, 0-1, using template unwrapping

% Echo selection
handles.Params.echoStart            = 3;           % Edit for selecting starting echo
handles.Params.echoStep             = 1;           % default, can edit, 0 is single echo
handles.Params.echoEnd              = 6;           % Edit for selecting end echo
if handles.Params.echoStep > 0
    handles.Params.echoNums         = (handles.Params.echoStart):(handles.Params.echoStep):(handles.Params.echoEnd);
else
    handles.Params.echoNums         = handles.Params.echoStart;
end

% for BrainMask
% handles.Params.FSLFolder          = '/usr/local/fsl/bin';   % in case needs to change
handles.Params.SaveEcho             = 1;              % default, edit, '1', or '[1,3]'  
handles.Params.FSLThreshold         = 0.4;            % default, edit, 0-1
handles.Params.ErodeRadius          = 1;              % default, edit, in mm
handles.Params.UnreliThreshold      = 2;              % default, edit, 0-2, threshold for detecting unreliable phase
handles.Params.FSLThreshold         = num2str(handles.Params.FSLThreshold);
handles.Params.unrelyPhase1_thresh  = 0.1;            % before BR, default 0.5
handles.Params.unrelyPhase2_thresh  = 0.1;            % after BR, default 0.5

% Echo Average
handles.Params.EchoAvg              = 1;              % default, edit, 0/1
handles.Params.AvgBeforeBR          = 1;              % default, edit, 0/1 

% Background Removal
handles.Params.BgRemoval            = handles.VarBgRemoval.Value;
handles.Params.SHARPradius          = 6;              % default, edit, in mm

% QSM
handles.Params.QSMSolver            = handles.VarQSMSolver.Value;
handles.Params.R2starFlag           = 0;              % default, edit, 0/1
handles.Params.AutoRefFlag          = 0;              % default, edit, 0/1

% % % --------------- Edit data filenames or data list
% % % --------------- Load in a test data
% PathName_1                          = '../data/exampledata/001';
% FileName_1                          = '001.par';
% handles = OpenFiles_cluster(fullfile(PathName_1, FileName_1), handles);

% % --------------- Load in Processing list for batch processing
handles.DataListFile                = fullfile('./DataTable_template.mat');
handles.TableDatasets.Data = [];
handles.CurrentDataset = 0;
handles = LoadDataList_cluster(handles.DataListFile, handles);