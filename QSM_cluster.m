function [status] = QSM_cluster(varargin)
% this is the cluster version with no gui interaction
% status = QSM_cluster(varargin)
% Example call: 
%   ParamsSetFile = 'ParamsSetting_cluster.m';
%   LogFile = 'logd.txt';
%   QSM_cluster(ParamsSetFile, LogFile)

% 2024-05-31, update with TableDataFile argument

if nargin < 1
    ParamsFile = 'ParamsSetting_cluster.m';
    logfile = fullfile(fileparts(mfilename('fullpath')), 'logd.txt');
    logtxt = 'Run default parameters setting script: ParamsSetting_cluster.m and log to logd.txt';
elseif nargin < 2
    ParamsFile = varargin{1};
    logfile = fullfile(fileparts(mfilename('fullpath')), 'logd.txt');
    logtxt = ['Run parameters setting script:', strrep(ParamsFile, '\', '\\'), ' and log to logd.txt'];
elseif nargin < 3
    ParamsFile = varargin{1};
    logfile = varargin{2};
    logtxt = ['Run parameters setting script:', strrep(ParamsFile, '\', '\\'), ' and log to ', strrep(logfile, '\', '\\')];
elseif nargin < 4
    ParamsFile = varargin{1};
    logfile = varargin{2};
    logtxt = ['Run parameters setting script:', strrep(ParamsFile, '\', '\\'), ' and log to ', strrep(logfile, '\', '\\')];

    TableDataFile = varargin{3};    % set dataTable separately from ParamsFile
    handles.DataListFile = TableDataFile;
end
disp(logtxt);

% Add modules
FilePath = fileparts(mfilename('fullpath'));
addpath(FilePath);
addpath(genpath(fullfile(FilePath, 'QSM_Modules')));
addpath(genpath(fullfile(FilePath, 'QSM_NIFTI')));
addpath(genpath(fullfile(FilePath, 'QSM_Utility')));

% Start logging
writelog(logfile, [logtxt, '\n'], 'w');     % discard previous logs

% Use structure "handles" instead of graphic object handles
handles.Params.cluster = 1;       % flag the cluster version without GUI
handles.logfile = logfile;        % save logfile

% setting up constants
Constants;

try
    run(ParamsFile);
catch ME
    writelog(logfile, (ME.message))
    status = 1;
    return;
end

% Multi-step QSM processing including 1. PerformUnwrapping; 2.
% CreateBrainMask; 3. RemoveBackground; 4. CalculateQSM;
StartMultiProcess;

% Done
logtxt = 'QSM_script processing completed.';
disp(logtxt);
writelog(logfile, '\n');
writelog(logfile, logtxt);
status = 0;
