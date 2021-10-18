%% Author: Jiri van Bergen and Xu Li
% Affiliation: Radiology @ JHU - Kirby Center
% Contact via xuli@mri.jhu.edu

% Updated by Xu Li, 2019-06-20
% Updated by Xu Li, 2021-06-23, called in cluster version

%% Show Image info in GUI or print
% Format
Params.voxSize = single(Params.voxSize);

TextDimensions = sprintf('Dimensions:  %d x %d x %d px', Params.sizeVol(1), Params.sizeVol(2), Params.sizeVol(3));
TextVoxelSize = sprintf('Voxel size:  %0.3g x %0.3g x %0.3g mm', Params.voxSize(1), Params.voxSize(2), Params.voxSize(3));
TextEchoes = sprintf('Echoes:  %d', Params.nEchoes);
TextDynamics = sprintf('Dynamics:  %d', Params.nDynamics);

% How many echoes?
if(length(Params.TEs) == 1)
    TextTEs = sprintf('TE:  %0.3g ms', Params.TEs(1)*1000);
else
    TextTEs = sprintf('TE''s: %0.3g - %0.3g ms  (delta-TE: %0.3g ms)', Params.TEs(1)*1000, Params.TEs(end)*1000, (Params.TEs(2)-Params.TEs(1))*1000);
end

if isfield(handles.Params, 'cluster')
    if handles.Params.cluster == 1
        handles.TextDimensions.String = TextDimensions;
        handles.TextVoxelSize.String = TextVoxelSize;
        handles.TextEchoes.String = TextEchoes;
        handles.TextDynamics.String = TextDynamics;
        handles.TextTEs.String = TextTEs;

        % now print infomation
        writelog(handles.logfile, [handles.TextDimensions.String, '\n']);
        writelog(handles.logfile, [handles.TextVoxelSize.String, '\n']);
        writelog(handles.logfile, [handles.TextEchoes.String, '\n']);
        writelog(handles.logfile, [handles.TextDynamics.String, '\n']);
        writelog(handles.logfile, [handles.TextTEs.String, '\n']);
    end
else 
    % Update values in the information table, GUI
    set(handles.TextDimensions, 'String', TextDimensions);
    set(handles.TextVoxelSize, 'String', TextVoxelSize);
    set(handles.TextEchoes, 'String', TextEchoes)
    set(handles.TextDynamics, 'String', TextDynamics)
    set(handles.TextTEs, 'String', TextTEs)
end