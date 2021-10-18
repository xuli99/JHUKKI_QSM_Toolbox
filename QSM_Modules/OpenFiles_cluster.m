function handles = OpenFiles_cluster(fullFileName, handles)

%% Authors: Jiri van Bergen and Xu Li
% Affiliation: Radiology @ JHU - Kirby Center
% Contact via xuli@mri.jhu.edu
% cluster version

if exist(fullFileName, 'file')
    handles.TableDatasets.Data = { fullFileName, handles.textReadyLoad };
    handles.CurrentDataset = 1;
end
