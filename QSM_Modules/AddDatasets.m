%% Author: Jiri van Bergen and Xu Li
% Affiliation: Radiology @ JHU - Kirby Center
% Contact via xuli@mri.jhu.edu

% Open file
[FileName,PathName,FilterIndex] = uigetfile(handles.Params.fileTypes, 'Select files (PAR/REC, DICOM, NIFTI)');

% Something?
if(FileName ~= 0)
    % Get data
    tableData   = get(handles.TableDatasets, 'Data');
    
    % Check if in there
    % BECAUSE STRFIND SUCKS FOR ARRAY OF CELLS!
    isNew = true;
    for c = 1:size(tableData,1)
        % Already in there
        if(strcmpi(tableData{c,1},[PathName FileName]))
            errordlg('Dataset already listed');
            isNew = false;
            return;
        end
    end
    
    % New?
    if(isNew)
        % New row to add
        newRow      =  { [PathName FileName], handles.textReadyLoad };
        tableData   = [tableData; newRow];
        
        % Update
        set(handles.TableDatasets, 'Data', tableData);
    end
end