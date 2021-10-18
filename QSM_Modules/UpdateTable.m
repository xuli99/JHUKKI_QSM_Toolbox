function handles = UpdateTable(handles, statustext )
    %% Author: Jiri van Bergen
    % Affiliation: Radiology @ JHU - Kirby Center
    % Contact via xuli@mri.jhu.edu
    % updated 2021-06-24, cluster version
    
    % Update table with new status for this dataset
    if isfield(handles.Params, 'cluster')
        if handles.Params.cluster == 1
            handles.TableDatasets.Data{handles.CurrentDataset, 2} = statustext;
        end
    else
        curTableData = get(handles.TableDatasets, 'Data');
        curTableData{handles.CurrentDataset, 2} = statustext;
        set(handles.TableDatasets, 'Data', curTableData);
    end
end

