%% Author: Jiri van Bergen and Xu Li
% Affiliation: Radiology @ JHU - Kirby Center
% Contact via xuli@mri.jhu.edu

% Updated 2019-07-09 X.L.
% Updated 2021-06-16 X.L., improved memory efficiency
% Updated 2021-06-24 X.L., cluster version
% Updated 2024-06-07 X.L., for cluster array

%% Start processing of multiple datasets
% Save Constants file
Params          = handles.Params;
if ~isfield(Params,'NoSaveQSMsetting')
    % default
    save(handles.Params.QSMSettingsFile, 'Params');
else
    if Params.NoSaveQSMsetting == 0
        save(handles.Params.QSMSettingsFile, 'Params');
    end
end

% Get table data
if isfield(handles.Params, 'cluster')
    if handles.Params.cluster == 1
        tableData = handles.TableDatasets.Data;
    end
else
    tableData = get(handles.TableDatasets, 'Data');
    % Disable buttons
    set([handles.ButtonEditEchoes handles.VarB0 handles.VarRadiusDisk handles.VarMaskEchoes handles.VarFSL handles.VarBgRemoval handles.VarSHARPradius handles.VarQSMSolver], 'Enable', 'Off');
end

% Go through all the data
for c = 1:size(tableData, 1)
    if(strcmpi(tableData{c,2}, handles.textReadyLoad) || strcmpi(tableData{c,2}, handles.textReconstructing))
        % Update number and text
        handles.CurrentDataset = c;
        handles = UpdateTable(handles, 'Reconstructing...');
        
        % Filename
        FileFullName = tableData{c,1};
        
        % Extract info
        [PathName, FileBaseName, FileExt] = fileparts(FileFullName);
        
        % Switch to that folder
        cd(PathName);
        
        % readin data
        [GREMag, GREPhase, Params, handles] = readerwrapper([PathName, filesep], [FileBaseName, FileExt], handles);
        
        % cluster only
        if isfield(handles.Params, 'cluster')
            writelog(handles.logfile, '\n');
            writelog(handles.logfile, ['Data ', num2str(c), ': ', strrep(handles.TableDatasets.Data{c,1}, '\', '\\'), ' \n']);
            if handles.Params.cluster == 1
                % if num of echoes is smaller than 2, cannot fit R2*
                if Params.nEchoes < 2
                    handles.Params.R2starFlag           = 0;
                end

                % for dynamic data, don't fit R2*
                if Params.nDynamics > 1
                    handles.Params.R2starFlag           = 0;
                end

                % check
                if Params.echoNums(end) > length(Params.TEs)
                    writelog(handles.logfile, ['Echoes available: ', num2str(length(Params.TEs)), '\n']);
                    error('Need to reset selected echoes in Params.echoNums.')
                end
            end
        end        
        
        % Save to workspace
        handles.GREPhase    = single(GREPhase);       
        handles.GREPhaseRaw = single(GREPhase);              % for backup raw phase
        handles.GREMag      = single(GREMag);
        handles.Params      = Params;        

        if ~isfield(handles.Params, 'cluster')  % GUI only
            handles.GREPhaseRaw1 = single(GREPhase(:,:,:,1,1));    % just for display using the first echo/dynamic, 3D
            handles.slice2disp  = floor(Params.sizeVol(3)/2);
            handles = LoadImage(hObject, handles, handles.GREPhaseRaw1, 'Phase (rad)');
            
            % Enable slider
            set(handles.ImageSlider,'Enable','On')
            set(handles.ImageSlider,'Value',handles.slice2disp)
            set(handles.ImageSlider,'Min',1)
            set(handles.ImageSlider,'Max',handles.Params.sizeVol(3))
        end

        % clear redundant
        clear GREMag GREPhase
        
        % Update info
        ShowImageInfo;
        
        % free some memory 
        handles.chi_res = [];
        handles.freqMap = [];
        handles.maskErode = [];
        if isfield(handles, 'DPWeight')
            handles.DPWeight = [];
        end
        
        %% LETS DO THISSSSS
        PerformUnwrapping;
        CreateBrainMask;
        RemoveBackground;
        CalculateQSM;
    end
    
end

