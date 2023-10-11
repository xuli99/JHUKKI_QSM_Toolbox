%% Authors: Jiri van Bergen and Xu Li
% Affiliation: Radiology @ JHU - Kirby Center
% Contact via xuli@mri.jhu.edu

% updated 2019-07-09, X.L.
% updated 2021-05-04, X.L. fixed a bug
% updated 2021-08-11, X.L. fixed mac issue
% updated 2023-10-11, X.L. moved parameter update code here from CreateFcns

% Open file
if ismac
    [FileName,PathName,FilterIndex] = uigetfile('*', 'Select files (PAR/REC, DICOM, MAT)');
else
    [FileName,PathName,FilterIndex] = uigetfile(handles.Params.fileTypes, 'Select files (PAR/REC, DICOM, MAT)');
end

% Reset panels
set([handles.PanelStart handles.PanelStep1 handles.PanelStep2 handles.PanelStep3 handles.PanelStep4 handles.PanelDatasets], 'HighlightColor', [0 0 0]); % Black

% Something?
if(FileName ~= 0)
    %% Update text and find extension
    set(handles.TextFileName, 'String', FileName);
    [~,FileBaseName,FileExt] = fileparts(FileName);
    
    % Switch to that folder
    cd(PathName);
   
    [GREMag, GREPhase, Params, handles] = readerwrapper(PathName, FileName, handles);
    
    %% Determine TE's 
    if(Params.nEchoes > 1)
        Params.echoStep = 1;    % default        
        if isfield(Params, 'B0')
            switch Params.B0
                case {7, '7'}
                    Params.echoNums = find(Params.TEs >= handles.TELowerLimit_7T & Params.TEs <= handles.TEUpperLimit_7T);
                case {3, '3'}
                    Params.echoNums = find(Params.TEs >= handles.TELowerLimit_3T & Params.TEs <= handles.TEUpperLimit_3T);
                case {1.5, '1.5'}
                    Params.echoNums = find(Params.TEs >= handles.TELowerLimit_1p5T & Params.TEs <= handles.TEUpperLimit_1p5T);
                case {11.7, '11.7'}
                    Params.echoNums = find(Params.TEs >= handles.TELowerLimit_11p7T & Params.TEs <= handles.TEUpperLimit_11p7T);
                otherwise
                    Params.echoNums = 1:length(Params.TEs);            
            end
        else   % use 3T as default
            Params.echoNums = find(Params.TEs >= handles.TELowerLimit_3T & Params.TEs <= handles.TEUpperLimit_3T);
        end
        
        set(handles.TextEchoInfo, 'String',  sprintf('Using echoes  %d - %d  (TE %0.3g - %0.3g ms)', min(Params.echoNums) , max(Params.echoNums), Params.TEs(min(Params.echoNums))*1000, Params.TEs(max(Params.echoNums))*1000));
    else
        Params.echoNums = 1;
        Params.echoStep = 0;    % default
        set(handles.TextEchoInfo, 'String',  sprintf('Using echo  1  (TE %0.3g ms)',Params.TEs(1)*1000));
    end
    
    % Fill-in initial Parameters if loaded from Constants otherwise default
    % (moved from CreateFcns)
    set(handles.CheckSaveData,'Value', handles.Params.saveOutput);
    set(handles.VarTemplateEcho, 'String', num2str(handles.Params.TemplateEcho))

    % Fill masking echo
    handles.Params.SaveEcho = Params.echoNums(1);
    set(handles.VarMaskEchoes,'String',sprintf('[%d]',Params.echoNums(1)));
    
    % Fill in FSL folders, masking and SHARP parameter
    set(handles.VarFSL, 'String', handles.Params.FSLFolder);
    set(handles.VarFSLThres,'String', num2str(handles.Params.FSLThreshold));
    set(handles.VarRadiusDisk,'String',num2str(handles.Params.ErodeRadius));
    set([handles.VarSHARPradius],'String',num2str(handles.Params.SHARPradius));

    % Find Method Strings Dict
    Params.UnwrappingMethodsDict    = get(handles.VarUnwrappingMethod, 'String');
    Params.BgRemovalMethodsDict     = get(handles.VarBgRemoval, 'String');
    Params.QSMSolverDict            = get(handles.VarQSMSolver, 'String');
    
    % Save to workspace
    handles.Params      = Params;
    handles.GREPhase    = single(GREPhase);
    handles.GREPhaseRaw = single(GREPhase);              % for backup raw phase
    handles.GREPhaseRaw1 = single(GREPhase(:,:,:,1,1));  % just for display using the first echo/dynamic, 3D
    handles.GREMag      = single(GREMag);
    handles.slice2disp  = floor(Params.sizeVol(3)/2);    

    % Show
    handles = LoadImage(hObject, handles, handles.GREPhaseRaw1, 'Phase (rad)');
    
    % Update info
    ShowImageInfo;
    
    % Enable slider
    set(handles.ImageSlider,'Enable','On')
    set(handles.ImageSlider,'Value',handles.slice2disp)
    set(handles.ImageSlider,'Min',1)
    set(handles.ImageSlider,'Max',handles.Params.sizeVol(3))
    
    % Show Panels
    set([handles.PanelStart handles.PanelStep1 handles.PanelStep2 handles.PanelStep3 handles.PanelStep4 handles.PanelDatasets],'Visible','On')
    % Show buttons
    set([handles.ButtonShowRawPhase handles.ButtonShowPhase handles.ButtonShowMask handles.ButtonShowBg handles.ButtonShowQSM],'Visible','On')
    
    % Update panel and disable open button
    set(handles.PanelFiles, 'HighlightColor', [0 0.5 0]);
    set(handles.ButtonSelectFile, 'Enable', 'Off')
    set(handles.ButtonShowRawPhase, 'Enable', 'On')
    set([handles.ButtonAddDataset handles.ButtonStartDatasets handles.ButtonLoadDataList], 'Enable', 'Off');
    
    % if num of echoes is smaller than 2, cannot fit R2*, cannot use
    % template unwrapping
    if Params.nEchoes < 2
        set([handles.checkboxR2star], 'Visible', 'Off')
        set([handles.Tag_TemplateEcho_text], 'Visible', 'off')
        set([handles.VarTemplateEcho], 'Visible', 'off')
        handles.Params.TemplateEcho = 0;
    end
    
    % for dynamic data, don't fit R2*
    if Params.nDynamics > 1
        set([handles.checkboxR2star], 'Visible', 'Off')
    end
    
    % Set table
    set(handles.TableDatasets, 'Data', { [PathName FileName], handles.textReadyLoad });
    handles.CurrentDataset = 1;
end

% Save
guidata(hObject, handles);