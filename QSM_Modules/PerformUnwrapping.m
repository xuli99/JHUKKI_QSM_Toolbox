%% Author: Jiri van Bergen and Xu Li
% Affiliation: Radiology @ JHU - Kirby Center
% Contact via xuli@mri.jhu.edu

% Updated by Xu Li, 2019-07-09
% Updated by Xu Li, 2019-09-07, Added LapPhaseCorrection
% Updated by Xu Li, 2020-10-16, Updated with new phase_unwrap_laplacian
% Updated 2021-06-27 X.L., cluster version

%% Perform Phase Unwrapping
% Remove open waitbars!
wh=findall(0,'tag','TMWWaitbar');
delete(wh);

Params = handles.Params;

if ~isfield(handles.Params, 'cluster')  % GUI only
    % ------------------------------------------------------------
    % Disable buttons
    set([handles.ButtonEditEchoes handles.VarFSLThres handles.VarB0 handles.VarRadiusDisk handles.VarMaskEchoes handles.VarFSL handles.VarBgRemoval handles.VarSHARPradius handles.VarQSMSolver], 'Enable', 'Off');

    % Update panel
    set([handles.PanelStep1 handles.PanelStep2 handles.PanelStep3 handles.PanelStep4], 'HighlightColor', [0 0 0]); % Black
    set(handles.PanelStart, 'HighlightColor', [0 0.5 0]); % Green
end

% Get variables
GREPhase    = handles.GREPhase;

% Output
StringApp1 = ['_PhaseUnwrapped_echo' num2str(Params.echoNums(1))]; 

if(length(Params.echoNums) > 1)
    if Params.echoStep > 1
        StringApp2  = ['-' num2str(Params.echoNums(end)), '_s', num2str(Params.echoStep)];        
    else
        StringApp2  = ['-' num2str(Params.echoNums(end))];
    end
else
    StringApp2  = '';
end

switch Params.UnwrappingMethodsDict{Params.UnwrappingMethod}
    case 'Path'
        StringApp3 = 'Path';        
    case 'Laplacian'
        StringApp3 = 'Lap';
    case 'NonlinearFit + Path'
        StringApp3 = 'NLFpath';
    otherwise
        error('Unknown unwrapping method.')
end

if Params.nDynamics == 1
    stringApp4 = '';    % default with just 1 dynamic
else
    stringApp4 = ['_dyn', num2str(Params.nDynamics)]; % 5D
end

outputFile = [Params.FileBaseName, StringApp1, StringApp2, StringApp3, stringApp4];

% Switch
cd(Params.PathName);

% Check if exists
if(exist(fullfile(cd, [outputFile, '.mat']), 'file') == 2 && Params.saveOutput)
    textWaitbar = 'Loading previously unwrapped phase...';
    if ~isfield(handles.Params, 'cluster')  % GUI only
        multiWaitbar(textWaitbar, 'Busy');
    else
        disp(textWaitbar)
        writelog(handles.logfile, [textWaitbar, '...']);
    end
    load([outputFile, '.mat'])      
    if ~isfield(handles.Params, 'cluster')  % GUI only 
        multiWaitbar('CloseAll');
    else
        disp('Done.')
        writelog(handles.logfile, 'Done. \n');
    end
else
    textWaitbar = ['Unwrapping ' num2str(length(Params.TEs)) ' echoes, ', num2str(Params.nDynamics), ' dynamics'];
    if ~isfield(handles.Params, 'cluster')  % GUI only 
        multiWaitbar(textWaitbar, 0, 'CanCancel', 'On');
    else
        disp(textWaitbar)
        writelog(handles.logfile, [textWaitbar, '...']);
    end
    
    switch Params.UnwrappingMethodsDict{Params.UnwrappingMethod} 
        case 'Laplacian'
            for dynamic_ind = 1:Params.nDynamics
                for echo_ind = 1:length(Params.TEs)
                    % Unwrapping all the echoes
                    if Params.phase2DprocFlag == 0  
                        % default is 3D
                        GREPhase(:,:,:,echo_ind,dynamic_ind) = phase_unwrap_laplacian(GREPhase(:,:,:,echo_ind,dynamic_ind), Params, [], 2);
                    else
%                         % 2D phase unwrapping slice by slice
%                         textWaitbar2D = 'Performing 2D phase unwrapping first';
%                         multiWaitbar(textWaitbar2D, 0, 'Color', 'b' );
%                         
%                         % 2D unwrapping not good
%                         for sliceii = 1:Params.sizeVol(3)
%                             GREPhase(:,:,sliceii,echo_ind,dynamic_ind) = phase_unwrap_laplacian_2D(GREPhase(:,:,sliceii,echo_ind,dynamic_ind), Params, 0, 2);  % no refVox
%                             multiWaitbar(textWaitbar2D, sliceii/Params.sizeVol(3), 'Color', 'b' );
%                         end            
                        GREPhase(:,:,:,echo_ind,dynamic_ind) = phase_unwrap_laplacian(GREPhase(:,:,:,echo_ind,dynamic_ind), Params, [], 2);
                    end
                    if ~isfield(handles.Params, 'cluster')  % GUI only 
                        % Waitbar        
                        hasCanceled = multiWaitbar( textWaitbar, (echo_ind/length(Params.TEs))*(dynamic_ind/Params.nDynamics));
                        HandleStopReconstruction;
                    else
                        disp([num2str(100*(echo_ind/length(Params.TEs))*(dynamic_ind/Params.nDynamics)), '% Done.']);
                    end
                end
            end                

        case 'NonlinearFit + Path'
            if length(Params.echoNums) > 2  % use nonlinear formula  
                GREMag = handles.GREMag;    % complex fitting need GREMag too   
                for dynamic_ind = 1:Params.nDynamics                     
                    [GREPhase(:,:,:,1,dynamic_ind), N_std, relres, ~] = Fit_ppm_complex(GREMag(:,:,:,Params.echoNums,dynamic_ind).*cos(GREPhase(:,:,:,Params.echoNums,dynamic_ind)) ...
                                                        + 1i*GREMag(:,:,:,Params.echoNums,dynamic_ind).*sin(GREPhase(:,:,:,Params.echoNums,dynamic_ind)));                                            
                    if ~isfield(handles.Params, 'cluster')  % GUI only
                        multiWaitbar( textWaitbar, 0.6*(dynamic_ind/Params.nDynamics)); 
                    end
                    GREPhase(:,:,:,1,dynamic_ind) = phase_unwrap_path_mex(double(GREPhase(:,:,:,1,dynamic_ind)));    % path based mex version, still in radian
                    if ~isfield(handles.Params, 'cluster')  % GUI only
                        hasCanceled = multiWaitbar(textWaitbar, dynamic_ind/Params.nDynamics);  
                    else
                        disp([num2str(100*dynamic_ind/Params.nDynamics), '% Done.']);
                    end
                    
                    DPWeight = 1./N_std;
                    DPWeight(isinf(DPWeight)) = 0;
                   
                    % May want to add SPUR as further unwrapping
                end
            else  % 1 or 2 echoes, otherwise just do pathbased unwrapping, change GREPhase to only one volume
                for dynamic_ind = 1:Params.nDynamics 
                    GREPhase(:,:,:,:,dynamic_ind) = double(GREPhase(:,:,:,:,dynamic_ind));
                    for echo_ind = 1:length(Params.echoNums)
                        GREPhase(:,:,:,Params.echoNums(echo_ind),dynamic_ind) = phase_unwrap_path_mex(GREPhase(:,:,:,Params.echoNums(echo_ind),dynamic_ind));
                        if ~isfield(handles.Params, 'cluster')  % GUI only
                            hasCanceled = multiWaitbar( textWaitbar, (echo_ind/length(Params.echoNums))*(dynamic_ind/Params.nDynamics) );
                        else
                            disp([num2str(100*(echo_ind/length(Params.echoNums))*(dynamic_ind/Params.nDynamics)), '% Done.']);
                        end
                    end
                    if length(Params.echoNums) > 1
                        GREPhase(:,:,:,1,dynamic_ind) = find_slopeintercept_phasevstime_fast(GREPhase(:,:,:,Params.echoNums,dynamic_ind), [0,1]);
                    else
                        GREPhase(:,:,:,1,dynamic_ind) = GREPhase(:,:,:,Params.echoNums,dynamic_ind);
                    end
                end
            end
            GREPhase = GREPhase(:,:,:,1,:); % echo dimention collapse to 1
                        
        case 'Path'
            GREPhase = double(GREPhase);
            N = size(GREPhase);
            RefVox = [floor(N(1)/2), floor(N(2)/2), floor(N(3)/2)]; 

            for dynamic_ind = 1:Params.nDynamics
                for echo_ind = 1:length(Params.TEs)
                    GREPhase(:,:,:,echo_ind,dynamic_ind) = phase_unwrap_path_mex(GREPhase(:,:,:,echo_ind,dynamic_ind));
                    if ~isfield(handles.Params, 'cluster')  % GUI only
                        hasCanceled = multiWaitbar( textWaitbar, (echo_ind/length(Params.TEs))*(dynamic_ind/Params.nDynamics));
                    else
                        disp([num2str(100*((echo_ind/length(Params.TEs))*(dynamic_ind/Params.nDynamics))), '% Done.']);
                    end
                    GREPhase(:,:,:,echo_ind,dynamic_ind) = GREPhase(:,:,:,echo_ind,dynamic_ind) -  GREPhase(RefVox(1), RefVox(2), RefVox(3), echo_ind,dynamic_ind);
                end
            end

        otherwise
             error('Unknown unwrapping method.')

    end
    
    % save data
    switch Params.UnwrappingMethodsDict{Params.UnwrappingMethod} 
        case {'Laplacian', 'Path'}
            save([outputFile '.mat'], 'GREPhase', 'Params');
        case {'NonlinearFit + Path'}
            save([outputFile '.mat'], 'GREPhase', 'Params', 'DPWeight');
        otherwise
            error('Unknown unwrapping method.')
    end
    
    if ~isfield(handles.Params, 'cluster')  % GUI only
        multiWaitbar('CloseAll');
    else
        disp('Done.')
        writelog(handles.logfile, 'Done. \n');
    end
end

% Save
switch Params.UnwrappingMethodsDict{Params.UnwrappingMethod} 
    case {'Laplacian', 'Path'}
        handles.GREPhase        = GREPhase;
    case {'NonlinearFit + Path'}
        handles.GREPhase        = GREPhase;
        handles.DPWeight        = DPWeight;
    otherwise
        error('Unknown unwrapping method.')
end

if ~isfield(handles.Params, 'cluster')  % GUI only
    guidata(hObject, handles);

    % Display
    handles = LoadImage(hObject, handles, handles.GREPhase, 'Phase (rad)');

    % Update panel
    set(handles.PanelStep1, 'HighlightColor', [0 0.5 0]);
    set(handles.ButtonShowPhase, 'Enable', 'On')
end

% Update table
handles = UpdateTable(handles, 'Completed 1 of 4');
