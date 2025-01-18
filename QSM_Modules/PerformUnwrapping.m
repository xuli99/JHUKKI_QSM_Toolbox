%% Author: Jiri van Bergen and Xu Li
% Affiliation: Radiology @ JHU - Kirby Center
% Contact via xuli@mri.jhu.edu

% Updated by Xu Li, 2019-07-09
% Updated by Xu Li, 2019-09-07, Added LapPhaseCorrection
% Updated by Xu Li, 2020-10-16, Updated with new phase_unwrap_laplacian
% Updated 2021-06-27 X.L., cluster version
% Updated 2023-04-04 X.L., added ROMEO option for cluster version
% Updated 2023-06-01 X.L., saving format update
% Updated 2023-10-09 X.L., adding template unwrapping for path-based method 
% Updated 2024-06-01 X.L., added phase_quality_map for path-based method 
% Updated 2024-06-24 X.L., added RefVox check for path-based method to increase robustness
% Updated 2025-01-18 X.L., bug fix

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
    case 'ROMEO'
        StringApp3 = 'ROMEO';
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
    
    if isfield(handles.Params, 'RefVox')
        RefVox = handles.Params.RefVox;
    else
        RefVox = [];
    end

    switch Params.UnwrappingMethodsDict{Params.UnwrappingMethod} 
        case 'Laplacian'
            for dynamic_ind = 1:Params.nDynamics
                for echo_ind = 1:length(Params.TEs)

%                     % Unwrapping all the echoes
%                     if Params.phase2DprocFlag == 0  
%                         % default is 3D
%                         GREPhase(:,:,:,echo_ind,dynamic_ind) = phase_unwrap_laplacian(GREPhase(:,:,:,echo_ind,dynamic_ind), Params, RefVox, 2);
%                     else
% %                         % 2D phase unwrapping slice by slice
% %                         textWaitbar2D = 'Performing 2D phase unwrapping first';
% %                         multiWaitbar(textWaitbar2D, 0, 'Color', 'b' );
% %                         
% %                         % 2D unwrapping not good
% %                         for sliceii = 1:Params.sizeVol(3)
% %                             GREPhase(:,:,sliceii,echo_ind,dynamic_ind) = phase_unwrap_laplacian_2D(GREPhase(:,:,sliceii,echo_ind,dynamic_ind), Params, 0, 2);  % no refVox
% %                             multiWaitbar(textWaitbar2D, sliceii/Params.sizeVol(3), 'Color', 'b' );
% %                         end            
%                         GREPhase(:,:,:,echo_ind,dynamic_ind) = phase_unwrap_laplacian(GREPhase(:,:,:,echo_ind,dynamic_ind), Params, RefVox, 2);
%                     end
                    
                    GREPhase(:,:,:,echo_ind,dynamic_ind) = phase_unwrap_laplacian(GREPhase(:,:,:,echo_ind,dynamic_ind), Params, RefVox, 2);

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
                    [GREPhase(:,:,:,1,dynamic_ind),~] = phase_unwrap_path_mex(double(GREPhase(:,:,:,1,dynamic_ind)));    % path based mex version, still in radian
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
                        [GREPhase(:,:,:,Params.echoNums(echo_ind),dynamic_ind), ~] = phase_unwrap_path_mex(GREPhase(:,:,:,Params.echoNums(echo_ind),dynamic_ind));
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
            
            if isempty(RefVox) 
                RefVox = [floor(N(1)/2), floor(N(2)/2), floor(N(3)/2)]; 
            end

            nEchoes = Params.nEchoes;
            
            % added phase_quality_map from path-based unwrapping (reciprocal of sum of 2nd diff)
            phase_quality_map = zeros([Params.sizeVol, Params.nDynamics]);

            for dynamic_ind = 1:Params.nDynamics
                GREPhase_Ref = zeros(nEchoes,1);    % do it for each dynamic separately
                % unwrap each echo
                for echo_ind = 1:nEchoes
                    [GREPhase(:,:,:,echo_ind,dynamic_ind), phase_quality_echo] = phase_unwrap_path_mex(GREPhase(:,:,:,echo_ind,dynamic_ind));
                    
                    if echo_ind == 1
                        phase_quality_map(:,:,:,dynamic_ind) = phase_quality_echo; 
                    else
                        phase_quality_map(:,:,:,dynamic_ind) = min(phase_quality_map(:,:,:,dynamic_ind), phase_quality_echo);
                    end

                    % display progress
                    if ~isfield(handles.Params, 'cluster')  % GUI only
                        hasCanceled = multiWaitbar( textWaitbar, (echo_ind/nEchoes)*(dynamic_ind/Params.nDynamics));
                    else
                        disp([num2str(100*((echo_ind/nEchoes)*(dynamic_ind/Params.nDynamics))), '% Done.']);
                    end

                    GREPhase_Ref(echo_ind) = GREPhase(RefVox(1), RefVox(2), RefVox(3), echo_ind, dynamic_ind);
                end
                
                if nEchoes > 1
                    % Adjust RefVox in cases of large jump (e.g. in veins) to avoid potential Error
                    % only for multi-echo data
                    while abs(GREPhase_Ref(2) - GREPhase_Ref(1)) > pi && all(RefVox(1:2) < floor(N(1:2)*0.6))
                        shiftstep = 3;  % in-slice shift
                        RefVox(1) = RefVox(1) + shiftstep;
                        RefVox(2) = RefVox(2) + shiftstep;
                        GREPhase_Ref = squeeze(GREPhase(RefVox(1), RefVox(2), RefVox(3), :, dynamic_ind));
                    end
    
                    if abs(GREPhase_Ref(2) - GREPhase_Ref(1)) > pi
                        % adjustment not successfuly
                        disp('RefVox may be problemetic.')
                    end
                end

                for echo_ind = 1:nEchoes
                    GREPhase(:,:,:,echo_ind,dynamic_ind) = GREPhase(:,:,:,echo_ind,dynamic_ind) -  GREPhase_Ref(echo_ind);
                end

                % using template unwrapping to increase robustness
                if Params.TemplateEcho > 0
                    disp('Doing template unwrapping.')
                    GREPhase(:,:,:,:,dynamic_ind) = phase_unwrap_template(GREPhase(:,:,:,:,dynamic_ind), Params.TEs, Params.TemplateEcho);
                end

                if isfield(Params, 'UnwrapCheck')
                    % unwrap echo with potential residual wraps, opt
                    disp('rerun unwrapping for double checking ...')
                    for echo_ind = 1:nEchoes
                        [GREPhase(:,:,:,echo_ind,dynamic_ind), ~] = phase_unwrap_path_mex(GREPhase(:,:,:,echo_ind,dynamic_ind));
                        GREPhase(:,:,:,echo_ind,dynamic_ind) = GREPhase(:,:,:,echo_ind,dynamic_ind) -  GREPhase(RefVox(1), RefVox(2), RefVox(3), echo_ind,dynamic_ind);
                    end
                end
            end

        case 'ROMEO'

            % need to setup ROMEO path and parameters in ParameteSetting files
            % make it availabe for cluster version only

            GREMag = handles.GREMag;    % ROMEO is better with GREMag, with template unwrapping, skip referencing
            romeo_parameters = handles.Params.romeo_parameters;

            % parameters set according to each data
            GREPhase_romeo_b0 = zeros([Params.sizeVol, Params.nDynamics]);
            phase_quality_map = zeros([Params.sizeVol, Params.nDynamics]);

            tempdir = './tmp';
            romeo_parameters.output_dir = fullfile(tempdir, 'romeo_tmp'); % if not set pwd() is used
            romeo_parameters.TE = Params.TEs*1e3;     % required for multi-echo, in the unit of ms
            romeo_parameters.voxel_size = Params.voxSize;
            mkdir(romeo_parameters.output_dir);
            
            for dynamic_ind = 1:Params.nDynamics

                romeo_parameters.mag = GREMag(:,:,:,:,dynamic_ind);

                [GREPhase(:,:,:,:,dynamic_ind), GREPhase_romeo_b0(:,:,:,dynamic_ind)] = ROMEO(GREPhase(:,:,:,:,dynamic_ind), romeo_parameters);
                quality_fn = fullfile(romeo_parameters.output_dir, 'quality.nii');
                
                phase_quality_map(:,:,:,dynamic_ind) = load_untouch_nii(quality_fn).img;

            end

            if romeo_parameters.cleanup == 1
                [status, msg] = rmdir(tempdir, 's'); % remove the temporary ROMEO output folder
            end

        otherwise
             error('Unknown unwrapping method.')

    end
    
    % save data
    switch Params.UnwrappingMethodsDict{Params.UnwrappingMethod} 
        case {'Laplacian'}
            save([outputFile '.mat'], 'GREPhase', 'Params', '-v7.3');
        case {'Path'}
            save([outputFile '.mat'], 'GREPhase', 'Params', 'phase_quality_map', '-v7.3');
        case {'NonlinearFit + Path'}
            save([outputFile '.mat'], 'GREPhase', 'Params', 'DPWeight', '-v7.3');
        case {'ROMEO'}
            save([outputFile '.mat'], 'GREPhase', 'Params', 'GREPhase_romeo_b0', 'phase_quality_map', '-v7.3');
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
handles.GREPhase        = GREPhase;
switch Params.UnwrappingMethodsDict{Params.UnwrappingMethod} 
    case {'Laplacian'}
        % do nothing
    case {'Path'}
        if exist("phase_quality_map", "var")
            % for compatibility
            handles.phase_quality_map = phase_quality_map;
        end
    case {'NonlinearFit + Path'}
        handles.DPWeight        = DPWeight;
    case {'ROMEO'}
        handles.phase_quality_map = phase_quality_map;
        handles.GREPhase_romeo_b0 = GREPhase_romeo_b0;
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
