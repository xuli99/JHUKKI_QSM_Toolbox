%% Author: Jiri van Bergen and Xu Li
% Affiliation: Radiology @ JHU - Kirby Center
% Contact via xuli@mri.jhu.edu

% Updated by Xu Li, 2019-07-09
% Updated by Xu Li, 2021-02-09, added phase offset removal for Echo Averaging
% Updated by Xu Li, 2021-05-25, optimized weighting for Echo Averaging
% Updated 2021-06-28, X.L., cluster version
% Updated 2021-08-20, X.L., bug fix for dyanmic data

%% Get variables
Params      = handles.Params;
GREPhase    = handles.GREPhase;
GREMag      = handles.GREMag;
maskErode   = handles.maskErode;
GREPhaseRaw = handles.GREPhaseRaw;

if Params.single_step_flag   
    % skip assigment of background removal naming
    outputFile = [Params.FileBaseName, '_tmp'];
else
    % if multi-step QSM, do background removal
    if strcmp(Params.UnwrappingMethodsDict{Params.UnwrappingMethod}, 'NonlinearFit + Path')
        StringApp_Avg = 'NLSlope';    
    elseif Params.EchoAvg > 0
        StringApp_Avg = 'avg';
    else
        StringApp_Avg = 'LSlope';
    end

    %% determined output file names
    switch Params.BgRemovalMethodsDict{Params.BgRemoval}
        case 'VSHARP'
            %% Use V-SHARP algorithm, Average afterward
            StringApp0 = 'SMV';
            StringApp1 = ['_', StringApp0, '_', num2str(Params.SHARPradius), '_', StringApp_Avg];
            WaitBarMsgloading = ['V-SHARP ', StringApp_Avg];

        case 'PDF'
            %% DIPOLE FITTING
            StringApp0 = 'PDF';
            StringApp1 = ['_', StringApp0, '_', StringApp_Avg];
            WaitBarMsgloading = ['PDF ', StringApp_Avg];

        case 'LBV+VSHARP'
            %% using LBV+VSHARP
            StringApp0 = 'LBVSMV';
            StringApp1 = ['_', StringApp0, '_', num2str(Params.SHARPradius), '_', StringApp_Avg]; 
            WaitBarMsgloading = ['LBV+SHARP ', StringApp_Avg];

        case 'iRSHARP'
            %% iRSHARP 
            StringApp0 = 'iRSHARP';
            StringApp1 = ['_', StringApp0, '_', num2str(Params.SHARPradius), '_', StringApp_Avg]; 
            WaitBarMsgloading = ['iRSHARP ', StringApp_Avg];

        otherwise
            error('unknown option.')
    end

    %% naming conventions
    if (length(Params.echoNums) > 1)            
        StringApp2 = '';
        StringApp4 = ['-', num2str(Params.echoNums(end))];
    else
        StringApp2 = '';
        StringApp4 = '';
    end

    StringApp3 = ['_echo', num2str(Params.echoNums(1))];

    if Params.echoStep > 1      % For Single Echo, Params.echoStep=0;
        StringApp5 = ['_s', num2str(Params.echoStep)];
    else
        StringApp5 = '';
    end

    outputFile = [Params.FileBaseName, StringApp1, StringApp2, StringApp3, StringApp4, StringApp5];
end

% normalized phase by 2*pi and TEs
if strcmp(Params.UnwrappingMethodsDict{Params.UnwrappingMethod}, 'NonlinearFit + Path')
    % normalized slope
    if length(Params.echoNums) > 1
        GREPhase = GREPhase./(2*pi*(Params.TEs(Params.echoNums(2)) - Params.TEs(Params.echoNums(1))));  % in Hz
    else  
        GREPhase = GREPhase./(2*pi*(Params.TEs(Params.echoNums)));  % in Hz
    end
    Params.echoNums = 1;   

elseif (Params.EchoAvg > 0)
    % check if need to recalculate
    if Params.single_step_flag || (exist(fullfile(cd, [outputFile, '.mat']), 'file') == 0 || Params.MaskWasMade) % if not saved or new mask
        if length(Params.echoNums) > 1
            % remove phase offset phi0, if with multiple echoes
            disp('removing phase offset for echo averaging ...')
            PhaseOffset = zeros([size(maskErode), 1, Params.nDynamics]);

            for dynamic_ind = 1:Params.nDynamics    
                [~,temp] = find_slopeintercept_phasevstime_fast(GREPhase(:,:,:,:,dynamic_ind), Params.TEs);            
                disp('polyfit offset field ...')
                PhaseOffset(:,:,:,1,dynamic_ind) = polyfitn3d(temp, maskErode, 5, 1);
            end
            GREPhase = GREPhase - PhaseOffset;
            disp('Done.')
        end

        % do weighted averaging
        GREPhase = GREPhase.*GREMag/(2*pi);
        denorm = sum(GREMag(:,:,:,Params.echoNums,:).*reshape(Params.TEs(Params.echoNums), [1,1,1,length(Params.echoNums)]), 4);
    %     denorm = zeros(size(GREMag(:,:,:,1,:)));
    %     for selectedEcho = Params.echoNums %1:length(Params.TEs)
    %         TEsSE = Params.TEs(selectedEcho);
    %         % GREPhase(:,:,:,selectedEcho,:) = GREPhase(:,:,:,selectedEcho,:)./(2*pi*TEsSE);  % ME in Hz
    %         denorm = denorm + TEsSE*GREMag(:,:,:,selectedEcho,:);
    %     end 
        GREPhase = GREPhase./(denorm + eps)*length(Params.echoNums);
    end    
else
    % if multi-echo, fit the linear slope
    normTEs = 0:(length(Params.echoNums)-1);
    for dynamic_ind = 1:Params.nDynamics
        GREPhase(:,:,:,1,dynamic_ind) = find_slopeintercept_phasevstime_fast(GREPhase(:,:,:,Params.echoNums,dynamic_ind), normTEs);
    end

    GREPhase = GREPhase(:,:,:,1,:); 
    GREPhase = GREPhase./(2*pi*(Params.TEs(Params.echoNums(2)) - Params.TEs(Params.echoNums(1))));
    Params.echoNums = 1;     
end

% Set selected echo for weighting
selectedEcho = Params.echoNums(1);

if Params.single_step_flag
    % total field in Hz
    if Params.EchoAvg > 0
        freqMap = mean(GREPhase, 4);
        Params.echoNums = 1; 
    else
        freqMap = GREPhase;     
    end
else
    %% checking file existance
    if(exist(fullfile(cd, [outputFile, '.mat']), 'file') == 2 && ...
            Params.saveOutput && ...
            ~Params.MaskWasMade) % Re-DO background removal for new mask
        textWaitbar = ['Loading previously ', WaitBarMsgloading, ' calculated field...'];
        if ~isfield(handles.Params, 'cluster')  % GUI only
            multiWaitbar( textWaitbar, 'Busy');
        else
            disp(textWaitbar)
            writelog(handles.logfile, [textWaitbar, '...']);
        end
        load([outputFile, '.mat']);
        if ~isfield(handles.Params, 'cluster')  % GUI only
            multiWaitbar('CloseAll');
        else
            disp('Done.')
            writelog(handles.logfile, 'Done. \n');
        end
        % Mark that we loaded and no new QSM calculation is required
        handles.Params.BgWasRemoved = 0;
    else
        %% Removing background gradient
        textWaitbar = ['Removing Background field for '  num2str(Params.nDynamics), ' dynamics.'];
        if ~isfield(handles.Params, 'cluster')  % GUI only
            multiWaitbar(textWaitbar, 0, 'CanCancel', 'On');
        else
            disp(textWaitbar)
            writelog(handles.logfile, [textWaitbar, '...']);
        end
        freqMap = zeros([Params.sizeVol, 1, Params.nDynamics]); 

        switch Params.BgRemovalMethodsDict{Params.BgRemoval}
            case 'VSHARP'

                %% Use V-SHARP algorithm
                clear GREMag
                for dynamic_ind = 1:Params.nDynamics
                    if prod(Params.sizeVol)*floor(Params.SHARPradius./min(Params.voxSize)) > 1e9  % in case ultra-high res
                        % V-SHARP MULTI - image space version
                        disp('VSHARP image space version.')
                        [freqMap(:,:,:,1,dynamic_ind), dpfield_fit] = SHARP_Adaptive_Multi(GREPhase(:,:,:,:,dynamic_ind), maskErode, Params.SHARPradius, Params.thresh_tsvd, Params, handles);
                    else            
                        % VSHARP k-space verion
                        radiusStep = min(Params.voxSize);
                        radiusArray = radiusStep:radiusStep:Params.SHARPradius;
                        [freqMap(:,:,:,1,dynamic_ind), dpfield_fit, mask_eval] = VSHARP_k(GREPhase(:,:,:,:,dynamic_ind), maskErode, radiusArray, Params.thresh_tsvd, Params, handles);
                    end

                    if Params.phase2DprocFlag == 1                      
                        if sum(mod(Params.sizeVol, 2)) > 0 
                            textWaitbar2D = 'Performing 2D V-SHARP';
                            if ~isfield(handles.Params, 'cluster')  % GUI only
                                multiWaitbar(textWaitbar2D, 0, 'Color', 'b' );
                            else
                                disp(textWaitbar2D);
                            end

                            for sliceii = 1:Params.sizeVol(3)
                                [freqMap(:,:,sliceii,1,dynamic_ind), dpfield_fit2D] = SHARP_Adaptive_Multi2D(freqMap(:,:,sliceii,:,dynamic_ind), maskErode(:,:,sliceii,:), ...
                                                        Params.SHARPradius, Params.thresh_tsvd, Params);
                                if ~isfield(handles.Params, 'cluster')  % GUI only
                                    multiWaitbar(textWaitbar2D, sliceii/Params.sizeVol(3), 'Color', 'b' );
                                end
                            end
                        else
                            [freqMap(:,:,:,1,dynamic_ind), ~, mask_eval] = VSHARP2D_k(freqMap(:,:,:,:,dynamic_ind), maskErode, radiusArray, Params.thresh_tsvd, Params, handles);                       
                        end
                    end    
                    if ~isfield(handles.Params, 'cluster')  % GUI only
                        hasCanceled = multiWaitbar(textWaitbar, (dynamic_ind/Params.nDynamics));
                        HandleStopReconstruction;
                    end
                end

                if exist('mask_eval', 'var')
                    maskErode = mask_eval;  
                end

            case 'PDF'
                %% DIPOLE FITTING       
                for dynamic_ind = 1:Params.nDynamics
                    if isfield(handles, 'DPWeight')
                        DPWeight = handles.DPWeight;
                    else
                        DPWeight = single(GREMag(:,:,:,selectedEcho,dynamic_ind));
                    end

                    if Params.EchoAvg > 0
                        GREPhase(:,:,:,1,dynamic_ind) = mean(GREPhase(:,:,:,Params.echoNums,dynamic_ind), 4);
                     end

                    [freqMap(:,:,:,1,dynamic_ind), dpfield_fit] = dipolefit_v2(GREPhase(:,:,:,1,dynamic_ind), maskErode, DPWeight, Params);

                    % LapPhaseCorrection, 2019-09-07, xl
                    if Params.LapPhaseCorrection == 1
                        temp = angle(exp(1i*freqMap(:,:,:,1,dynamic_ind).*(2*pi*Params.TEs(end))));
                        temp = phase_unwrap_laplacian(temp, Params, 0);     % no Ref
                        freqMap(:,:,:,1,dynamic_ind) = temp./(2*pi*Params.TEs(end)).*maskErode;
                    end    

                    if ~isfield(handles.Params, 'cluster')  % GUI only
                        hasCanceled = multiWaitbar(textWaitbar, (dynamic_ind/Params.nDynamics));
                        HandleStopReconstruction;
                    end
                end

            case 'LBV+VSHARP'
                clear GREMag

                for dynamic_ind = 1:Params.nDynamics

                    Params.lbv.tol = 0.001;     % default 0.1
                    Params.lbv.depth = -1;      % default
                    Params.lbv.peel = 0;        % similar to mask erosion

                    % do LBV echo by echo
                    for echoii = 1:length(Params.echoNums)
                        GREPhase(:,:,:,Params.echoNums(echoii),dynamic_ind) = LBV(GREPhase(:,:,:,Params.echoNums(echoii),dynamic_ind), maskErode, Params.sizeVol, Params.voxSize, ...
                                                                Params.lbv.tol, Params.lbv.depth, Params.lbv.peel);
                    end

                    % VSHARP k-space verion
                    radiusStep = min(Params.voxSize);
                    radiusArray = radiusStep:radiusStep:Params.SHARPradius;
                    [freqMap(:,:,:,1,dynamic_ind), dpfield_fit, mask_eval] = VSHARP_k(GREPhase(:,:,:,:,dynamic_ind), maskErode, radiusArray, Params.thresh_tsvd, Params, handles); 

                    if ~isfield(handles.Params, 'cluster')  % GUI only
                        hasCanceled = multiWaitbar(textWaitbar, (dynamic_ind/Params.nDynamics));
                        HandleStopReconstruction;
                    end
                end
                maskErode = mask_eval;

            case 'iRSHARP'
                clear GREMag
                for dynamic_ind = 1:Params.nDynamics
                    Params.iRSHARP_C = 0.25;    % add in Params to pass to iRHSARPv1 
                    [freqMap(:,:,:,1,dynamic_ind), dpfield_fit] = iRSHARPv1(GREPhase(:,:,:,:,dynamic_ind), GREPhaseRaw(:,:,:,:,dynamic_ind), maskErode, Params, handles);
                    if ~isfield(handles.Params, 'cluster')  % GUI only
                        hasCanceled = multiWaitbar(textWaitbar, (dynamic_ind/Params.nDynamics));
                        HandleStopReconstruction;
                    end
                end

            otherwise
                error('unknown option.')
        end

        % Save MAT file
        save([outputFile '.mat'], 'freqMap', 'Params', 'maskErode');     

        % Save NIFTI file
        saveNII(squeeze(freqMap), outputFile, Params, 1);    

        % Ready
        if ~isfield(handles.Params, 'cluster')  % GUI only
            multiWaitbar('CloseAll');            
        else
            disp('Done.')
            writelog(handles.logfile, 'Done. \n');
        end
        % Mark that we made a new calculation and so we need a new QSM
        handles.Params.BgWasRemoved = 1;

    end
end    

%% Save and update display
handles.freqMap        = freqMap;
handles.maskErode      = maskErode;

if ~isfield(handles.Params, 'cluster')  % GUI only
    guidata(hObject, handles);
    % Display
    handles = LoadImage(hObject, handles, handles.freqMap, 'Frequency (Hz)');
    
    % Update panel
    set(handles.PanelStep3, 'HighlightColor', [0 0.5 0]);
    set(handles.ButtonShowBg, 'Enable', 'On')
end

% Update table
UpdateTable(handles, 'Completed 3 of 4');