%% Author: Jiri van Bergen and Xu Li
% Affiliation: Radiology @ JHU - Kirby Center
% Contact via xuli@mri.jhu.edu

% Updated by Xu Li, 2019-07-09
% Updated by Xu Li, 2019-12-19
% Updated 2021-06-28 X.L., cluster version
% Updated 2021-07-08 X.L.
% Updated 2021-08-05 X.L., added unreliable phase mask2
% Updated 2022-03-22, X.L. added option to use RAS NIFTI
% Updated 2022-05-03, X.L., added switch for compatibility with some old mask options
% Updated 2023-04-04 X.L., added ROMEO option for cluster version
% Updated 2023-06-12, X.L., added wsl support for pc
% Updated 2024-06-01, X.L., added phase_quality_map option for path-based method

%% Get variables
Params      = handles.Params;
GREPhase    = handles.GREPhase;
GREMag      = handles.GREMag;
GREPhaseRaw = handles.GREPhaseRaw;

if strcmpi(Params.QSMSolverDict{Params.QSMSolver}, 'TFI')
    single_step_flag = true;
    handles.Params.single_step_flag = single_step_flag;
    outputFile  = [Params.FileBaseName '_mask_whole_head'];
else
    single_step_flag = false;
    handles.Params.single_step_flag = single_step_flag;
    if(length(Params.SaveEcho) > 1)
        % if using multiple echoes for masking
        outputFile  = [Params.FileBaseName '_brain_mask_echo' num2str(Params.SaveEcho(1)) '-' num2str(Params.SaveEcho(2)) '_r' num2str(Params.ErodeRadius) '_t' num2str(Params.FSLThreshold)];
    else
        outputFile  = [Params.FileBaseName '_brain_mask_echo' num2str(Params.SaveEcho) '_r' num2str(Params.ErodeRadius) '_t' num2str(Params.FSLThreshold)];
    end
end

%% Check if exists
if(exist(fullfile(cd, [outputFile, '.mat']), 'file') == 2 && Params.saveOutput)
    textWaitbar = 'Loading previously created mask...';
    if ~isfield(handles.Params, 'cluster')  % GUI only
        multiWaitbar(textWaitbar, 'Busy');
    else
        disp(textWaitbar)
        writelog(handles.logfile, [textWaitbar, '...']);
    end
    load([outputFile, '.mat'])
    if ~isfield(handles.Params, 'cluster')  % GUI only
        multiWaitbar( 'CloseAll' );
    else
        disp('Done.')
        writelog(handles.logfile, 'Done. \n');
    end
    % Mark that we loaded mask
    handles.Params.MaskWasMade = 0;
else
    if ~single_step_flag
        %% DOES FSL BET EXISITS?
        if ispc
            [status, cmdoutput] = system(['wsl ', handles.Params.FSLFolder, 'flirt -version']);
        else
            cmdoutput=[];
        end
        
        if(exist(Params.FSLFolder, 'dir') == 7) || (ispc && ~isempty(cmdoutput))
            %% New mask
            textWaitbar = 'Creating brain mask...';
            if ~isfield(handles.Params, 'cluster')  % GUI only
                multiWaitbar( textWaitbar, 0 );
            else
                disp(textWaitbar);
                writelog(handles.logfile, [textWaitbar, '...']);
            end

            % First save the magnitude data
            for i = 1:length(Params.SaveEcho)
                fileName = [ Params.FileBaseName, '_GREMag', num2str(Params.SaveEcho(i))];

                if ~isfield(handles.Params, 'FSLBETskip') % default, skip for non_brain data
                    % Remove previous files of FSL output OTHERWSIE FSL WILL NOT RUN
                    delete([fileName '_brain.*']);
                    delete([fileName '_brain_mask.*']);
                end
                
                % Save new ones            
                saveNII(GREMag(:,:,:,Params.SaveEcho(i), 1).*1, fileName, Params, 1);  
                
            end

            %% Using BET to extract BrainMask
            fname1 = [ Params.FileBaseName '_GREMag', num2str(Params.SaveEcho(1))];

            if ~isfield(handles.Params, 'FSLBETskip') % default, skip for non_brain data
                inputstring1 = [Params.FSLFolder, 'bet2 ', fname1, '.nii ' fname1, '_brain', ' -f ', Params.FSLThreshold, ' -g 0 -m' ];

                % added feature to run fsl using wsl on pc
                if ispc
                    [~, wslpwd] = system('wsl pwd');  % with \n at end
                    fname1_wsl = [wslpwd(1:end-1), '/', fname1];
                    wslsetenv = 'export FSLOUTPUTTYPE=NIFTI_GZ; ';
                    inputstring1 = ['wsl ', wslsetenv, Params.FSLFolder, 'bet2 ', fname1_wsl, ' ', fname1_wsl, '_brain', ' -f ', Params.FSLThreshold, ' -g 0 -m' ];
                end

                % Additional optional BET Flag
                if isfield(handles.Params, 'FSLBETAdditionalFlag')
                    inputstring1 = [inputstring1, ' ', handles.Params.FSLBETAdditionalFlag];
                end
                
                if ~isfield(handles.Params, 'cluster')  % GUI only
                    multiWaitbar( textWaitbar, 0.2 );
                end

                system(inputstring1);
            else
                disp('skipping FSL BET, assuming mask already created.')
            end

            % Save
            output = [fname1, '_brain_mask.nii.gz'];

            % Did it do anything?
            if(exist(fullfile(cd, output), 'file') ~= 2)
                errordlg('FSL BET did not run and no previous mask was found.. Reconstruction stopped');
                if ~isfield(handles.Params, 'cluster')  % GUI only
                    multiWaitbar( 'CloseAll' );
                end
                return;
            end

            nii = load_untouch_nii(output); 
            % for RAS NIFTI
            nii = nii_img_update(nii, Params);

            maskBET1 = permute(nii.img, [2,1,3]);       
            maskBET1 = maskBET1 > 0;
            maskErode = maskBET1;

            if length(Params.SaveEcho) > 1
                % Second mask at later echo, needs to be more conservative
                % using BET
                fname2 = [ Params.FileBaseName '_GREMag', num2str(Params.SaveEcho(2))];

                if ~isfield(handles.Params, 'FSLBETskip') % default, otherwise skip for non_brain data
                    inputstring1 = [Params.FSLFolder, 'bet2 ', fname2, '.nii ' fname2, '_brain', ' -f ', Params.FSLThreshold, ' -g 0 -m' ];

                    % added feature to run fsl using wsl on pc
                    if ispc
                        fname2_wsl = [wslpwd(1:end-1), '/', fname2];
                        inputstring1 = ['wsl ', wslsetenv, Params.FSLFolder, 'bet2 ', fname2_wsl, ' ', fname2_wsl, '_brain', ' -f ', Params.FSLThreshold, ' -g 0 -m' ];
                    end                

                    if ~isfield(handles.Params, 'cluster') 
                        multiWaitbar( textWaitbar, 0.5 );
                    end                    
                    system(inputstring1);
                end

                % Combine BET masks
                if ~isfield(handles.Params, 'cluster')  % GUI only
                    multiWaitbar( textWaitbar, 0.7 );
                end
                output = [fname2, '_brain_mask.nii.gz'];

                nii = load_untouch_nii(output);
                % for RAS NIFTI
                nii = nii_img_update(nii, Params);

                maskBET2 = permute(nii.img, [2,1,3]);
                maskBET2 = maskBET2 > 0;

                % V2
                zstart = 1;
                zend = 0.8*size(maskErode, 3);    
                ystart = 1; 
                yend1 = floor(size(maskErode, 1)*0.8); 

                xstart = floor(size(maskErode,2)*0.1);
                xend = floor(size(maskErode,2)*0.9);

                maskErode(ystart:yend1, xstart:xend,zstart:zend) = ...
                   maskErode(ystart:yend1, xstart:xend, zstart:zend) & maskBET2(ystart:yend1, xstart:xend, zstart:zend);
            end
            
%             % Hacking mask with Animal data
%             fname2 = [ Params.FileBaseName '_GREMag', num2str(Params.SaveEcho(1)), '_Manual'];
%             nii = load_untouch_nii([fname2, '.nii.gz']); 
%             nii = nii_img_update(nii, Params)
%             % may need to check RAS NIFTI
%             maskErode = permute(nii.img, [2,1,3]); maskErode=maskErode>0;            
                        
            % MaskOut Unreliable Phase
            if Params.nEchoes > 1    
                temp = GREPhaseRaw(:,:,:,:,1);  
                mask_intrinsic = (sum(abs(temp - min(temp(:))) < 10*eps, 4) >= (Params.nEchoes-1)); % air/bone mask
            else
                mask_intrinsic = zeros(size(maskErode));
            end

            disp('estimating unreliable phase ...')
            TV = TVOP;

            switch Params.UnwrappingMethodsDict{Params.UnwrappingMethod}
                case {'Path', 'NonlinearFit + Path', 'ROMEO'} 

                    mask_unrelyPhase = zeros(size(maskErode));

                    if strcmp(Params.UnwrappingMethodsDict{Params.UnwrappingMethod}, 'NonlinearFit + Path')
                        echoNumPick = 1;  
                    else                                    
                        echoNumPick = Params.echoNums;
                    end

                    if ~isfield(handles.Params, 'cluster')
                        UnrelThreshold = 2;                             % for residual wraps, 0-2
                    else
                        UnrelThreshold = handles.Params.UnreliThreshold;
                    end
                    
                    for echoNumPickii = 1:length(echoNumPick)                          
                        tempGrad = TV*(GREPhase(:,:,:,echoNumPick(echoNumPickii),1).*maskErode);

                        mask_unrelyPhase_temp = sum(abs(tempGrad) >= pi, 4) > UnrelThreshold;   
                        mask_unrelyPhase_temp = imdilate(mask_unrelyPhase_temp, strel('disk', 1));
                        tempMask = ~mask_unrelyPhase_temp.*maskErode;              

                        CC = bwconncomp(tempMask, 6);       
                        numPixels = cellfun(@numel, CC.PixelIdxList);
                        [biggest, idx] = max(numPixels);
                        tempMask(CC.PixelIdxList{idx}) = 10;        
                        mask_unrelyPhase = mask_unrelyPhase | (tempMask ~= 10);                    
                    end

                    if strcmp(Params.UnwrappingMethodsDict{Params.UnwrappingMethod}, 'ROMEO')
                        mask_unrelyPhase = mask_unrelyPhase | (handles.phase_quality_map < Params.romeo_phasequality_thresh); 
                    end

                    if strcmp(Params.UnwrappingMethodsDict{Params.UnwrappingMethod}, 'Path') && isfield(Params, 'phasequality_thresh')
                        mask_unrelyPhase = mask_unrelyPhase | (handles.phase_quality_map < Params.phasequality_thresh); 
                    end

                    if isfield(Params, 'maskHs')
                        if Params.maskHs == 1
                            %%% ----- hacking mask 1, to remove more unreliable phase
                            [~, echoNumPickhk] = min(abs(Params.TEs - 30e-3)); 
                            unrelyPhaseThresh = pi/2;               % empirical number, pi/2
                            mask_unrelyPhase = mask_unrelyPhase | create_mask_unrelyPhase(GREPhase(:,:,:,echoNumPickhk,1), unrelyPhaseThresh, maskErode);  
                            %%% ------------------
                        end
                    end

                case 'Laplacian' % Laplacian based
                    % if using iRSHARP can keep the original BET mask
                    switch Params.BgRemovalMethodsDict{Params.BgRemoval}
                        case {'VSHARP','PDF','LBV+VSHARP'}

                            if ~isfield(handles.Params, 'cluster')
                                unrelyPhaseThresh = pi/2; % empirical number
                            else
                                unrelyPhaseThresh = handles.Params.unrelyPhase0_thresh;
                            end

                            if Params.B0 == 3 && length(Params.TEs)>1
                                [~, echoNumPick] = min(abs(Params.TEs - 30e-3));       
                                
                            elseif Params.B0 == 7 && length(Params.TEs)>1
                                [~, echoNumPick] = min(abs(Params.TEs - 18e-3));       
                                                                
                            elseif Params.B0 == 9.4 && length(Params.TEs)>1
                                [~, echoNumPick] = min(abs(Params.TEs - 12e-3));       
                            else
                                echoNumPick = 1;            
                                unrelyPhaseThresh = pi/4;    
                            end             
                            mask_unrelyPhase = create_mask_unrelyPhase(GREPhase(:,:,:,echoNumPick,1), unrelyPhaseThresh, maskErode);       
                        
                            % add
                        
                        case {'iRSHARP'}
                            mask_unrelyPhase = zeros(size(maskErode));
                        otherwise
                            error('unknown option.')
                    end

                otherwise
                    error('Unknown unwrapping method.');
            end
            
            % extra step to remove unreliable phase around brain boundary
            if isfield(handles.Params, 'unrelyPhase1_thresh') && isfield(handles.Params, 'cluster')
                unrelyPhase_thresh = handles.Params.unrelyPhase1_thresh;
            else
                unrelyPhase_thresh = 0.5; % for 3T, TE~20ms
            end
            mask_unrelyPhase2 = create_mask_unrelyPhase2(GREPhase(:,:,:,echoNumPick(1),1).*maskErode, maskErode, Params, unrelyPhase_thresh);
            mask_unrelyPhase2 = imdilate(mask_unrelyPhase2, strel('disk', 1));

            % combine unreliable phase mask
            mask_unrelyPhase = (mask_unrelyPhase | mask_unrelyPhase2 | mask_intrinsic);        
            maskErode = (maskErode.*(maskErode - mask_unrelyPhase) > 0);

            % normal procedure
            maskErode = imerode3dslice(maskErode, strel('disk', 1));
            maskErode = imdilate3dslice(maskErode, strel('disk', 1)); 

            if ~isfield(Params, 'maskHs')
                maskErode = imfill3(maskErode);  

            elseif Params.maskHs == 1
                % for hacking mask 1
                maskErode(floor(0.3*size(maskErode, 1)):end,:,1:floor(0.7*size(maskErode, 3))) = ...
                            imfill3(maskErode(floor(0.3*size(maskErode, 1)):end,:,1:floor(0.7*size(maskErode, 3))));            
           
            elseif Params.maskHs == 2
                % For hacking mask 2, do two imfill3, for hemarrage
                % patients? add two-pass option?
                maskErode = imfill3(maskErode);   
                maskErode = imdilate3dslice(maskErode, strel('disk', 1));
                maskErode = imfill3(maskErode);   
                maskErode = imerode3dslice(maskErode, strel('disk', 1));
            else
                % skip imfill3
            end

            disp('Done')   
            
            if ~isfield(handles.Params, 'cluster')  % GUI only
                multiWaitbar(textWaitbar, 0.8);        
            end
            
            % Erosion
            if max(Params.voxSize) > 4*min(Params.voxSize)
                maskErode = imerode3dslice(maskErode, strel('disk', double(floor(Params.ErodeRadius./min(Params.voxSize)))));          
            else
                maskErode = imerode3(maskErode, floor(Params.ErodeRadius./min(Params.voxSize)), 1);          
            end
            maskErode = maskErode.*(maskErode - mask_intrinsic) > 0;

            % save maskBET
            if ~isfield(handles.Params, 'cluster')  % GUI only
                multiWaitbar(textWaitbar, 0.9 );
            end
            save([outputFile '.mat'], 'maskErode')

            % Save NIFTI
            saveNII(maskErode.*1, outputFile, Params, 1);    
            clear maskBET1 maskBET2 maskBET nii
            
            if ~isfield(handles.Params, 'cluster')  % GUI only
                multiWaitbar( 'CloseAll' );
            else
                disp('Done.')
                writelog(handles.logfile, 'Done. \n');
            end

            % Mark that we made new mask and so we need to do new BG/QSM
            handles.Params.MaskWasMade = 1;
        else
            %% Backup in case FSL does not work!!
            % Is this correct threshold?
            isCorrect = false;
            while(~isCorrect)
                % Do the first masking
                [maskErode,~,~] = create_mask_erode(GREMag(:,:,:,Params.SaveEcho(1),1).*1, double(floor(Params.ErodeRadius./min(Params.voxSize))), Params.MaskThreshold);

                % Show
                if ~isfield(handles.Params, 'cluster')  % GUI only
                    handles.CurrentImage    = GREPhase(:,:,:,Params.SaveEcho(1),1).*maskErode;
                    isNewImage = 1;
                    ShowImage;
                    drawnow;

                    % Ask
                    correctCheck = questdlg(sprintf('FSL BET not found.. Using back-up function.\nMASKING QUALITY IS BEST WHEN USING FSL BET!\nIs the brain correctly extracted or should the threshold be adapted?'),...
                        'Masking','Correct','Change threshold','Stop reconstuction','Correct');
                    
                    % Switch
                    switch correctCheck
                        case 'Correct'
                            isCorrect = true;
                        case 'Change threshold'
                            newThres = inputdlg('Enter new threshold (1-100):', 'Change threshold', [1,50], {num2str(Params.MaskThreshold)});
                            Params.MaskThreshold = str2double(newThres{:});
                        case 'Stop reconstuction'
                            hasCanceled = true;
                            HandleStopReconstruction;
                        otherwise
                            hasCanceled = true;
                            HandleStopReconstruction;
                    end
                else
                    disp('FSL BET not found. Back-up masking used, check with caution...')
                end
            end
            % Mark that we made new mask and so we need to do new BG/QSM
            handles.Params.MaskWasMade = 1;
            % Not saving..
        end
    
    else
        % getting mask for the whole head: single_step_flag == true
        int_mag_max = max(reshape(GREMag(:,:,:,1), 1, [])); % use the first echo with shortest TE 
        int_mag_threshold = 0.18;                            % may need to adjust
        maskErode = GREMag(:,:,:,1) > int_mag_threshold*int_mag_max;
        
        % save mask
        save([outputFile '.mat'], 'maskErode')

        % Save NIFTI
        saveNII(maskErode.*1, outputFile, Params, 1); 

        % Mark that we made new mask and so we need to do new BG/QSM
        handles.Params.MaskWasMade = 1;
    end
end

% Save to workspace
handles.maskErode = maskErode;

if ~isfield(handles.Params, 'cluster')  % GUI only
    % Show
    if Params.EchoAvg > 0
        handles = LoadImage(hObject, handles, GREPhase(:,:,:,Params.echoNums(1),1).*maskErode, 'Phase (rad)');
    else
        handles = LoadImage(hObject, handles, GREPhase(:,:,:,1,1).*maskErode, 'Phase (rad)');
    end
    
    % Save
    guidata(hObject, handles);    

    % Update panel
    set(handles.PanelStep2, 'HighlightColor', [0 0.5 0]);
    set(handles.ButtonShowMask, 'Enable', 'On')
end

% Update table
UpdateTable(handles, 'Completed 2 of 4');