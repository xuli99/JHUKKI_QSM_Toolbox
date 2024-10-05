%% Author: Jiri van Bergen and Xu Li
% Affiliation: Radiology @ JHU - Kirby Center
% Contact via xuli@mri.jhu.edu

% Updated by Xu Li, 2019-07-09
% Updated by X.L., add in msSFCR for ARIC processing, 2019-08-21
% Updated by X.L., removed msSFCR, added NDI & TFI & FANSI for testing
% Updated by X.L., 2021-05-26, added nSFCR and wrapper
% Updated 2021-06-28 X.L., cluster version
% Updated 2021-08-03 X.L., add extra phase reliability mask
% Updated 2021-08-20, X.L., bug fix for dyanmic data
% Updated 2021-09-12, X.L., bug fix for padding in SFCR
% Updated 2022-03-22, X.L. added option to use RAS NIFTI
% Updated 2022-06-24, X.L. added Params check for kernel calculation
% Updated 2024-05-28, X.L., fixed weighting bug for MEDI
% Updated 2024-06-06, X.L., fixed QSMSettingsFile problem for cluster array 
% Updated 2024-10-03, X.L., added option for load CSF mask 

%% Get variables
Params      = handles.Params;
GREMag      = handles.GREMag;
maskErode   = handles.maskErode;
freqMap     = handles.freqMap;

if prod(Params.sizeVol) > 1e7 % in case of ultra high resolution
    datatype = 'single';
else
    datatype = 'double';
end

% extra step to remove unreliable phase around brain boundary
if isfield(handles.Params, 'unrelyPhase2_thresh') && isfield(handles.Params, 'cluster')
    unrelyPhase2_thresh = handles.Params.unrelyPhase2_thresh;
else
    unrelyPhase2_thresh = 0.5;
end
disp('2nd check for unreliable phase ...')
% freqMap may have multiple dyanmics, use dyanmic 1
mask_unrelyPhase2 = create_mask_unrelyPhase2(freqMap(:,:,:,1,1)*2*pi*mean(Params.TEs), maskErode, Params, unrelyPhase2_thresh);
mask_unrelyPhase2 = imdilate(mask_unrelyPhase2, strel('disk', 1));
maskErode2 = (maskErode.*(maskErode - mask_unrelyPhase2) > 0);
% normal procedure
maskErode2 = imerode3dslice(maskErode2, strel('disk', 1));
maskErode2 = imdilate3dslice(maskErode2, strel('disk', 1)); 
maskErode = maskErode.*imfill3(maskErode2); freqMap = freqMap.*maskErode;
disp('Done.')

% estimate noise level here
noise_level = calfieldnoise(GREMag(:,:,:,:,1), handles.GREPhaseRaw(:,:,:,1,1), maskErode);

% Calculate base values
deltaB          = freqMap./(Params.gamma*Params.B0)*1e6;   % in ppm, field map
clear freqMap GREPhase GREPhaseRaw GREPhaseRaw1

deltaB = cast(deltaB, datatype);
GREMag = cast(GREMag, datatype);

%% switch for getting R2* and AutoRef (SFCR+0)
if Params.R2starFlag == 1
    % fitting R2*
    R2starFile = [Params.FileBaseName, '_R2star'];         
    if ~exist(fullfile(cd, [R2starFile, '.mat']), 'file') && size(GREMag,4)>=2         
        
        if Params.single_step_flag
            BrainMaskFilename = [Params.FileBaseName, '_mask_whole_head.nii.gz'];
        else
            fname1 = [Params.FileBaseName '_GREMag', num2str(Params.SaveEcho(1))];
            BrainMaskFilename = [fname1, '_brain_mask.nii.gz'];
        end
        
        nii = load_untouch_nii(BrainMaskFilename);  
        % for RAS NIFTI
        nii = nii_img_update(nii, Params);

        maskBET = permute(nii.img, [2,1,3]);            
        maskBET = maskBET > 0;
        
        % ------------- calculate R2* Maps and save file
        disp('Fitting R2* ...')
        tic
            R2starmap = R2star_ARLO(GREMag, Params.TEs, maskBET);        % ARLO faster            
        toc
        disp('Done.')
        
        if isfield(handles.Params, 'cluster')  % cluster only
            writelog(handles.logfile, ['R2* fitting Done. \n']);
        end
        
        R2starMax = 500;
        R2starmap(R2starmap < 0) = 0;
        R2starmap(R2starmap > R2starMax) = R2starMax;
        
        save(fullfile(cd, [R2starFile, '.mat']), 'R2starmap', 'Params', 'maskBET');
        saveNII(R2starmap, R2starFile, Params);

        clear maskBET nii fname1 R2starmap S0map R2starFile
    else
        if exist(fullfile(cd, [R2starFile, '.mat']), 'file')
            logtxt = 'R2* fitting has been done already.';
            disp(logtxt)
            if isfield(handles.Params, 'cluster')  % cluster only
                writelog(handles.logfile, [logtxt, ' \n']);
            end
        else
            error('Cannot fit R2* with single echo data.')
        end
    end
    
    if (Params.AutoRefFlag == 1) && (Params.B0 == 3)    % if AutoRef and 3T
       
        fileName = [ Params.FileBaseName, '_GREMag1'];
        if ~exist(fullfile(cd, [fileName, '.nii']), 'file') || ~exist(fullfile(cd, [fileName, '.nii.gz']), 'file')
            saveNII(GREMag(:,:,:,1).*1, fileName, Params, 1);
        end

        BrainMaskFilename = [fileName, '_brain_mask.nii.gz'];
        if ~exist(fullfile(cd, BrainMaskFilename), 'file')
            inputstring1 = [Params.FSLFolder, 'bet ', fileName, '.nii ' fileName, '_brain', ' -f 0.5 -g 0 -m' ];
            disp(inputstring1)
            system(inputstring1);
        end
        disp('bet done.')

        % do fast on GREMag 1 to improve CSF mask
        BETFilename = [fileName, '_brain'];   
        FASTFilename = [fileName, '_brain_seg.nii.gz'];
        if ~exist(fullfile(cd, FASTFilename), 'file') 
            inputstring1 = [Params.FSLFolder, 'fast -t 1 -n 2 -H 0.1 -I 4 -l 20.0 -o ', BETFilename, ' ', BETFilename];
            disp(inputstring1)
            system(inputstring1);
        end
        disp('fast done.')
        
        if isfield(handles.Params, 'cluster')  % cluster only
            writelog(handles.logfile, ['FSL FAST Done. \n']);
        end
        clear BETFilename FASTFilename BrainMaskFilename
    end
 
end

%% load in data
% -------------------- if AutoRef or single_step_flag
if Params.AutoRefFlag == 1 || Params.single_step_flag
    % -------------  load in R2* and fsl fast segmetnation
    R2starFile = [Params.FileBaseName, '_R2star'];     
    if (exist(fullfile(cd, [R2starFile, '.mat']), 'file') == 2)
        disp('R2* map exist.')
        temp = load(R2starFile);
        R2starMap = temp.R2starmap;
        clear temp
        R2starMap = cast(R2starMap, datatype);
    end
end

% % -------------------- if single_step_flag, get an extra BET mask
% if Params.single_step_flag
%     fileName = [ Params.FileBaseName, '_GREMag1'];
%     if ~exist(fullfile(cd, [fileName, '.nii.gz']), 'file')
%         saveNII(GREMag(:,:,:,1).*1, fileName, Params, 1);
%     end
% 
%     BrainMaskFilename = [fileName, '_brain_mask.nii.gz'];
%     if ~exist(fullfile(cd, BrainMaskFilename), 'file')
%         inputstring1 = [Params.FSLFolder, 'bet ', fileName, '.nii ' fileName, '_brain', ' -f 0.5 -g 0 -m' ];
%         disp(inputstring1)
%         system(inputstring1);
%     end
%     disp('bet done.')
%     
%     nii = load_untouch_nii(BrainMaskFilename); 
%     nii = nii_img_update(nii, Params);
%     maskBET = permute(nii.img, [2,1,3]);            
%     maskBET = maskBET > 0;
% end

% -------------------- if AutoRef
if Params.AutoRefFlag == 1
    GREMagSegFile = [Params.FileBaseName, '_brain_mixeltype.nii.gz'];             % old version
    GREMagSegFile2 = [Params.FileBaseName, '_GREMag1_brain_mixeltype.nii.gz'];    % new version
    GREMagSegFile3 = [Params.FileBaseName, '_Seg_brain_mixeltype.nii.gz'];        % other option

    CSFmaskFileFlag = 0;
    CSFmaskFile = [];
    if isfield(Params, 'FileCSFmask_App')
        CSFmaskFile = [Params.FileBaseName, Params.FileCSFmask_App];              % load CSFmask directly, if without R2*
    elseif isfield(Params, 'FileCSFmask')
        CSFmaskFile = Params.FileCSFmask;
    end

    if ~isempty(CSFmaskFile)
        if (exist(fullfile(cd, CSFmaskFile), 'file') == 2)
            CSFmaskFileFlag = 1;
        else
            disp([CSFmaskFile, ' not found.'])
        end
    end

    GREMagSegFileFlag = 1;
    if (exist(fullfile(cd, GREMagSegFile), 'file') == 2) 
        GREMagSegFileTarget = GREMagSegFile;            % .nii.gz
    
    elseif (exist(fullfile(cd, GREMagSegFile(1:end-3)), 'file') == 2)
        GREMagSegFileTarget = GREMagSegFile(1:end-3);   % .nii

    elseif (exist(fullfile(cd, GREMagSegFile2), 'file') == 2) 
        GREMagSegFileTarget = GREMagSegFile2;

    elseif (exist(fullfile(cd, GREMagSegFile2(1:end-3)), 'file') == 2) 
        GREMagSegFileTarget = GREMagSegFile2(1:end-3);   % .nii

    elseif (exist(fullfile(cd, GREMagSegFile3), 'file') == 2) 
        GREMagSegFileTarget = GREMagSegFile3;            % .nii.gz
    
    else
        GREMagSegFileFlag = 0;
    end

    if GREMagSegFileFlag == 1
        disp('GREMagSeg file exists.')
        nii = load_untouch_nii(GREMagSegFileTarget);
        % for RAS NIFTI
        nii = nii_img_update(nii, Params);        
        GREMagSeg = permute(nii.img, [2,1,3]);   
    end

    if CSFmaskFileFlag == 1
        disp('CSFmask file exists. Load directly.')
        nii = load_untouch_nii(CSFmaskFile);
        % for RAS NIFTI
        nii = nii_img_update(nii, Params);        
        CSFmask1 = cast(permute(nii.img, [2,1,3]), 'double');
        CSFmask1 = CSFmask1 > 0;
    end
end

% ----------------------
GREMag = GREMag./noise_level;

if ~strcmpi(Params.QSMSolverDict{Params.QSMSolver}, 'TFI')
    deltaB          = deltaB.*maskErode;                       
end

%% make DPWeight if needed
switch Params.QSMSolverDict{Params.QSMSolver} 
    case {'MEDI', 'SFCR', 'nSFCR', 'FANSI', 'NDI', 'TFI'}
        DPWeight = zeros([Params.sizeVol, 1, Params.nDynamics]);
        if isfield(handles, 'DPWeight')
            DPWeight = handles.DPWeight.*maskErode;            % loaded weight may not be normalized
            switch Params.QSMSolverDict{Params.QSMSolver} 
                case {'MEDI', 'SFCR', 'NDI', 'TFI'}         
                    DPWeight = DPWeight./mean(DPWeight(maskErode>0));   % mean normalization, v3.0
                case {'FANSI','nSFCR'}
                    DPWeight = DPWeight./max(DPWeight(maskErode>0));    % max normalization
            end
        else
            for dynamic_ind = 1:Params.nDynamics
                DPWeight(:,:,:,1,dynamic_ind) = ...
                    sqrt(sum(single(GREMag(:,:,:,Params.echoNums,dynamic_ind).*reshape(Params.TEs(Params.echoNums), [1,1,1,length(Params.echoNums)])).^2, 4));
                DPWeight(:,:,:,1,dynamic_ind) = DPWeight(:,:,:,1,dynamic_ind).*maskErode; 
                temp = DPWeight(:,:,:,1,dynamic_ind);

                switch Params.QSMSolverDict{Params.QSMSolver} 
                    case {'MEDI', 'SFCR', 'NDI', 'TFI'}  
                        DPWeight(:,:,:,1,dynamic_ind) = DPWeight(:,:,:,1,dynamic_ind)./mean(temp(maskErode>0));  % mean normalization, v3.0
                    case {'FANSI','nSFCR'}    
                        DPWeight(:,:,:,1,dynamic_ind) = DPWeight(:,:,:,1,dynamic_ind)./max(temp(maskErode>0));   % max normalization
                end
            end
        end
        clear GREMag  
    otherwise
        clear GREMag
end

%% determin output file name
if Params.AutoRefFlag == 1 && contains(Params.QSMSolverDict{Params.QSMSolver}, ['SFCR'])
    StringApp1 = ['_chi_', Params.QSMSolverDict{Params.QSMSolver}, '+0'];
    WaitBarMsgloading = [Params.QSMSolverDict{Params.QSMSolver}, '+0'];
else
    StringApp1 = ['_chi_', Params.QSMSolverDict{Params.QSMSolver}];
    WaitBarMsgloading = Params.QSMSolverDict{Params.QSMSolver};
end

% naming convention
if strcmp(Params.UnwrappingMethodsDict{Params.UnwrappingMethod}, 'NonlinearFit + Path')
    StringApp2 = '_NLSlope';
elseif Params.EchoAvg > 0
    StringApp2 = '_Avg';
else
    StringApp2 = '_Slope';
end

if (length(Params.echoNums) > 1)            
    StringApp4 = ['-', num2str(Params.echoNums(end))];
else
    StringApp4 = '';
end

StringApp3 = ['_echo', num2str(Params.echoNums(1))];

if Params.echoStep > 1      % For Single Echo, Params.echoStep=0;
    StringApp5 = ['_s', num2str(Params.echoStep)];
else
    StringApp5 = '';
end

outputFile = [Params.FileBaseName, StringApp1, StringApp2];     % version 1
% outputFile = [Params.FileBaseName, StringApp1, StringApp2, StringApp3, StringApp4, StringApp5];

%% QSM solver
if(exist(fullfile(cd, [outputFile, '.mat']), 'file') == 2 && ...
        Params.saveOutput && ...
        ~Params.BgWasRemoved) % Re-DO QSM if new background removal was performed!!
    textWaitbar = ['Loading previously ', WaitBarMsgloading,  ' calculated QSM...'];
    if ~isfield(handles.Params, 'cluster')  % GUI only
        multiWaitbar( textWaitbar, 'Busy');
    else
        disp(textWaitbar);
        writelog(handles.logfile, [textWaitbar, '...']);
    end
    load(outputFile)
    if ~isfield(handles.Params, 'cluster')  % GUI only
        multiWaitbar('CloseAll');
    else
        disp('Done.')
        writelog(handles.logfile, 'Done. \n');
    end
else            
    textWaitbar = ['Calculating QSM for '  num2str(Params.nDynamics), ' dynamics.'];
    if ~isfield(handles.Params, 'cluster')  % GUI only
        multiWaitbar(textWaitbar, 0, 'CanCancel', 'On');
    else
        disp(textWaitbar);
        writelog(handles.logfile, [textWaitbar, '...']);
    end
    chi_res = zeros(size(deltaB));
    % checking Params for kernel
    Params.fov = Params.fov(:)';
    Params.sizeVol = Params.sizeVol(:)';
    Params.voxSize = Params.voxSize(:)';
    
    switch Params.QSMSolverDict{Params.QSMSolver} 
        case 'iLSQR'        
            D = conv_kernel_rot_c0(Params, Params.TAng, datatype); 
            
            for dynamic_ind = 1:Params.nDynamics
                
                % initial estimation using LSQR2
                disp('initial LSQR...')
                tol = 0.01;            
                chi_res0 = delta2chi_lsqr2(deltaB(:,:,:,1,dynamic_ind), D, tol, maskErode);                   
                chi_res0 = chi_res0.*maskErode;

                disp('estimating a-priori using iTKD ...')
                tol = 0.04;       
                [chi_ap] = delta2chi_iTKD(deltaB(:,:,:,1,dynamic_ind), D, maskErode, tol);    

                % further remove possible streaking artifact  if necessary            
                thresh_SAR = 0.15;
                [chi_ap, chi_SA_ap] = delta2chi_SAR(deltaB(:,:,:,1,dynamic_ind), D, chi_ap, maskErode, thresh_SAR, chi_ap);  

                % Final estimate streaking artifact in LSQR2
                disp('final SAR ...')
                thresh_SAR = 0.15;    
                [chi_res(:,:,:,1,dynamic_ind), chi_SA] = delta2chi_SAR(deltaB(:,:,:,1,dynamic_ind), D, chi_res0, maskErode, thresh_SAR, chi_ap);
                
                if ~isfield(handles.Params, 'cluster')  % GUI only
                    hasCanceled = multiWaitbar(textWaitbar, (dynamic_ind/Params.nDynamics));
                    HandleStopReconstruction;
                end
            end

        case 'TKD' 
            thresh_tkd = 0.2;
            D = conv_kernel_rot_c0(Params, Params.TAng, datatype);
            
            for dynamic_ind = 1:Params.nDynamics
                chi_res(:,:,:,1,dynamic_ind) = delta2chi_tso(deltaB(:,:,:,1,dynamic_ind), D, thresh_tkd);
                chi_res(:,:,:,1,dynamic_ind) = chi_res(:,:,:,1,dynamic_ind).*maskErode;
                if ~isfield(handles.Params, 'cluster')  % GUI only
                    hasCanceled = multiWaitbar(textWaitbar, (dynamic_ind/Params.nDynamics));
                    HandleStopReconstruction;
                end
            end

        case 'iTKD'    
            D = conv_kernel_rot_c0(Params, Params.TAng, datatype);            

            tol = 0.04;       
            for dynamic_ind = 1:Params.nDynamics
                chi_res0 = delta2chi_iTKD(deltaB(:,:,:,1,dynamic_ind), D, maskErode, tol);   
                thresh_SAR = 0.15;
                [chi_res(:,:,:,1,dynamic_ind), chi_SA] = ...
                    delta2chi_SAR(deltaB(:,:,:,1,dynamic_ind), D, chi_res0, maskErode, thresh_SAR, chi_res0);      
                if ~isfield(handles.Params, 'cluster')  % GUI only
                    hasCanceled = multiWaitbar(textWaitbar, (dynamic_ind/Params.nDynamics));
                    HandleStopReconstruction;
                end
            end
            
        case 'MEDI'     
            if isfield(Params, 'QSM_MEDIlambda')
                lambda = Params.QSM_MEDIlambda;
            else
                lambda = 2500;  % default
            end
            edgePer = 3*0.3;    % Edge voxel percentage
            merit = 1;          % fine tuning
            
            for dynamic_ind = 1:Params.nDynamics
                
                [chi_res(:,:,:,1,dynamic_ind), regv, datav] = ...
                    delta2chi_MEDI(deltaB(:,:,:,1,dynamic_ind), Params, DPWeight(:,:,:,1,dynamic_ind), maskErode, lambda, merit, edgePer);  
                if ~isfield(handles.Params, 'cluster')  % GUI only
                    hasCanceled = multiWaitbar(textWaitbar, (dynamic_ind/Params.nDynamics));
                    HandleStopReconstruction;               
                end
            end
            
        case {'SFCR'}          %% modifed SFCR, SFCR+0
            % padding
            padsize = [0, 0, 0];
            % padsize = [12, 12, 12]; % for reduced FOV use
            Params.sizeVol = Params.sizeVol + 2*padsize;
            Params.fov = Params.sizeVol.*Params.voxSize;
            D = conv_kernel_rot_c0(Params, Params.TAng, datatype);
            D = ifftshift(D);
            
            % Parameters set for modified SFCR code (2&3)
            MagEdgePer = 0.3;  % Edge voxel percentage for magnitude 
            chiEdgePer = 0.35;
            lambdaSet.beta = 10;
            lambdaSet.lambda1_M = 500;      
            lambdaSet.lambda2_M = 0;
            lambdaSet.lambda1_S = 1000;     
            lambdaSet.lambda2_S = 0;
            
            if (Params.AutoRefFlag == 1)
                if (exist('R2starMap', 'var') == 1)
                    % Use R2* to get CSF masking
                    lambdaSet.R2sThresh = 5;        % 5 Hz for extracting central CSF region for automatic CSF referencing               
                    [CSFmask1] = CSFmaskThresh(R2starMap, lambdaSet.R2sThresh, maskErode, Params.voxSize);                
                    
                    if (exist('GREMagSeg', 'var') == 1)
                        CSFmask2 = (GREMagSeg == 0) & maskErode;
                        disp('updating CSF mask based on R2* with FSL Segmentation')
                        % --- Good choice if with small ventricles (e.g. in RLS study),
                        % use with caution for other cases
                        CSFmask1 = R2starMap < lambdaSet.R2sThresh;
    
                    else
                        CSFmask2 = maskErode;
                    end
                    lambdaSet.maskSS = CSFmask1 & CSFmask2;
                    lambdaSet.lambda2_M = lambdaSet.lambda1_M./5;
                    lambdaSet.lambda2_S = lambdaSet.lambda1_S./5;

                elseif CSFmaskFileFlag == 1
                    lambdaSet.maskSS = CSFmask1;
                    lambdaSet.lambda2_M = lambdaSet.lambda1_M./5;
                    lambdaSet.lambda2_S = lambdaSet.lambda1_S./5;

                else
                    disp('No R2* map or CSF mask, cannot do AutoRef to CSF.')
                    if isfield(handles.Params, 'cluster')  % if cluster version, skip the current case and continue
                        disp('For cluster version, skip the current case and continue ...')
                        writelog(handles.logfile, 'AutoRef Error, no R2* or CSF mask. QSM Skipped. \n');
                        handles = UpdateTable(handles, 'Error with AutoRef. Skipped.');
                        return
                    end

                end

            end
            
            deltaB = padarray(deltaB, padsize);
            maskErode = padarray(maskErode, padsize);
            DPWeight = padarray(DPWeight, padsize);
            chi_res = padarray(chi_res, padsize);
            if isfield(lambdaSet, 'maskSS')
                lambdaSet.maskSS = padarray(lambdaSet.maskSS, padsize);
            end
            
            % extra parameters
            thre_tkd = 0.18;                                % 
            lambdaSet.datatype      = datatype;
            lambdaSet.gamma         = Params.gamma;         % constants
            lambdaSet.B0            = Params.B0;
            lambdaSet.TEs           = Params.TEs;
            lambdaSet.codeswitch    = 2;                    % adjust
            % 2 for original version; 0 for datawith large susceptibiltiy sources
            
            for dynamic_ind = 1:Params.nDynamics
                if strcmp(Params.QSMSolverDict{Params.QSMSolver}, 'SFCR')
                    [chi_res(:,:,:,1,dynamic_ind), SB_residual, SB_regularization] = ...
                        delta2chi_SFCR(deltaB(:,:,:,1,dynamic_ind), D, maskErode, DPWeight(:,:,:,1,dynamic_ind), ...
                        MagEdgePer, chiEdgePer, lambdaSet, thre_tkd, Params.AutoRefFlag);
                end
                if ~isfield(handles.Params, 'cluster')  % GUI only
                    hasCanceled = multiWaitbar(textWaitbar, (dynamic_ind/Params.nDynamics));
                    HandleStopReconstruction;
                end
            end
            
            chi_res = chi_res(padsize(1)+1:end-padsize(1), padsize(2)+1:end-padsize(2), padsize(3)+1:end-padsize(3), :, :);
            maskErode = maskErode(padsize(1)+1:end-padsize(1), padsize(2)+1:end-padsize(2), padsize(3)+1:end-padsize(3), :, :);
            if isfield(lambdaSet, 'maskSS')
                lambdaSet.maskSS = lambdaSet.maskSS(padsize(1)+1:end-padsize(1), padsize(2)+1:end-padsize(2), padsize(3)+1:end-padsize(3), :, :);
            end
            Params.sizeVol = Params.sizeVol - 2*padsize;
            Params.fov = Params.sizeVol.*Params.voxSize;
            
        case 'nSFCR'
            % cluster setup may change these parameters
            if isfield(handles.Params, 'nSFCRparams') && isfield(handles.Params, 'cluster')
                if isfield(handles.Params.nSFCRparams, 'L1orL2')
                    nSFCRparams.L1orL2 = handles.Params.nSFCRparams.L1orL2;
                end
                
                if isfield(handles.Params.nSFCRparams, 'nlM')
                    nSFCRparams.nlM = handles.Params.nSFCRparams.nlM;
                end
                
                if isfield(handles.Params.nSFCRparams, 'TV')
                    nSFCRparams.TV = handles.Params.nSFCRparams.TV;
                end
                
            else
                % defaults
                nSFCRparams.L1orL2 = 1;          % 1 for L1 and 2 for L2
                nSFCRparams.nlM = 1;             % 0 use old linear M-step, 1 use nSFCR
                nSFCRparams.TV  = 0;             % 0 use both M-step and F-step as in SFCR, 1 use TV and only S-step
            end
                
            nSFCRparams.padsize = [12, 12, 12]; % pad for smaller FOV
            nSFCRparams.Params = Params;
            nSFCRparams.datatype = datatype;
            nSFCRparams.lambda2 = 0;            % AutoRef L2 penalty, default
            nSFCRparams.maskRef = maskErode;
            
            if (Params.AutoRefFlag == 1) 
                if exist('R2starMap', 'var') == 1
                    nSFCRparams.R2sThresh = 5;        % 5 Hz for extracting central CSF region for automatic CSF referencing               
                    [CSFmask1] = CSFmaskThresh(R2starMap, nSFCRparams.R2sThresh, maskErode, Params.voxSize);                

                    if (exist('GREMagSeg', 'var') == 1)
                        CSFmask2 = (GREMagSeg == 0) & maskErode;
                        disp('updating CSF mask based on R2* with FSL Segmentation')
                        CSFmask1 = R2starMap < nSFCRparams.R2sThresh;
                    else
                        CSFmask2 = maskErode;
                    end
                    nSFCRparams.maskRef = CSFmask1 & CSFmask2;
                    nSFCRparams.lambda2 = 0.5;               
                    
                else
                    if (exist('GREMagSeg', 'var') == 1)
                        CSFmask = (GREMagSeg == 0) & maskErode;
                        disp('loading CSF mask from FSL Segmentation')
                        
                        GREMag  = handles.GREMag;
                        maskExclude = GREMag < median(GREMag(GREMagSeg > 0)); % T2/T2* based for Heme OR Cal not seen in T1 with low T2 signal
                        CSFmask = CSFmask.*~maskExclude;
                        clear GREMag    
                    end

                    nSFCRparams.maskRef = CSFmask;
                    nSFCRparams.lambda2 = 0.5;               
                end

            end
            
            for dynamic_ind = 1:Params.nDynamics
                [chi_res(:,:,:,1,dynamic_ind), chi_res_M] = delta2chi_nSFCR(deltaB(:,:,:,1,dynamic_ind), maskErode, DPWeight(:,:,:,1,dynamic_ind), nSFCRparams);
                if ~isfield(handles.Params, 'cluster')  % GUI only
                    hasCanceled = multiWaitbar(textWaitbar, (dynamic_ind/Params.nDynamics));
                    HandleStopReconstruction;
                end
            end
            
        case 'FANSI'
            % NOTE: has to be even dimension
            N = size(maskErode); dimsOdd = mod(N, 2);
            deltaB = padarray(deltaB, dimsOdd, 'replicate', 'post');
            maskErode = padarray(maskErode, dimsOdd, 'replicate', 'post');
            DPWeight = padarray(DPWeight, dimsOdd, 'replicate', 'post');
            chi_res = padarray(chi_res, dimsOdd, 'replicate', 'post');
            
            Params.sizeVol = Params.sizeVol + dimsOdd;
            Params.fov = Params.sizeVol.*Params.voxSize;
    
            AlphaSweep = 1;
            for ka = 1:length(AlphaSweep)               
 
                L1orL2 = 1;             % 1 for L1 (QSM) and 2 for L2 (FANSI)
                if L1orL2 == 1
                    alpha1 = 5e-3;    % AlphaSweep(ka); % for nlL1TV, better to do L-curve searching/tuning of alpha1 & mu1 & mu2=1
                    FANSIparams.merit = 0;
                    FANSIparams.tol_update = 1;                % default 1
                    FANSIparams.maxOuterIter = 50;             % default 150
                elseif L1orL2 == 2
                    alpha1 = 2.5e-4;   % AlphaSweep(ka); % for nlTV, default 2.5e-4, better to do L-curve searching/tuning of alpha1 & mu1 & mu2=1
                    FANSIparams.merit = 1;                     % merit makes hemorrage looks larger
                    FANSIparams.tol_update = 0.1;              % default 0.1
                    FANSIparams.maxOuterIter = 150;            % default 150
                end
                
                FANSIparams.D = ifftshift(conv_kernel_rot_c0(Params, Params.TAng, datatype));
                FANSIparams.alpha1 = alpha1;             % gradient L1 penalty
                FANSIparams.mu1 = 100*alpha1;            % gradient consistency weight
                FANSIparams.mask = maskErode;
                phasescale = 2*pi*Params.gamma*Params.B0*mean(Params.TEs)*1e-6; % ppm to radian

                for dynamic_ind = 1:Params.nDynamics
                    FANSIparams.input = deltaB(:,:,:,1,dynamic_ind)*phasescale;
                    FANSIparams.weight = DPWeight(:,:,:,1,dynamic_ind);
                    % FANSIparams.regweight = gradient_mask_all(DPWeight, maskErode, 0.3);
                    if L1orL2 == 1 
                        out = delta2chi_FANSI_nlL1TV(FANSIparams);  % nonlinear L1 data fidelity
                        % outputFile = [outputFile, '_L1'];           % 
                    elseif L1orL2 == 2
                        out = delta2chi_FANSI_nlTV(FANSIparams);      % nonlinear L2 data fidelity 
                        % outputFile = [outputFile, '_L2'];           % 
                    end            
                    chi_res(:,:,:,1,dynamic_ind) = out.x/phasescale.*maskErode;
                    
                    if ~isfield(handles.Params, 'cluster')  % GUI only
                        hasCanceled = multiWaitbar(textWaitbar, (dynamic_ind/Params.nDynamics));
                        HandleStopReconstruction;
                    end
                end

                % change back
                chi_res = chi_res(1:N(1), 1:N(2), 1:N(3), :, :);
                Params.sizeVol = Params.sizeVol - dimsOdd;
                Params.fov = Params.sizeVol.*Params.voxSize;
                
                % % for parameter sweep 
                % save([outputFile, '_L', num2str(L1orL2), '_alpha', num2str(ka), '.mat'], 'out', 'FANSIparams'); 
             end
            
        case 'NDI'
            NDIParams.step_size = 1;       % gradient descent step size
            NDIParams.num_iter = 20; 
            NDIParams.tik  = 1;            % adding tikhonov regularization
            NDIParams.lambda = 1e-3;       % make if small
            NDIParams.tol = 1;           % for 1 orientation
            NDIParams.datatype = datatype;
            NDIParams.gamma = Params.gamma;          % constants
            NDIParams.B0 = Params.B0;
            NDIParams.TEs = Params.TEs;
            
            D = conv_kernel_rot_c0(Params, Params.TAng, datatype);
            D = ifftshift(D);

            for dynamic_ind = 1:Params.nDynamics                
                [chi_res(:,:,:,1,dynamic_ind), flag, relres, iter] = ...
                    delta2chi_NDI(deltaB(:,:,:,1,dynamic_ind), D, DPWeight(:,:,:,1,dynamic_ind), maskErode, NDIParams);  
                
                if ~isfield(handles.Params, 'cluster')  % GUI only
                    hasCanceled = multiWaitbar(textWaitbar, (dynamic_ind/Params.nDynamics));
                    HandleStopReconstruction;  
                end
            end
            
        case 'TFI'
            % test version, only match with NLFpath phase unwrapping
            
            % deltaB is the total field (in ppm), maskErode is whole head mask M
            TFI_params.R2star_thresh = 30;             % Hz, threshold for detecting soft tissue, for old TFI
            TFI_params.Ps            = 30;             % Preconditioning value
            TFI_params.P_map = ones(size(deltaB)).*TFI_params.Ps;     % default is 1 in M, otherwise estimated according to chi_est
            
            if (exist('R2starMap', 'var') == 1) 
                %% for old TFI
                % TFI_params.P_map(maskErode & (R2starMap < TFI_params.R2star_thresh)) = 1;  % soft-tissue, excluding hemorrhage etc
                
                %% for auto TFI
                % Auto estimation of TFI P_map in maskErode and ~maskErode
                % with continuous values

                % 1. get D map according to maskErode (soft tissue mask here)
                D_map = bwdist(maskErode); % distance assuming isotropic voxels (need to change for anisotropic voxels)         
                % 2. get chi_est in 1-M using PDF 
                disp('doing fast PDF ...')
                [~, chi_est0] = dipolefit_v2(deltaB, maskErode, DPWeight, Params, 0.001, 5, 0);      % 5cg, no cropping
                % 3. get chi_est in M (LBV + TKD)
                % do LBV
                disp('doing LBV')
                Params.lbv.tol = 0.01;      % default
                Params.lbv.depth = -1;      % default
                Params.lbv.peel = 3;        % similar to mask erosion
                deltaB_local_est = LBV(deltaB, maskErode, Params.sizeVol, Params.voxSize, ...
                                            Params.lbv.tol, Params.lbv.depth, Params.lbv.peel);
                mask_LBV = abs(deltaB_local_est) > 2*eps;
                
                % do TKD
                disp('doing TKD')
                thresh_tkd = 0.2;
                D = conv_kernel_rot_c0(Params, Params.TAng, datatype);
                chi_est1 = delta2chi_tso(deltaB_local_est, D, thresh_tkd);
                chi_est1 = chi_est1.*mask_LBV;

                chi_est = chi_est0.*~maskErode + chi_est1.*maskErode;
                clear chi_est0 chi_est1
                TFI_params.P_map = TFI_autofit(D_map, chi_est, R2starMap, maskErode, mask_LBV);
                
            else
                %% if without R2*
                TFI_params.P_map(maskErode) = 1;
            end
             
            % other TFI parameters
            TFI_params.apEdge   = 0.4;   % in maskErode
            TFI_params.lambda1  = 5e-4;  % WTV regularization
            TFI_params.lambda2  = 0.1;   % if AutoRef, not done yet
            TFI_params.lambda3  = 1e-4;  % Tikhonov regularization as in LN-QSM, on maskBET
            TFI_params.epsilon  = 1e-6*TFI_params.P_map.^2;
            TFI_params.tol      = 0.01;  % tol for GNCG
            TFI_params.itermax  = 300;
            TFI_params.cgtol    = 0.01;  % cg tol
            TFI_params.cgitermax= 100;   % main iteration
            TFI_params.verbose  = 0;
            
            % calculate TFI with TFI_params.P_map
            D = conv_kernel_rot_c0(Params, Params.TAng, datatype);
            D = ifftshift(D);
            [chi_res, TFI_params] = delta2chi_TFI(deltaB, D, DPWeight, maskErode, TFI_params);
             
            % for parameter sweeping
            % outputFile = [outputFile, 'lambda3_', num2str(TFI_params.lambda3)];
        otherwise
            error('unknown selection of QSM method.')
    end     


    %% save
    switch Params.QSMSolverDict{Params.QSMSolver} 
        case 'iLSQR'
            save([outputFile '.mat'], 'chi_res', 'chi_SA', 'chi_ap', 'tol', 'thresh_SAR', 'maskErode');  
        case 'TKD'
            save([outputFile '.mat'], 'chi_res', 'thresh_tkd', 'maskErode');
        case 'iTKD'
            save([outputFile '.mat'], 'chi_res', 'chi_SA', 'chi_res0', 'thresh_SAR', 'tol', 'maskErode');
        case 'MEDI'
            save([outputFile '.mat'], 'chi_res', 'maskErode', 'Params', 'lambda', 'merit');  
        case 'SFCR'
            save([outputFile '.mat'], 'chi_res', 'maskErode', 'Params', 'MagEdgePer', ...
                    'chiEdgePer', 'lambdaSet', 'SB_residual', 'SB_regularization');
        case 'nSFCR'
            save([outputFile, '.mat'], 'chi_res', 'chi_res_M', 'maskErode', 'Params', 'nSFCRparams');
        case 'FANSI'
            save([outputFile, '.mat'], 'chi_res', 'alpha1', 'maskErode', 'Params');            
        case 'NDI'
            save([outputFile, '.mat'], 'chi_res', 'flag', 'relres', 'iter', 'maskErode', 'NDIParams');
        
        case 'TFI'
            save([outputFile, '.mat'], 'chi_res', 'TFI_params', 'maskErode');
        otherwise
            error('unknown selection of QSM method.')
    end

    if ~isfield(handles.Params, 'cluster')  % GUI only
        multiWaitbar('CloseAll');
    else
        disp('QSM Done.')
        writelog(handles.logfile, 'QSM Done. \n');
    end
    saveNII(squeeze(chi_res), outputFile, Params, 1);
    
    % option to save the final QSM mask
    if isfield(Params, 'SaveQSMmask')
        if Params.SaveQSMmask == 1
            saveNII(maskErode.*1, [outputFile, '_QSMmask'], Params, 1);
        end
    end
end

% Save
handles.chi_res        = chi_res;

if ~isfield(handles.Params, 'cluster')  % GUI only
    guidata(hObject, handles);
    % Display

    handles = LoadImage(hObject, handles, handles.chi_res, 'Susceptibility (ppm)');
    % Update panel
    set(handles.PanelStep4, 'HighlightColor', [0 0.5 0]);
    set(handles.ButtonShowQSM, 'Enable', 'On')
end

% Save Constants file
Params                 = handles.Params;
if ~isfield(Params,'NoSaveQSMsetting')
    % default
    save(handles.Params.QSMSettingsFile, 'Params');
else
    if Params.NoSaveQSMsetting == 0
        save(handles.Params.QSMSettingsFile, 'Params');
    end
end

% Update table
handles = UpdateTable(handles, 'Completed 4 of 4 - Finished!');

if ~isfield(handles.Params, 'cluster')  % GUI only
    % Update buttons
    set([handles.ButtonAddDataset handles.ButtonStartDatasets handles.ButtonLoadDataList], 'Enable', 'On')

    % Enable buttons
    set([handles.ButtonEditEchoes handles.VarB0 handles.VarRadiusDisk handles.VarMaskEchoes handles.VarFSL handles.VarBgRemoval handles.VarSHARPradius handles.VarQSMSolver handles.VarFSLThres], 'Enable', 'On');
end