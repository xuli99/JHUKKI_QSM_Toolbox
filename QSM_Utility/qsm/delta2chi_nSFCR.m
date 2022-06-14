function [chi_res, chi_res_M] = delta2chi_nSFCR(deltaB, maskErode, DPWeight, nSFCRparams)
% chi_res = delta2chi_nSFCR(deltaB, maskErode, DPWeight, nSFCRparams)
% delta2chi_nSFCR main wrapper of nSFCR
% which calls delta2chi_nSFCRL2 & delta2chi_nSFCRL1
% Ref: Bao et al, IEEE TMI, 2016 35(9):2040-50
%      Chen et al, J Neurosci Res, 2019, 97:467:479
%      Milovic, MRM, 2018, 80:814–821 (FANSI)
%      FANSI toolbox:  https://gitlab.com/cmilovic/FANSI-toolbox 
% 
% Authors: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
%
% INPUT:
% deltaB: relative field shift
% D: convolution kernel
% maskErode: brain mask
% DPWeight: normalized magnitude data, SNR weighting
% nSFCRparams: parameters set
%             nSFCRparams.L1orL2 = 2;             % 1 for L1 and 2 for L2
%             nSFCRparams.nlM    = 1;             % 0 use old linear M-step + k-space combination, 1 use nSFCR
%             nSFCRparams.TV     = 0;             % 0 use both M-step and F-step as in SFCR, 1 use TV and only S-step 
%             nSFCRparams.padsize = [12, 12, 12]; 
%             nSFCRparams.Params                  % original Params in the toolbox
%             nSFCRparams.maskRef
%             nSFCRparams.lambda2                 % for AutoRef
% 
% Updated 2021-07-17, X.L., fixed the bug with odd dimension data (for FANSI_L1 as in M-step)
% Updated 2021-07-31, X.L., added in TV option to use TV

    padsize = nSFCRparams.padsize;
    Params = nSFCRparams.Params;
    
    % padding
    N = size(deltaB); dimsOdd = mod(N, 2);      % N is original volume size
    deltaB = padarray(padarray(deltaB, dimsOdd, 'replicate', 'post'), padsize);
    maskErode = padarray(padarray(maskErode, dimsOdd, 'replicate', 'post'), padsize);
    DPWeight = padarray(padarray(DPWeight, dimsOdd, 'replicate', 'post'), padsize);
    nSFCRparams.maskRef = padarray(padarray(nSFCRparams.maskRef, dimsOdd, 'replicate', 'post'), padsize);
    
    Params.sizeVol = Params.sizeVol + dimsOdd + 2*padsize;
    Params.fov = Params.sizeVol.*Params.voxSize;

    % precalculation
    nSFCRparams.D = ifftshift(conv_kernel_rot_c0(Params, Params.TAng, nSFCRparams.datatype));
    nSFCRparams.mask = maskErode;
    phasescale = 2*pi*Params.gamma*Params.B0*mean(Params.TEs)*1e-6; % ppm to radian
    nSFCRparams.input = deltaB*phasescale;
    nSFCRparams.weight = DPWeight;    
    
    if nSFCRparams.TV == 0
        % for M-step
        if nSFCRparams.nlM == 0
            lambdaSet.beta = 10;
            lambdaSet.lambda1_M = 500;      
            lambdaSet.lambda2_M = 0;
            lambdaSet.datatype      = nSFCRparams.datatype;
            lambdaSet.gamma         = Params.gamma;         % constants
            lambdaSet.B0            = Params.B0;
            lambdaSet.TEs           = Params.TEs;
            lambdaSet.codeswitch    = 2;                     % adjust
            if Params.AutoRefFlag == 1
                lambdaSet.lambda2_M = lambdaSet.lambda1_M./5;
            end
        end

        if isfield(nSFCRparams, 'lambda_M')
            lambda1 = nSFCRparams.lambda_M;
        else
            if nSFCRparams.L1orL2 == 1
                lambda1 = 5e-3;     % AlphaSweep(ka); % for nlL1TV, better to do L-curve searching/tuning of alpha1 & mu1 & mu2=1
            elseif nSFCRparams.L1orL2 == 2
                lambda1 = 2.5e-4;   % AlphaSweep(ka); % for nlTV, default 2.5e-4, better to do L-curve searching/tuning of alpha1 & mu1 & mu2=1
            end                
        end
    
        nSFCRparams.lambda1 = lambda1;          % M_step gradient L1 penalty
        nSFCRparams.mu1 = 100*lambda1;          % gradient consistency weight (10-100)
        nSFCRparams.mu2 = 1;                    % data fidelity consistency weight
        nSFCRparams.mu3 = 1;
        nSFCRparams.mu_adap = 0;
        nSFCRparams.merit = 0;
        nSFCRparams.tol_update = 1;             % 1% with FANSI_L1, TV 
        nSFCRparams.maxOuterIter = 50;          %           
        MagEdgePer = 0.3;
        chiEdgePer = 0.3;
        wG1 = gradient_mask_all(nSFCRparams.weight, nSFCRparams.mask, MagEdgePer); % Mag

        disp('First step: SFCR-M: ')
        if nSFCRparams.nlM == 0
            % use old linear M-step
            nSFCRparams.kthresh = 0.18;              % tkd threshod
            [out_m.x, ~] = SFCR_m_SB(nSFCRparams.input, nSFCRparams.D, DPWeight, maskErode, ...
                nSFCRparams.maskRef, wG1, lambdaSet, 2, nSFCRparams.kthresh);
            wG3 = gradient_mask_all(out_m.x.*maskErode, nSFCRparams.mask, chiEdgePer); % chi_m, no straking, smoothed
            nSFCRparams.wG = wG3;

        else
            if nSFCRparams.L1orL2 == 1
                nSFCRparams.alpha1 = nSFCRparams.lambda1;       % using FANSI_L1_TV
                out_m = delta2chi_FANSI_nlL1TV(nSFCRparams);    % FANSI needs even dimen
            elseif nSFCRparams.L1orL2 == 2
                nSFCRparams.alpha1 = nSFCRparams.lambda1;
                out_m = delta2chi_FANSI_nlTV(nSFCRparams);    % with L2 data fidelity
            end
            wG2 = gradient_mask_all(delta2chi_tso(nSFCRparams.input, fftshift(nSFCRparams.D), 1).*maskErode, nSFCRparams.mask, chiEdgePer); % fastQSM, may have streaking 
            wG3 = gradient_mask_all(out_m.x.*maskErode, nSFCRparams.mask, chiEdgePer); % chi_m, no streaking, smoothed
            nSFCRparams.wG = (wG1 + wG2 + wG3) > 1;
            % nSFCRparams.x0 = out_m.x;
        end
        
    else
        disp('Skipping SFCR-M with TV.')
        nSFCRparams.wG = ones([size(deltaB), 3]);       % with TV for chi_m
        out_m.x = 0;
    end
    
    disp('Second step: SFCR-F: ') 
    if nSFCRparams.L1orL2 == 1
        if isfield(nSFCRparams, 'lambda_F')
            lambda1 = nSFCRparams.lambda_F;         % if defined outside
        else
            lambda1 = 2e-2;                         % default 2e-2
        end
        nSFCRparams.lambda1 = lambda1;
        nSFCRparams.mu1 = 100*lambda1;
        nSFCRparams.tol_update = 1;            % default 1
        nSFCRparams.maxOuterIter = 50;
        out_s = delta2chi_nSFCRL1(nSFCRparams); % with L1 data fidelity, sift data inconsistency
    elseif nSFCRparams.L1orL2 == 2
        out_s = delta2chi_nSFCRL2(nSFCRparams); % with L2 data fidelity
    end

    if nSFCRparams.nlM == 0
        chik_res = fftn(out_s.x);
        chik_M0 = fftn(out_m.x);
        inx1 = abs(nSFCRparams.D) <= nSFCRparams.kthresh - 0.05;
        inx2 = abs(nSFCRparams.D) > nSFCRparams.kthresh + 0.05;
        chik_res(~inx1 & ~inx2) = 0.5*(chik_M0(~inx1 & ~inx2) + chik_res(~inx1 & ~inx2));           
        chik_res(inx2) = chik_M0(inx2);
        chi_res = real(ifftn(chik_res)).*maskErode;
    else
        chi_res = out_s.x;
    end

    % automatic referencing 
    if Params.AutoRefFlag == 1
        chi_csf_ref = mean(chi_res(nSFCRparams.maskRef>0));
        chi_res = (chi_res - chi_csf_ref).*maskErode;
    end
    
    chi_res_M = (out_m.x/phasescale).*maskErode;
    chi_res = (chi_res/phasescale).*maskErode;
    chi_res_M = chi_res_M(padsize(1)+1:padsize(1)+N(1), padsize(2)+1:padsize(2)+N(2), padsize(3)+1:padsize(3)+N(3), :,:);
    chi_res = chi_res(padsize(1)+1:padsize(1)+N(1), padsize(2)+1:padsize(2)+N(2), padsize(3)+1:padsize(3)+N(3), :,:);
