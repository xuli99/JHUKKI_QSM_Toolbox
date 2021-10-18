function [chi_res, SB_residual, SB_regularization, chi_M] = delta2chi_SFCR(deltaB, D, maskErode, DPWeight, MagEdgePer, chiEdgePer, lambdaSet, thred, AutoRefFlag)
% [chi_res, SB_residual, SB_regularization, chi_M] = delta2chi_SFCR(deltaB, D, maskErode, DPWeight, MagEdgePer, chiEdgePer, lambdaSet, thred, AutoRefFlag)
% delta2chi_SFCR
% Ref: Bao et al, IEEE TMI, 2016 35(9):2040-50
%      Chen et al, J Neurosci Res, 2019, 97:467:479
% Authors: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
%
% INPUT:
% deltaB: relative field shift
% D: convolution kernel
% maskErode: brain mask
% DPWeight: normalized magnitude data, SNR weighting
% MagEdgePer: edge percentage in Magnitude data
% chiEdgePer: edege percentage in susceptibilty data (intermediate estimate)
% lambdaSet: L1 regularization parameter in the SFCR formulation, including
%         lambdaSet.beta
%         lambdaSet.lambda1_M
%         lambdaSet.lambda2_M
%         lambdaSet.lambda1_S 
%         lambdaSet.lambda2_S
%         lambdaSet.maskSS
% thred: k-space threshold for TKD and k-space combination
% AutoRefFlag:    Flag for Auto Reference to maskSS
% chi_M:  optional input for the M-step

% Update Note: 2016-07-19, Added code switch for k-space based SB solver
% update Note: 2018-03-09, X.L., added in SFCR+0 based on R2* map with
%        automatic referencing, use code switch 2
% updated 2019-07-09

% updated 2020-06-17, added option to use NDI for M-step

if nargin < 8
    thred = 0.2;        % default
    AutoRefFlag = 0;
elseif nargin < 9
    AutoRefFlag = 0;
end

% code switch
% 1 for using original code
% 2 for using modifed SFCR with new SB solver or SFCR+0 if doing auto ref.
% 3 for using SFCR with k-space based solver

% precalculations:
TV = TVOP; 

% input for SFCR_m
if AutoRefFlag == 1 && isfield(lambdaSet, 'maskSS')
    maskSS = lambdaSet.maskSS;  
    maskSSV = maskSS;
    codeswitch = lambdaSet.codeswitch;  
else
    maskSS = maskErode;    % 
    maskSSV = maskSS;      % Saved for future use  
    codeswitch = lambdaSet.codeswitch;        % SFCR k-space version, faster
end

if codeswitch == 3
    wG = gradient_mask_all(DPWeight, maskErode, MagEdgePer, 1); 
else
    wG = gradient_mask_all(DPWeight, maskErode, MagEdgePer);
end

%% ---------------------
% M-Step        
disp('First step: SFCR-M: ')   

switch codeswitch
    case 1
        chi_M = SFCR_m(deltaB, D, maskErode, maskSS, wG, lambdaSet);    %  Original SFCR code
    case {0;2}
       [chi_M0, ~] = SFCR_m_SB(deltaB, D, DPWeight, maskErode, maskSS, wG, lambdaSet, 2, thred);   % SFCR+0
       chi_M = chi_M0.*maskErode;
    case 3
       [chi_M, chik_tkd] = SFCR_m_SB_k(deltaB, D, maskErode, maskSS, wG, lambdaSet, 2, thred);  % k-space version of SFCR_m_SB      
    otherwise
        return;
end

% S-Step
disp('Seond step: SFCR-S: ')       
% extract a priori from chi_M
if codeswitch == 3
    wG = gradient_mask_all(chi_M, maskErode, chiEdgePer, 1);    % using chi_M to get a priori, use k-space based, faster  
elseif codeswitch == 2
    wG = gradient_mask_all(chi_M, maskErode, chiEdgePer);          % using chi_M to get a priori, use image based
elseif codeswitch == 0
    % use the same wG, % good if with long TR
end

switch codeswitch
    case 1
        chi_res = SFCR_s(deltaB, D, DPWeight, maskErode, maskSSV, wG, lambdaSet);     %  Oriignal SFCR code    
    case {0; 2}        
        chi_res = SFCR_s_SB(deltaB, D, DPWeight, maskErode, maskSSV, wG, lambdaSet);  % SFCR+0  
    case 3        
        chi_res = SFCR_s_SB_k(deltaB, D, DPWeight, maskErode, maskSSV, wG, lambdaSet);  % S fitting in image space, k-space version
    otherwise
        return;
end

switch codeswitch
    case {0; 2}
        chik_res = fftn(chi_res);
        chik_M0 = fftn(chi_M0);
%         inx = abs(D) <= thred;
%         chik_res(~inx) = chik_M0(~inx);           
        inx1 = abs(D) <= thred - 0.05;
        inx2 = abs(D) > thred + 0.05;
        chik_res(~inx1 & ~inx2) = 0.5*(chik_M0(~inx1 & ~inx2) + chik_res(~inx1 & ~inx2));           
        chik_res(inx2) = chik_M0(inx2);

        chi_res = real(ifftn(chik_res)).*maskErode;
        chi_res = chi_res.*maskErode;
        
        % automatic referencing 
        chi_csf_ref = mean(chi_res(maskSS>0));
        chi_res = (chi_res - chi_csf_ref).*maskErode;
        
    case 3
        inx = abs(D) <= thred;
        chik_res = fftn(chi_res);
        chik_res(~inx) = chik_tkd(~inx);  
        chi_res = real(ifftn(chik_res)).*maskErode;
        chi_res = chi_res.*maskErode;
    otherwise
        % do nothing
end

% estimate residuals of data fidelity and regularization terms 
SB_residual = DPWeight.*(real(ifftn(fftn(chi_res).*D)).*maskErode - deltaB);
SB_residual = norm(SB_residual(:));

SB_regularization = TV*(chi_res).*wG;       % L1 norm
SB_regularization = norm(SB_regularization(:), 1);

fprintf('Final Redidual: %4.2g;  L1 norm: %4.2g \n', SB_residual, SB_regularization)

end
