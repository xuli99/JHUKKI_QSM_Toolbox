function [chi_res] = delta2chi_iTKD(deltaB, D, maskErode, tol)     
% function [chi_res] = delta2chi_iTKD(deltaB, D, maskErode, tol)     
% Input:
% deltaB: relative field shift
% D: convolution kernel
% maskErode: brain mask
% tol: converge tolerance
%
% Authors: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
%
% Ref: Aggarwal, Li et. al. JMRI, 2018, 47(2):554-564
%
% updated 2019-07-09

D = ifftshift(D);

% filtering window Wfs    
a = prctile(abs(D(:)), [1, 30]);
W0 = (abs(D) - a(1))./(a(2)-a(1));            
Wfs = W0;
Wfs(W0<0) = 0;              % k-space regions where averaging filter takes place
Wfs(W0>1) = 1;

% modified Dinv
Dinv = sign(D);             % TKD, with threshold of 1

psf = real(ifftn(D.*Dinv));
c = psf(1);

% step 1, tkd
chi_res1k = Dinv.*(fftn(deltaB));       % fitting what is left       
chi_res1k = chi_res1k./c;               % psf correction

% step 2, filter
chi_res1k_filt = smvfilter(chi_res1k, 3);           % smv filter
chi_res1 = real(ifftn(((chi_res1k.*Wfs + chi_res1k_filt.*(1-Wfs)))));

% step 3, mask and filter
chi_res1 = chi_res1.*maskErode;  
chi_res1k = (fftn(chi_res1));

chi_res1k_filt = smvfilter(chi_res1k, 3);           % smv filter
chi_res1 = real(ifftn(((chi_res1k.*Wfs + chi_res1k_filt.*(1-Wfs)))));         

chi_res1 = chi_res1.*maskErode;                   
deltaBSim = real(ifftn(fftn(chi_res1).*D)).*maskErode;
chi_res = chi_res1;                             % initial fast QSM

% iteration setup
resTol = tol;                                   % fitting residual
solTol = tol;
maxit = 50;

disp_flag = 1;

deltaBnorm = norm(deltaB(:));
ResINdeltaB = norm(deltaBSim(:) - deltaB(:))./deltaBnorm;   % residual 
ReInSol = 1;
iter = 0;

while ((ResINdeltaB > resTol) && ReInSol > solTol && (iter < maxit))
            
    % step 1, tkd
    chi_res1k = Dinv.*(fftn(deltaB - deltaBSim));    % fitting what is left       
    
    % step 2, filter
    chi_res1k_filt = smvfilter(chi_res1k, 3);           % smv filter
    chi_res1 = real(ifftn(((chi_res1k.*Wfs + chi_res1k_filt.*(1-Wfs)))));

    % step 3, mask and filter
    chi_res1 = chi_res1.*maskErode;  
    chi_res1k = (fftn(chi_res1));

    chi_res1k_filt = smvfilter(chi_res1k, 3);           % smv filter
    chi_res1 = real(ifftn(((chi_res1k.*Wfs + chi_res1k_filt.*(1-Wfs)))));         

    chi_res1 = chi_res1.*maskErode;
    
    chi_res = chi_res + chi_res1;                      
    
    ReInSol = norm(chi_res1(:))./(norm(chi_res(:)));   % Updata ReInSol
    
    deltaBSim = real(ifftn(fftn(chi_res).*D)).*maskErode;  
    ResINdeltaB = norm(deltaBSim(:) - deltaB(:))./deltaBnorm;  
            
    iter = iter + 1;    
    if disp_flag == 1
        disp(['iteration: ', num2str(iter)]);
        disp(['fitting relative residual: ', num2str(ResINdeltaB)]);
        disp(['relative update in solution: ', num2str(ReInSol)]);        
    end
    
end
