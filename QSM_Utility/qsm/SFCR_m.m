function chi_res=SFCR_m(deltaB, D_k, maskErode, maskSS, wG, lambdaSet)
% chi_M = SFCR_m(deltaB, D, maskErode, maskSS, wG, lambdaSet);
% SFCR M-step 
% adopted from Lijun's Code
% used to get initial M-step estimate using the SFCR method
% using CS in k-space with L1 regularization with a priori info from Mag gradient
%
% maskSS is a mask for Strong chi source such as iron rich GM ()
%
% Ref: Bao et al, IEEE TMI, 2016 35(9):2040-50
% updated 2019-06

%% calculate chi with TKD
thred=0.2;  % result of 0.2 is better TKD threshold for getting initial estimate and CS thresthod

%% calculate chi_k with TKD
inx=abs(D_k)<thred;     % ill-conditioned region
D_kinv = 1./D_k;
D_kinv(inx)=sign(D_k(inx))*(1./thred);
chik_tkd=fftn(deltaB).*D_kinv;     % solution of TDK

clear D_kinv

%% D_k->H matrix in CS QSM model
H=~inx;                     % well-conditioned region
x = zeros(size(deltaB));    
beta = 4;

% %% ------------------------- adopted from Lijun's Original Code
TV = TVOP;                  % total varition class (3D finite difference operator)

tic
for titer = 1:2
    
    % geting gradient of chi_res: chi_res_grad
    % First get b (alpha here)    
    alpha  = SoftThresh(wG.*(TV*(x)),1/beta);

    % Second get chi_res (cg)
    right = lambdaSet.lambda1_M*ifftn(H.*(chik_tkd)) + beta*(TV'*(wG.*alpha));  % 
    left = @(z) (lambdaSet.lambda1_M*ifftn(H.*(fftn(z))) + lambdaSet.lambda2_M*(~maskSS).*z + beta*(TV'*(wG.*(TV*(z))))); 

    [xp,cgres] = cgsolve(left,right,1e-3,100,20);  %
    if (cgres >= 1/2)
        disp('error')
        return
    end

    beta = beta*2;
    x = xp;
    x = x.*maskErode; 

end
toc

chi_res = real(x);

end    
