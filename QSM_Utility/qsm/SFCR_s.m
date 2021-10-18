function chi_res=SFCR_s(deltaB, D_k, DPWeight, maskErode, maskSSV, wG, lambdaSet, x)
% chi_S = SFCR_s(deltaB, D, DPWeight, maskErode, maskU, wG, lambdaSet);
% SFCR S-step 
% adopted from Lijun's Code
% Ref: Bao et al, IEEE TMI, 2016 35(9):2040-50
% updated 2019-06


if nargin < 8
    x = zeros(size(deltaB));    
end
beta = 4;

% %% ------------------------- adopted from Lijun's Original Code
TV = TVOP;                  % total varition class (3D finite difference operator)
A = @(z)ifftn(D_k.*fftn(z));  %

titermax = 4;
titer = 1;
resnorm = 1;
tol = 0.2;

tic
while (resnorm>tol) && (titer <= titermax)
    

    alpha  = SoftThresh(wG.*(TV*(x)),1/beta);

    right = lambdaSet.lambda1_S*A(DPWeight.*deltaB) + beta*(TV'*(wG.*alpha));  % 
    left = @(z) (lambdaSet.lambda1_S*A(DPWeight.*A(z))+ beta*(TV'*(wG.*(TV*(z))))) + lambdaSet.lambda2_S*(maskSSV).*z ;     
    
    [xp,cgres] = cgsolve(left,right,1e-2,100,20);  %    
    
    if (cgres >= 1/2)
        disp(' cg not converging ...')
        % return
    end
    
    resnorm = norm(x(:) - xp(:))./norm(xp(:));
    fprintf('SFCR_s relative change: %4.2g \n', resnorm);

    beta = beta*2;
    
    x = xp;
    x = x.*maskErode; 
        
    titer = titer + 1;
end

chi_res = real(x);

end  