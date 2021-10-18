function chi_S = SFCR_s_SB_k(deltaB, D_k, DPWeight, maskErode, maskSSV, wG, lambdaSet, x)
% chi_S = SFCR_s_SB_k(deltaB, D_k, DPWeight, maskErode, maskSSV, wG, lambdaSet, x)
% SFCR S-step using Split-Bregmen L1 solver, k-space version
% 
% Input: deltaB, D_k, DPWeight, maskErode, maskSSV, wG, lambdaSet, x
% 
% DPWeight is kind of SNR weighting
% maskSSV is a mask with possible artifacts, reserved for future use
%
% Ref: Chen et al, J Neurosci Res, 2019, 97:467:479
% updated 2019-06

datatype = class(deltaB);
if nargin < 8
    x = zeros(size(deltaB), datatype);        % initial value
end

if lambdaSet.lambda2_S ~= 0
    error('current version assumes lambda2_S to be zero')
end

verbose = 1;

N = size(deltaB);

% gradient operator in k-space
[k2,k1,k3] = meshgrid(0:N(2)-1, 0:N(1)-1, 0:N(3)-1);
k1 = cast(k1, datatype);
k2 = cast(k2, datatype);
k3 = cast(k3, datatype);

fd1 = -1 + exp(2*pi*1i*k1/N(1));    % dimention 1, delta(x+1) - delta(x), exp(j*2pi*k) - 1
fd2 = -1 + exp(2*pi*1i*k2/N(2));    % dimention 2
fd3 = -1 + exp(2*pi*1i*k3/N(3));    % dimention 3
E2 = abs(fd1).^2 + abs(fd2).^2 + abs(fd3).^2; 

clear k1 k2 k3
% E conjugate
cfd1 = conj(fd1);
cfd2 = conj(fd2);
cfd3 = conj(fd3);

beta = lambdaSet.beta;
betaVFlag = 1;   % beta variation Flag

% Adptive beta V2
betamu = 0.7;    % for adaptive beta
betatau = 2;
rkp1norm = Inf;    % inital residual norm

tol = 0.01;

titer = 0;
resnorm = 1;
titermax = 20;      % ususally converged in less than 20

dx = zeros(N, datatype);   
dy = zeros(N, datatype);  
dz = zeros(N, datatype);  

bx = zeros(N, datatype);  
by = zeros(N, datatype);
bz = zeros(N, datatype);

cgcount = 0;

tic
while (resnorm>tol) && (titer < titermax)
    
    right = lambdaSet.lambda1_S*D_k.*fftn(DPWeight.*deltaB) + beta*(cfd1.*fftn( (dx - bx).*wG(:,:,:,1) ) ...
                    + cfd2.*fftn( (dy - by).*wG(:,:,:,2) ) + cfd3.*fftn( (dz - bz).*wG(:,:,:,3) ));

    SB_precond = 1./(eps + lambdaSet.lambda1_S*D_k.^2 + beta*E2);
                
    [xp,cpflag,cgrelres, cgiter] = pcg(@afun,right(:),1e-2, 100, @mfun, [], x(:));
    % disp(['step cg iteration: ',  num2str(cgiter)])
    
    cgcount = cgcount + cgiter;    
    if (cpflag) || (cgrelres >= 1/2)
        disp(['pcg  not converging, with pcgflag: ', num2str(cpflag)])
    end    
    
    % 2nd get d_k+1 using soft thresholding (d)    
    xp = reshape(xp, N);

    phix = ifftn(fd1.*xp).*wG(:,:,:,1);
    phiy = ifftn(fd2.*xp).*wG(:,:,:,2);
    phiz = ifftn(fd3.*xp).*wG(:,:,:,3);          
    
    % testing    
    dxp  = SoftThresh(phix+bx,1/beta);
    dyp  = SoftThresh(phiy+by,1/beta);
    dzp  = SoftThresh(phiz+bz,1/beta);
    
    % 3rd update residual b_k+1 (scaled dual variable)
    bx = bx + (phix - dxp);          
    by = by + (phiy - dyp);
    bz = bz + (phiz - dzp);    
        
    % estimating residual primary and dual V2
    skp1norm = rkp1norm;             % old primary residual
    rkp1norm = norm([(phix(:)-dxp(:)); (phiy(:) - dyp(:)); (phiz(:) - dzp(:))]);   % primary residual r_k+1
    
    % check convergence
    resnorm = norm(x(:) - xp(:))./norm(xp(:));
    
    if verbose
        disp(['primary residual:   ', num2str(rkp1norm)]);   
        disp(['SB relative change: ', num2str(resnorm)]);
    end
    
    if (resnorm < tol)
        break;
    end
         
    % update iteration 
    x = xp;        
    dx = dxp;
    dy = dyp;
    dz = dzp;
    titer = titer + 1;
    
    % update adaptive beta
    if betaVFlag == 1
        if (rkp1norm > betamu*skp1norm) 
            beta = beta*betatau;            % beta increase
            bx = bx./betatau;               % change dual variable accordingly
            by = by./betatau;
            bz = bz./betatau;
        end

    end    
end
toc

disp(['total cg iteration: ',  num2str(cgcount)])
chi_S = real(ifftn(xp));               % Note, k-space version solves F(chi_k+1)

function y = afun(x)    
    x = reshape(x, N);
    y = lambdaSet.lambda1_S*D_k.*fftn(DPWeight.*ifftn(D_k.*x))  + beta*(fftn(wG(:,:,:,1).*ifftn(fd1.*x)).*cfd1 + ...
            fftn(wG(:,:,:,2).*ifftn(fd2.*x)).*cfd2 + ...
            fftn(wG(:,:,:,3).*ifftn(fd3.*x)).*cfd3); 
    y = y(:);
end

function y = mfun(x)       
   y = SB_precond(:).*x(:);
end

end  