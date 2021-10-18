function [chi_M, chik_tkd] =SFCR_m_SB_k(deltaB, D_k, maskErode, maskSS, wG, lambdaSet, itermax, thred)
% [chi_M, chik_tkd] =SFCR_m_SB_k(deltaB, D_k, maskErode, maskSS, wG, lambdaSet, itermax, thred)
% SFCR M-step using Split-Bregmen L1 solver, usin k-space based version
%
% Input: deltaB, D_k, maskErode, maskSS, wG, lambdaSet, itermax, thred
% updated 2019-06
%
% Ref: Chen et al, J Neurosci Res, 2019, 97:467:479

if nargin < 7    
    itermax = 20;
    thred = 0.2;
elseif nargin < 8
    thred = 0.2;
end

verbose = 1;

if lambdaSet.lambda2_M ~= 0
    error('current version assumes lambda2_M to be zero')
end

N = size(deltaB);

datatype = class(deltaB);

%% calculate chi_k with TKD
inx = abs(D_k) <= thred;     % ill-conditioned region
D_kinv = 1./D_k;
D_kinv(inx)=sign(D_k(inx));        % TKD v2
chik_tkd=fftn(deltaB).*D_kinv;     % solution of TKD in k space
clear D_kinv

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

%% D_k->H matrix in CS QSM model
H=~inx;                     % well-conditioned region

beta = lambdaSet.beta;      % initial beta
betaVFlag = 1;   % beta variation Flag

% Adptive beta V2
betamu = 0.7;    % for adaptive beta
betatau = 2;
rkp1norm = Inf;    % inital residual norm

tol = 0.01;

titer = 0;
resnormx = 1;                

% Hard start
x = zeros(size(deltaB), datatype);    % initial value

dx = zeros(N, datatype);              
dy = zeros(N, datatype);  
dz = zeros(N, datatype);  

bx = zeros(N, datatype);              
by = zeros(N, datatype);
bz = zeros(N, datatype);

cgcount = 0;

tic
while (titer < itermax) && (resnormx > tol)
    
    % 1st using cg to get chi_k+1 (xp)
    right = lambdaSet.lambda1_M*H.*chik_tkd + beta*(cfd1.*fftn( (dx - bx).*wG(:,:,:,1) ) ...
                    + cfd2.*fftn( (dy - by).*wG(:,:,:,2) ) + cfd3.*fftn( (dz - bz).*wG(:,:,:,3) ));
    
    SB_precond = 1./(eps + lambdaSet.lambda1_M*H + beta*E2);
                
    [xp,cpflag,cgrelres, cgiter] = pcg(@afun,right(:),1e-3, 100, @mfun, [], x(:));
    
    cgcount = cgcount + cgiter;
    
    if (cpflag) || (cgrelres >= 1/2)
        disp(['error: pcg  not converging, with pcgflag: ', num2str(cpflag)])
        return
    end    
    
    % 2nd get d_k+1 using soft thresholding (d)
    xp = reshape(xp, N);
    phix = ifftn(fd1.*xp).*wG(:,:,:,1);
    phiy = ifftn(fd2.*xp).*wG(:,:,:,2);
    phiz = ifftn(fd3.*xp).*wG(:,:,:,3);          
    
    dx  = SoftThresh(phix+bx,1/beta);
    dy  = SoftThresh(phiy+by,1/beta);
    dz  = SoftThresh(phiz+bz,1/beta);
    
    % 3rd update residual b_k+1
    bx = bx + (phix - dx);
    by = by + (phiy - dy);
    bz = bz + (phiz - dz);

    % estimating residual primary and dual V2
    skp1norm = rkp1norm;             % old primary residual
    rkp1norm = norm([(phix(:)-dx(:)); (phiy(:) - dy(:)); (phiz(:) - dz(:))]);   % primary residual r_k+1
    
    % check convergence
    resnormx = norm(x(:) - xp(:), 2)./norm(xp(:), 2);

    if verbose
        disp(['primary residual:   ', num2str(rkp1norm)]);        
        disp(['SB relative change in x: ', num2str(resnormx)]);
    end
    
    if (resnormx < tol)
        break;
    end
    
    % update iteration 
    x = xp;       
    titer = titer + 1;
    
    % update adaptive beta
    if betaVFlag == 1
        if (rkp1norm > betamu*skp1norm) 
            beta = beta*betatau;            % beta increase
            bx = bx./betatau;              
            by = by./betatau;
            bz = bz./betatau;
        end
    end        
    
end
toc

disp(['total cg iteration: ',  num2str(cgcount)]);
chi_M = real(ifftn(xp)).*maskErode;         % Note, k-space version solves F(chi_k+1)

function y = afun(x)    
    x = reshape(x, N);
    y = lambdaSet.lambda1_M*H.*x  + beta*(fftn(wG(:,:,:,1).*ifftn(fd1.*x)).*cfd1 + ...
            fftn(wG(:,:,:,2).*ifftn(fd2.*x)).*cfd2 + ...
            fftn(wG(:,:,:,3).*ifftn(fd3.*x)).*cfd3); 
    y = y(:);
end

function y = mfun(x)    
   y = SB_precond(:).*x(:);
end

end    
