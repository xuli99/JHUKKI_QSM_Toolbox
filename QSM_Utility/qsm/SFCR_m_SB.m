function [chi_M, chik_tkd] =SFCR_m_SB(deltaB, D_k, DPWeight, maskErode, maskSS, wG, lambdaSet, itermax, thred)
% [chi_M, chik_tkd] =SFCR_m_SB(deltaB, D_k, maskErode, maskSS, wG, lambdaSet, itermax, thred)
% SFCR M-step using Split-Bregmen L1 solver
%
% Input: deltaB, D_k, maskErode, maskSS, wG, lambdaSet, itermax, thred
% thred is the threshold used in TKD
%
% Authors: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
% 
% updated 2018-03-07, X.L, changed the maskSS to be CSF mask similar like
% MEDI+0 option giving SFCR+0 M-step solution
% updated 2019-06

if nargin < 7    
    itermax = 20;
    thred = 0.2;
elseif nargin < 8
    thred = 0.2;
end

N = size(deltaB);
datatype = class(deltaB);
if isa(deltaB, 'single')
    verbose = 1;
    cgverbose = 50;
else
    verbose = 1;
    cgverbose = 0;
end

NwG = size(wG);
cM = 1./sum(maskSS(:) > 0);   
oneT = ones(1, prod(N), datatype);

%% calculate chi_k with TKD or NDI

codeswitch = lambdaSet.codeswitch;   % 2 for TKD, 0 for NDI

inx = abs(D_k) <= thred;    

if codeswitch == 2
    D_kinv = 1./D_k;
    D_kinv(inx)=sign(D_k(inx));        
    chik_tkd=fftn(deltaB).*D_kinv;     
    clear D_kinv
else
    NDIParams.step_size = 0.5;       % gradient descent step size
    NDIParams.num_iter = 100; 
    NDIParams.tik  = 1;            % adding tikhonov regularization
    NDIParams.lambda = 1e-3;        % 0.1%
    NDIParams.tol = 0.5;
    NDIParams.datatype = lambdaSet.datatype;
    NDIParams.gamma = lambdaSet.gamma;          % constants
    NDIParams.B0 = lambdaSet.B0;
    NDIParams.TEs = lambdaSet.TEs;
    NDIParams.verbose = 0;
    [chi_M, ~, ~, ~] = delta2chi_NDI(deltaB, D_k, DPWeight, maskErode, NDIParams);  
    chik_tkd = fftn(chi_M);
end


%% D_k->H matrix in CS QSM model
H=~inx;                     
n1 = abs(chik_tkd(H));

TV = TVOP;   

beta = lambdaSet.beta;     
betaVFlag = 1;  

betamu = 0.7;    
betatau = 2;
rkp1norm = Inf;    

tol = 0.1;

titer = 0;
resnormx = 1;  

x = zeros(size(deltaB), datatype);  
d = zeros(NwG, datatype);             
b = zeros(NwG, datatype);             

cgcount = 0;

tic
while (titer < itermax) && (resnormx > tol)
    
    right = lambdaSet.lambda1_M*ifftn(H.*(chik_tkd)) + beta*(TV'*(wG.*(d - b)));  %
    [xp,cgrelres, cgiter] = cgsolve(@leftFun,right,1e-3,100, cgverbose, x);  % cgsolve    
    cgcount = cgcount + cgiter;
    
    if (cgrelres >= 1/2)
        disp('error: cg not converging ...')
        return
    end    
    
    xp = reshape(xp, N);
    phi = wG.*(TV*(xp));
    d  = SoftThresh(phi+b,1/beta);
    
    b = b + (phi - d);

    skp1norm = rkp1norm;            
    rkp1norm = norm(phi(:) - d(:)); 
    
    resnormx = norm(x(:) - xp(:), 2)./norm(xp(:), 2);
    
    if verbose
        disp(['primary residual:   ', num2str(rkp1norm)]);            
        disp(['SB relative change in x: ', num2str(resnormx)]);    
        e1 = abs((fftn(xp) - chik_tkd)); e1 = e1(H); 
        disp(['fidelity error: ', num2str(norm(e1(:))./norm(n1(:)))])
    end
    
    if (resnormx < tol)
        break;
    end
    
    x = xp;       
    titer = titer + 1;
    
    if betaVFlag == 1
        if (rkp1norm > betamu*skp1norm) 
            beta = beta*betatau;    
            b = b./betatau;               
        end
    end        
    
end
toc

disp(['total cg iteration: ',  num2str(cgcount)])    

chi_M = real(xp);

function y = leftFun(x)    
        DfRgterm = lambdaSet.lambda1_M*ifftn(H.*fftn(x)) + beta*(TV'*(wG.*(TV*(x))));
        temp = maskSS.*(x - cM*maskSS(:)'*x(:));
        L2term = lambdaSet.lambda2_M*(temp - maskSS.*cM*(oneT*temp(:)));
        y = DfRgterm + L2term;
end


end    
