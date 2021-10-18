function chi_S = SFCR_s_SB(deltaB, D_k, DPWeight, maskErode, maskSSV, wG, lambdaSet, titermax)
% chi_S = SFCR_s_SB(deltaB, D_k, DPWeight, maskErode, maskSSV, wG, lambdaSet, titermax)
% SFCR S-step using Split-Bregmen L1 solver
% 
% Input: deltaB, D_k, DPWeight, maskErode, maskSSV, wG, lambdaSet, x
% DPWeight is kind of SNR weighting
%
% Authors: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
%
% updated 2018-03-07, X.L, changed the maskSSV to be CSF mask similar like
% MEDI+0 option, giving SFCR+0 S-step solution
% 2019-05-02: added choice of merit

% cleanned, 2019-06
% updated, 2020-06-17

datatype = class(deltaB);
if isa(deltaB, 'single')
    verbose = 0;
    cgverbose = 50;
else
    verbose = 1;
    cgverbose = 0;
end

if nargin < 8
    titermax = 20; 
end

x = zeros(size(deltaB), datatype);  

N = size(deltaB);
NwG = size(wG);
cM = 1./sum(maskSSV(:) > 0);   
oneT = ones(1, prod(N), datatype);

beta = lambdaSet.beta;
betaVFlag = 1;  
merit = 1;     

betamu = 0.7;  
betatau = 2;
rkp1norm = Inf;  

tol = 0.01;  

titer = 0;
resnorm = 1;

d = zeros(NwG, datatype);     
b = zeros(NwG, datatype);     

TV = TVOP;                   
A = @(z)ifftn(D_k.*fftn(z));   

cgcount = 0;

tic
while (resnorm>tol) && (titer < titermax)
    
    right = lambdaSet.lambda1_S*A(DPWeight.*deltaB) + beta*(TV'*(wG.*(d - b)));      
    [xp,cgrelres, cgiter] = cgsolve(@leftFun,right,1e-2,100, cgverbose, x);      % cgsolve
    % disp(['step cg iteration: ',  num2str(cgiter)])
    
    cgcount = cgcount + cgiter;
    
    if (cgrelres >= 1/2)
        disp(' cg not converging ...')
    end

    xp = reshape(xp, N);
    phi = wG.*(TV*(xp));
        
    dp  = SoftThresh(phi+b,1/beta);
    
    b = b + (phi - dp);
    
    resnorm = norm(x(:) - xp(:))./norm(xp(:));    
    if (resnorm < tol)
        break;
    end
    
    skp1norm = rkp1norm;           
    rkp1norm = norm(phi(:)-dp(:));    
    
    if verbose
        disp(['primary residual:   ', num2str(rkp1norm)]);
        disp(['SB relative change: ', num2str(resnorm)]);
    end
    
    x = xp;     
    d = dp;
    titer = titer + 1;
 
    if betaVFlag == 1
        if (rkp1norm > betamu*skp1norm) 
            beta = beta*betatau;           
            b = b./betatau;                 
        end

    end  
    
    if merit
        wres = DPWeight.*(exp(1i*real(A(x))) - exp(1i*deltaB));
        % wres = DPWeight.*(A(x)-deltaB);        
        wres = wres - mean(wres(maskErode>0));  % normalized residual
        a = wres(maskErode>0);
        factor = std(abs(a))*6;                 % std of abs(res) * 6 sigma
        wres = abs(wres)/factor;                % divided 6 sigma
        wres(wres<1) = 1;                       % within 6 sigma ==> 1, does not change
        t = 2;
        DPWeight= DPWeight./wres.^t;
    end
    
end

toc

disp(['total cg iteration: ',  num2str(cgcount)])
chi_S = real(xp);

function y = leftFun(x)    
        DfRgterm = lambdaSet.lambda1_S*A(DPWeight.*A(x)) + beta*(TV'*(wG.*(TV*(x))));
        temp = maskSSV.*(x - cM*maskSSV(:)'*x(:));
        L2term = lambdaSet.lambda2_S*(temp - maskSSV.*cM*(oneT*temp(:)));
        y = DfRgterm + L2term;
end

end  