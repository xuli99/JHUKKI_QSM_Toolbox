function out = delta2chi_nSFCRL1(params)
% out = delta2chi_nSFCRL1(params)
% delta2chi_nSFCRL1
% Ref: Bao et al, IEEE TMI, 2016 35(9):2040-50
%      Chen et al, J Neurosci Res, 2019, 97:467:479
%      Milovic, MRM, 2018, 80:814–821 (FANSI)
%      Milovic, 2019, QSMworkshop 
%      FANSI toolbox:  https://gitlab.com/cmilovic/FANSI-toolbox 
%
% Authors: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
%
% INPUT:
% params.input:         local phase
% params.weight:        SNR weight
% params.D:             kernel
% params.wG:            gradient weighting a prori
% params.lambda1:       TV or gradient L1 penalty
% params.mu1:           gradient consistency weight (10-100)
% params.mu2:           data fidelity consistency weight
% params.mu3:           L1 data fidelity consistency weight
% params.mask:          brain mask
% params.maskRef:       mask for AutoReference
% params.lambda2:       AutoRef L2 penalty
% params.maxOuterIter:  maximum iteration
% 
% Output: out - structure with the following fields:
% out.x: calculated susceptibility map
% out.iter: number of iterations needed
% out.time: total elapsed time (including pre-calculations)
% 
% Adapted from code from FANSI toolbox: https://gitlab.com/cmilovic/FANSI-toolbox
%
% Updated 2021-07-21, added warm start option

tic
    % main params
    lambda1 = params.lambda1;
    lambda2 = params.lambda2;
    D_k = params.D;
    phase = params.input;
    wG = params.wG;             % binary
    mask = params.mask;
    maskRef = params.maskRef;
    
    N = size(phase);
    NwG = size(wG);
    datatype = class(phase);

    % functions
    TV = TVOP;                   
    A = @(z)real(ifftn(D_k.*fftn(z)));   
    AA = @(z)real(ifftn(D_k.^2.*fftn(z)));    
    
    % optional parameters
    if isfield(params,'maxOuterIter')
        itermax = params.maxOuterIter;
    else
        itermax = 50;
    end

    if isfield(params,'tol_update') % for outer iteration
       tol_update  = params.tol_update;
    else
       tol_update = 1;
    end

    if ~isfield(params,'delta_tol')
        delta_tol = 1e-6;
    else
        delta_tol = params.delta_tol;
    end

    if isfield(params,'mu1')
        mu1 = params.mu1;
    else
        mu1 = 100*lambda1;
    end
    
    if isfield(params,'mu2')
        mu2 = params.mu2;
    else
        mu2 = 1.0;
    end
    
    if isfield(params,'mu3')
        mu3 = params.mu3;
    else
        mu3 = 1.0;
    end
    
    if isfield(params,'precond')
        precond = params.precond;
    else
        precond = 1;
    end
    
    if isfield(params,'merit')
        merit = params.merit;
    else
        merit = 0;
    end
    
    if isfield(params,'mu_adap')
        mu_adap = params.mu_adap;
        mu1lambda = 0.75;  
        mu1gamma = 2;  
    else
        mu_adap = 0;
    end
    
    % initialization
    rkp1norm = Inf;  % r_k+1_norm
            
    IS = exp(1i*phase);    
    if isfield(params,'weight')
       W  = params.weight;
    else
       W = ones(N);
    end    
    
    if isfield(params, 'x0')
        x = params.x0;
        z1 = wG.*(TV*(x)) + eps;
    else
        x = zeros(N, datatype);  
        z1 = zeros(NwG, datatype);  
    end
    s1 = zeros(NwG, datatype);  

    if precond
        z2 =  W.*phase./(W+mu2); % start with something similar to the input phase, weighted to reduce noise
    else
        z2 = zeros(N, datatype);
    end
    s2 = zeros(N, datatype); 
    
    z3 = zeros(N,'single');
    s3 = zeros(N,'single');
    
    cM = 1./sum(maskRef(:) > 0);   
    % oneT = ones(1, prod(N), datatype); % 1XN
    
    cgcount = 0;
    titer = 0;
    x_update = 100;     % percent

    lm = lambda1/mu1;
    debug = 0;
    if debug == 1
        verbose = 1;
        cgverbose = 1;
    else
        verbose = 0;
        cgverbose = 0;
    end
    
    % Main iteration
    while (x_update>tol_update) && (titer < itermax)
        titer = titer + 1;
        
        % sp 1
        right = mu2*A(z2 - s2) + mu1*(TV'*(wG.*(z1 - s1)));      
        [xp,cgrelres, cgiter] = cgsolve(@leftFun,right,1e-2,100, cgverbose, x);      % cgsolve
        cgcount = cgcount + cgiter;
        if (cgrelres >= 1/2)
            disp(' cg not converging ...')
        end
        
        xp = reshape(xp, N);
        phi = wG.*(TV*(xp));

        z1  = SoftThresh(phi+s1, lm);
        s1 = s1 + (phi - z1);

        x_update = 100* norm(x(:) - xp(:))./norm(xp(:));
        disp(['Iter: ', num2str(titer), '   Update: ', num2str(x_update)])
        if (x_update < tol_update) || isnan(x_update)
            break;
        end

        skp1norm = rkp1norm;           
        rkp1norm = norm(phi(:)-z1(:));    

        if verbose
            disp(['primary residual:   ', num2str(rkp1norm)]);
        end

        % update z3 first
        z3  = SoftThresh(exp(1i*z2)-IS+s3, W./(mu3+eps));

        % sp 2
        Dxp = A(xp);
        rhs_z2 = mu2*(Dxp + s2);
        z2 =  rhs_z2 ./ mu2 ;   % z2 with xp

        % Newton-Raphson method
        delta = inf;
        inn = 0;        
        yphase = angle( IS+z3-s3 );
        ym = abs(IS+z3-s3);
        
        while (delta > delta_tol && inn < 50)
            inn = inn + 1;
            norm_old = norm(z2(:));
            
            update = ( mu3 .* sin(z2 - yphase-1i*log(ym)) + mu2*z2 - rhs_z2 )./( mu3 .* cos(z2 - yphase-1i*log(ym)) + mu2 +eps);   
            % update = ( mu3 .* ym.*sin(z2 - yphase) + mu2*z2 - rhs_z2 )./( mu3 .* ym.*cos(z2 - yphase) + mu2 +eps);   
                  
            z2 = z2 - update;     
            delta = norm(update(:)) / norm_old;
        end        
        if verbose
            disp(delta)
        end
        s2 = s2 + Dxp - z2;

        % sp 3
        s3 = s3 + exp(1i*z2)-IS - z3;
        
        x = xp;     
        
        if mu_adap == 1
            % mu1, s1
            if (rkp1norm > mu1lambda*skp1norm) && (mu1 < 1)
                disp('update mu1 to')
                mu1 = mu1*mu1gamma;
                disp(mu1)
                s1 = s1./mu1gamma;                 
            end
            % mu2, s2
        end  

        if merit
            wres = params.weight.*(exp(1i*A(x)) - exp(1i*phase));
            wres = wres - mean(wres(mask>0));       % normalized residual
            a = wres(mask>0);
            factor = std(abs(a))*6;                 % std of abs(res) * 6 sigma
            wres = abs(wres)/factor;                % divided 6 sigma
            wres(wres<1) = 1;                       % within 6 sigma ==> 1, does not change
            wres_t = 2;
            W = (params.weight./(wres.^wres_t)).^2;
        end

    end
    
out.time = toc;toc
out.x = real(xp).*mask;
out.iter = titer;

chi_ref = mean(out.x(maskRef>0));
out.x = (out.x - chi_ref).*mask;

    function y = leftFun(x)    
        DfRgterm = mu2*AA(x) + mu1*(TV'*(wG.*(TV*(x))));
        temp = maskRef.*(x - cM*maskRef(:)'*x(:));
        L2term = lambda2*(temp - maskRef.*cM*(sum(temp(:))));
        y = DfRgterm + L2term;
    end
end
