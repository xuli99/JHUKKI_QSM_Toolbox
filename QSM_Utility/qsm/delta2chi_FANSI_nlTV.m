function out = delta2chi_FANSI_nlTV(params)
% Nonlinear QSM and Total Variation regularization 
% with spatially variable fidelity and regularization weights.
% This uses ADMM to solve the functional.
%
% Parameters: params - structure with 
% Required fields:
% params.input: local field map
% params.D: dipole kernel in the frequency space
% params.alpha1: gradient penalty (L1-norm) or regularization weight
% Optional fields:
% params.mu1: gradient consistency weight (ADMM weight, recommended = 100*alpha1)
% params.mu2: fidelity consistency weight (ADMM weight, recommended value = 1.0)
% params.maxOuterIter: maximum number of iterations (recommended = 150)
% params.tol_update: convergence limit, update ratio of the solution (recommended = 0.1)
% params.weight: data fidelity spatially variable weight (recommended = magnitude_data). 
% params.regweight: regularization spatially variable weight.
% params.precond: preconditionate solution by smart initialization
%
% Output: out - structure with the following fields:
% out.x: calculated susceptibility map
% out.iter: number of iterations needed
% out.time: total elapsed time (including pre-calculations)
%
% Based on the code by Bilgic Berkin at http://martinos.org/~berkin/software.html
% Modified by Carlos Milovic in 2017.06.05
% Last modified by Carlos Milovic in 2020.07.07
% 
% Adapted from code from FANSI toolbox: https://gitlab.com/cmilovic/FANSI-toolbox
% 
% modified by Xu Li for the JHUKKI_QSM_Toolbox, 2021.05.26
% added merit option

tic

    % Required parameters
    lambda = params.alpha1;
    kernel = params.D;
    phase = params.input;

    % Optional parameters
    if isfield(params,'mu1')
         mu = params.mu1;
    else
        mu = 100*lambda;
    end
    
    if isfield(params,'mu2')
         mu2 = params.mu2;
    else
        mu2 = 1.0;
    end
    N = size(params.input);

    if isfield(params,'maxOuterIter')
        num_iter = params.maxOuterIter;
    else
        num_iter = 150;
    end
    
    if isfield(params,'tol_update')
       tol_update  = params.tol_update;
    else
       tol_update = 0.1;
    end

    
    if isfield(params,'regweight')
        regweight = params.regweight;
        if length(size(regweight)) == 3
            regweight = repmat(regweight,[1,1,1,3]);
        end
    else
        regweight = ones([N 3]);
    end
    
    if ~isfield(params,'delta_tol')
        delta_tol = 1e-6;
    else
        delta_tol = params.delta_tol;
    end
    
    if isfield(params,'merit')
        merit = params.merit;
    else
        merit = 1;
    end
    
    W = params.weight.*params.weight; % for computational efficiency
    
    % Variable initialization
    z_dx = zeros(N, 'single');
    z_dy = zeros(N, 'single');
    z_dz = zeros(N, 'single');

    s_dx = zeros(N, 'single');
    s_dy = zeros(N, 'single');
    s_dz = zeros(N, 'single');

    x = zeros(N, 'single');

    if isfield(params,'precond')
        precond = params.precond;
    else
        precond = true;
    end
    
    if precond
        z2 =  W.*phase./(W+mu2); % start with something similar to the input phase, weighted to reduce noise
    else
        z2 = zeros(N,'single');
    end
    s2 = zeros(N,'single'); 



% Define the operators
[k1, k2, k3] = ndgrid(0:N(1)-1,0:N(2)-1,0:N(3)-1);

E1 = 1 - exp(2i .* pi .* k1 / N(1));
E2 = 1 - exp(2i .* pi .* k2 / N(2));
E3 = 1 - exp(2i .* pi .* k3 / N(3));

E1t = conj(E1);
E2t = conj(E2);
E3t = conj(E3);

EE2 = E1t .* E1 + E2t .* E2 + E3t .* E3;
K2 = abs(kernel).^2;

%tic
    ll = lambda/mu;
for t = 1:num_iter
    % update x : susceptibility estimate
    tx = E1t .* fftn(z_dx - s_dx);
    ty = E2t .* fftn(z_dy - s_dy);
    tz = E3t .* fftn(z_dz - s_dz);
    
    x_prev = x;
    Dt_kspace = conj(kernel) .* fftn(z2-s2);
    x = real(ifftn( (mu * (tx + ty + tz) + mu2*Dt_kspace) ./ (eps + mu2*K2 + mu * EE2) ));

    x_update = 100 * norm(x(:)-x_prev(:)) / norm(x(:));
    disp(['Iter: ', num2str(t), '   Update: ', num2str(x_update)])
    
    if x_update < tol_update || isnan(x_update)
        break
    end
    
    if t < num_iter
        % update z : gradient varible
        Fx = fftn(x);
        x_dx = real(ifftn(E1 .* Fx));
        x_dy = real(ifftn(E2 .* Fx));
        x_dz = real(ifftn(E3 .* Fx));
        
        z_dx = max(abs(x_dx + s_dx) - regweight(:,:,:,1)*ll, 0) .* sign(x_dx + s_dx);
        z_dy = max(abs(x_dy + s_dy) - regweight(:,:,:,2)*ll, 0) .* sign(x_dy + s_dy);
        z_dz = max(abs(x_dz + s_dz) - regweight(:,:,:,3)*ll, 0) .* sign(x_dz + s_dz);
    
        % update s : Lagrange multiplier
        s_dx = s_dx + x_dx - z_dx;
        s_dy = s_dy + x_dy - z_dy;            
        s_dz = s_dz + x_dz - z_dz;  
        
        
        rhs_z2 = mu2*real(ifftn(kernel.*Fx)+s2  );
        z2 =  rhs_z2 ./ mu2 ;

        % Newton-Raphson method
        delta = inf;
        inn = 0;
        while (delta > delta_tol && inn < 50)
            inn = inn + 1;
            norm_old = norm(z2(:));
            
            update = ( W .* sin(z2 - phase) + mu2*z2 - rhs_z2 ) ./ ( W .* cos(z2 - phase) + mu2 );            
        
            z2 = z2 - update;     
            delta = norm(update(:)) / norm_old;
        end        
        disp(delta)
        
        s2 = s2 + real(ifftn(kernel.*Fx)) - z2;
    end
    
    if merit
        wres = params.weight.*(exp(1i*real(ifftn(kernel.*Fx))) - exp(1i*phase));
        wres = wres - mean(wres(params.mask>0));  % normalized residual
        a = wres(params.mask>0);
        factor = std(abs(a))*6;                 % std of abs(res) * 6 sigma
        wres = abs(wres)/factor;                % divided 6 sigma
        wres(wres<1) = 1;                       % within 6 sigma ==> 1, does not change
        wres_t = 2;
        W = (params.weight./wres.^wres_t).^2;
    end
end
% Extract output values
out.time = toc;toc
out.x = x;
out.iter = t;


end
