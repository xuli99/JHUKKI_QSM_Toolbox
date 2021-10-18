function [chi_res, flag, relres, iter] = delta2chi_NDI(deltaB, kernel, DPWeight, maskErode, NDIParams)
% [chi_res, flag, relres, iter] = delta2chi_NDI(deltaB, kernel, DPWeight, maskErode, NDIParams)
% 
%% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
% 
% Quantitative Susceptibility Mapping (QSM) using
% with NDI method, Ref: Polak et al, NMRB, 2020; e4271
%                       Kames et al, ISMRM, 2021, 3981, HANDI, not working yet
%
% Input:  
%       NDIParams includes fields:
%             NDIParams.step_size = 1;       % gradient descent step size
%             NDIParams.num_iter = 50;      % number of iteration
%             NDIParams.tik  = 1;            % adding/no tikhonov regularization
%             NDIParams.lambda = 1e-3;       % 
%             NDIParams.tol = 1;           % tolerance
%               NDIParams.gamma = Params.gamma;          % constants
%               NDIParams.B0 = Params.B0
%               NDIParams.TEs = Params.TEs;
%
% Output:
%     chi_res_res: isotropic susceptibility map
%     flag: convergence flag: 0 means converged (similar as lsqr)
%     relres: final relative residual in terms of gradient
%     iter: total iteration number
%
%  updated 2020-06-17
%  updated 2021-05-06
%  updated 2021-05-24, Tried Hessian Accelaration

N = size(deltaB);       % for only 1 orientation
scale_factor = NDIParams.gamma*NDIParams.B0*2*pi*mean(NDIParams.TEs)*1e-6;

num_iter = NDIParams.num_iter;
step_size = NDIParams.step_size;
tik = NDIParams.tik;             % adding tikhonov regularization
lambda = NDIParams.lambda;      % 
tol = NDIParams.tol;

flag = 1;                           % not converging in num_iter
if isfield(NDIParams, 'verbose')
    verbose = NDIParams.verbose;
else
    verbose = 1;
end

% % use precalculated kernel
% kernel = conv_kernel_rot_c0(Params, Params.TAng, NDIParams.datatype);
% kernel = ifftshift(kernel);

M2 = DPWeight.^2;                   % weighting square
phs_use = deltaB*scale_factor;       % in radian

precond = 0;
if precond == 1
    chi_res = M2.*phs_use;
else
    chi_res = zeros(N);
end
grad_prev = 0;

disp('solving QSM inverse problem using ndi ...')

tic
for iter = 1:num_iter
    Dx = ifft(ifft(ifft(kernel.*fftn(chi_res), [], 1), [], 2), [], 3);
    temp = M2 .* sin(Dx - phs_use);
    grad_f = 2 * sum(ifft(ifft(ifft(kernel.*fft(fft(fft(temp, [], 1), [], 2), [], 3), [], 1), [], 2), [], 3), 4);
    
%     temp1 = M2 .* cos(Dx - phs_use);
%     denorm = abs(2 * sum(ifft(ifft(ifft(kernel.^2.*fft(fft(fft(temp1, [], 1), [], 2), [], 3), [], 1), [], 2), [], 3), 4));
%     mu = sqrt(max(abs(grad_f(:))));
        
    if tik == 1
        grad_f = (grad_f + 2*lambda.*chi_res);
        % grad_f = (grad_f + 2*lambda.*chi_res)./(denorm + 2*lambda + mu);
    else
        % grad_f = grad_f./(denorm + mu);
    end
    
    chi_res = chi_res - step_size * real(grad_f);

    update_grad = rmse(grad_prev, grad_f);

    if verbose
        disp(['iter: ', num2str(iter), '   grad update:', num2str(update_grad)])
    end
    
    if update_grad < tol
        flag = 0;
        break
    end

    grad_prev = grad_f;
end
toc

disp(['final iter: ', num2str(iter), '   grad update:', num2str(update_grad)])

res = ifft(ifft(ifft(kernel.*fftn(chi_res), [], 1), [], 2), [], 3) - phs_use;
relres = norm(res(:).*maskErode(:))./norm(phs_use(:).*maskErode(:));

chi_res = chi_res.*maskErode/scale_factor;

