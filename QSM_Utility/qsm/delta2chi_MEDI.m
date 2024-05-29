function [x, cost_reg_history, cost_data_history] =  delta2chi_MEDI(deltaB, Params, m, Mask, lambda, merit, edgePer, smv_rad)  
% [x, cost_reg_history, cost_data_history] =  delta2chi_MEDI(deltaB, Params, m, maskErode, lambda, merit, edgePer, smv_rad)   
%   QSM using MEDI_L1
%   deltaB: field change in ppm
%   Params: parameters structures
%   m: data term weighting
%   maskErode: brain Mask
%   lambda: for regularization
%   merit: fine tuning
%   edgePer: percentage of voxel as structure edges
%   smv_rad: smv radius for smv filtered option, default 0
% 
% modified according to MEDI_L1.m by Tian Liu, Ref: T. Liu et al. MRM 2013;69(2):467-76
% http://weill.cornell.edu/mri/pages/qsm.html
% 
% modified by Xu Li, 2021-06-11
% Added in SMV filtered option as in MSDI
% Acosta-Cabronero J et. al, Neuroimage 2018
% 2021-09-15, X.L. bug fix
% 2024, X.L. parameter fix

%%%%%%%%%%%%%%% weights definition %%%%%%%%%%%%%%
if nargin < 5
    lambda = 1000;
    merit = 0;
    edgePer = 0.9;
    smv_rad = 0;
elseif nargin < 6
    merit = 0;
    edgePer = 0.9; 
    smv_rad = 0;
elseif nargin < 7
    edgePer = 0.9;
    smv_rad = 0;
elseif nargin < 8
    smv_rad = 0;
end

% internal parameters for L1 solver
cg_max_iter = 100;              % for cg loop, or 100
cg_tol = 0.01;                  % cg tolerance, or 0.01

max_iter = 10;                  % outer loop
tol_norm_ratio = 0.1;

gradient_weighting_mode = 1;

grad = @cgrad;
div = @cdiv;

% kernel
D = conv_kernel_rot_c0(Params, Params.TAng);     
D = ifftshift(D);

% smv filtered option
if smv_rad > 0
    N = size(Mask);
    SMV = gen_SMVkernel_voxel_scaled(Params, N, smv_rad);
    D = D.*SMV;
    deltaB = real(ifftn(SMV.*fftn(deltaB))).*Mask;
end

% b0 and weighting
b0 = m.*exp(1i*deltaB);                         % change RDF to deltaB in ppm
wG = gradient_mask(gradient_weighting_mode, m, Mask, grad, Params.voxSize, edgePer);

iter=0;
x = zeros(Params.sizeVol); 
res_norm_ratio = Inf;
cost_data_history = zeros(1,max_iter);
cost_reg_history = zeros(1,max_iter);

e=0.000001; %a very small number to avoid /0

% main loop
while (res_norm_ratio>tol_norm_ratio)&&(iter<max_iter)
tic
    iter=iter+1;
    Vr = 1./sqrt(abs(wG.*grad(real(x),Params.voxSize)).^2+e);       % 1/abs(M*G*x_n)
    
    % complex
    w = m.*exp(1i*real(ifftn(D.*fftn(x))));                           % w: W*exp(i*D*x_n)
    
    % regularization term
    reg = @(dx) div(wG.*(Vr.*(wG.*grad(real(dx),Params.voxSize))),Params.voxSize);  
    
    % data fidelity term
    fidelity = @(dx)2*lambda*real(ifftn(D.*fftn(conj(w).*w.*real(ifftn(D.*fftn(dx))))));

    A =  @(dx) reg(dx) + fidelity(dx);       
    b = reg(x) + 2*lambda*real(ifftn(D.*fftn( conj(w).*conj(1i).*(w-b0))));

    dx = real(cgsolve(A, -b, cg_tol, cg_max_iter, 10));      % CG solver, solve A*dx = b
    res_norm_ratio = norm(dx(:))/norm(x(:));
    x = x + dx;

    wres=m.*exp(1i*(real(ifftn(D.*fftn(x))))) - b0;         % weighted residual, should be small

    cost_data_history(iter) = norm(wres(:),2);
    cost=abs(wG.*grad(x));
    cost_reg_history(iter) = sum(cost(:));

    
    if merit
        wres = wres - mean(wres(Mask(:)==1));   % unbias
        a = wres(Mask(:)==1);
        factor = std(abs(a))*6;     % 6 sigma
        
        wres = abs(wres)/factor;    % normalization        
        wres(wres<1) = 1;
        
        m = m./(wres.^2);        
        b0 = m.*exp(1i*deltaB);
        
    end
    
    fprintf('iter: %d; res_norm_ratio:%8.4f; cost_L2:%8.4f; cost_Reg:%8.4f.\n',iter, res_norm_ratio,cost_data_history(iter), cost_reg_history(iter));
toc
    
end

x = x.*Mask;

function SMV = gen_SMVkernel_voxel_scaled(Params, N, smv_rad)
  
    % Create k-space kernel with radius smv_rad 
    if N ~= Params.sizeVol
        disp('check Params.sizeVol')
    end
    % under N coordinate
    [Y,X,Z] = meshgrid(-floor(N(2)/2):ceil(N(2)/2-1),-floor(N(1)/2):ceil(N(1)/2-1),-floor(N(3)/2):ceil(N(3)/2-1));

    X = X * Params.voxSize(1);
    Y = Y * Params.voxSize(2);
    Z = Z * Params.voxSize(3);
    
    smv = (X.^2 + Y.^2 + Z.^2) <= smv_rad^2;
    smv = smv / sum(smv(:));                    % normalized 

    smv_kernel = zeros(size(X));
    smv_kernel(1+floor(end/2),1+floor(end/2),1+floor(end/2)) = 1;    
    smv_kernel = smv_kernel - smv;              % delta - smv

    SMV = fftn(fftshift(smv_kernel));           % kernel in k-space 

end


end

