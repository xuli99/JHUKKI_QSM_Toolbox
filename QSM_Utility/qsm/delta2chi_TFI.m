function [x, TFI_params] = delta2chi_TFI(deltaB, D, DPWeight, mask, TFI_params)
% [x, TFI_params] = delta2chi_TFI(deltaB, D, DPWeight, mask, TFI_params)
%
% Auto-TFI method using Gaussian Newton's method (NGCG)
% with auto determined P_map, plus Tikhonov regularization as in LN-QSM
% Ref: Liu et. al, 2020, 83, 271 (Auto-TFI)
%      Sun et al, 2018, 179, 166  (LN-QSM)
%      should add merit, also consider using mcTFI
%
% input of deltaB, DPWeight, D, TFI_params
% DPWeight is the apriori maps used to extract tissue boundary information
%
%     TFI_params.R2star_thresh = 30;             % Hz, threshold for detecting soft tissue, for old TFI
%     TFI_params.Ps            = 30;             % Preconditioning value
%     TFI_params.P_map                           % default is 1 in M, Ps otherwise, Auto estimated according to chi_est
%     TFI_params.lambda1  = 1e-3;  % regularization
%     TFI_params.lambda2  = 0.1;   % if AutoRef
%     TFI_params.lambda3  = 1e-5;  % Tikhonov regularization as in LN-QSM, on maskBET
%     TFI_params.epsilon  = 1e-6*TFI_params.P_map.^2;
%     TFI_params.tol      = 0.01;  % tol for GNCG
%     TFI_params.itermax  = 300;
%     TFI_params.cgtol    = 0.01;  % cg tol
%     TFI_params.cgitermax= 100;   % main iteration
%     TFI_params.verbose  = 1;
%
% % last update 2017-07-15 X.L.
% % Using Newton's method, Updated 2018-02-16, X.L.
% % updated 2020-05-16, updated for auto-determined continuous P_map X.L.
% % updated 2020-05-29, updated to add Tikhonov reg

TV = TVOP;              % total varition class (3D finite difference operator)
N = size(deltaB);       % [Nx, Ny, Nz]

% getting wG edge weighting
wG = gradient_mask_all(DPWeight, mask, TFI_params.apEdge);    % using Magnitude to get a priori, use image based
wG = wG.*mask;                                                % with implicit expansion

% TFI_params
P           = TFI_params.P_map;   % preconditioner
lambda1     = TFI_params.lambda1;
lambda2     = TFI_params.lambda2;   % for AutoRef
lambda3     = TFI_params.lambda3;   % for Tikhonov                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            = TFI_params.lambda3;   % for Tikhonov Reg
W           = DPWeight;
max_iter    = TFI_params.itermax;         % outer loop
tol         = TFI_params.tol;  % outer loop
cg_max_iter = TFI_params.cgitermax; 
cg_tol      = TFI_params.cgtol;
epsi        = TFI_params.epsilon;         % very small number to avoid /0, now adaptive
verbose     = TFI_params.verbose;

% initialization
b0 = deltaB;            % f
iter = 0;
x = zeros(N);           % yn, final x = P*yn
res_norm_ratio = Inf;
cost_data_history = zeros(1,max_iter);
cost_reg_history = zeros(1,max_iter);

% main loop
while (res_norm_ratio > tol) && (iter < max_iter)
tic
    iter=iter+1;
    Vr = 1./sqrt((wG.*(TV*(P.*real(x)))).^2 + epsi);    % 1./sqrt(abs(wG*G*Pyn).^2 + epsi), NwG, Nx3
        
    % regularization term
    reg = @(dx) lambda1*P.*(TV'*(wG.*(Vr.*(wG.*(TV*(P.*real(dx))))))) + lambda3*P.*mask.*P.*(dx);
    
    % data fidelity term
    fidelity = @(dx) P.*real(ifftn(D.*fftn(W.^2.*real(ifftn(D.*fftn(P.*dx)))))); 

    A =  @(dx) reg(dx) + fidelity(dx);       
    b = -reg(x) + P.*real(ifftn(D.*fftn(W.^2.*(b0 - real(ifftn(D.*fftn(P.*x)))))));

    dx = (cgsolve(A, b, cg_tol, cg_max_iter, verbose));      % CG solver
    res_norm_ratio = norm(dx(:))/norm(x(:));
    x = x + dx;

    wres=W.*(real(ifftn(D.*fftn(P.*x))) - b0);

    cost_data_history(iter) = norm(wres(:));
    cost = lambda1*abs(wG.*(TV*(P.*x)));
    cost_reg_history(iter) = sum(cost(:));
       
    fprintf('iter: %d; res_norm_ratio:%8.4f; cost_DF:%8.4f; cost_Reg:%8.4f.\n',iter, res_norm_ratio,cost_data_history(iter), cost_reg_history(iter));
toc
    
end

x = (P.*x);   % final solution
TFI_params.iter = iter;
TFI_params.cost_data_history = cost_data_history;
TFI_params.cost_reg_history = cost_reg_history;
TFI_params.res = res_norm_ratio;


end  