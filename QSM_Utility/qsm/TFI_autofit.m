function P_map = TFI_autofit(D_map, chi_est, R2starMap, maskErode, mask_LBV)
% P_map = TFI_autofit(D_map, chi_est, R2starMap, maskErode, mask_LBV)
%% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
%
% Purpose:
% fit the TFI P_map functino parameters, i.e. sigma0 ,r0, sigma1, sigma2, s1, s2
% With two models,  1. cubic decay model for ~maskErode
%                   2. Logistic model for maskErode
% Ref. Liu et. al, MRM, 2020, 83:271
% 
% Input: D_map, chi_est, R2starMap, maskErode
% Output: P_map

debug_flag = false;
options = optimoptions(@fmincon,'Display', 'notify');

%% fitting the cubic decay model
% learning (X, Y) of (D_unique_bin, chi_med0_bin)
[D_unique_bin, chi_med0_bin] = extract_training_data(D_map, chi_est, ~maskErode, 1);
nanindex = isnan(chi_med0_bin);  % x may have gaps
D_unique_bin(nanindex) = [];
chi_med0_bin(nanindex) = [];

lb = 10*eps*ones(2,1);
ub = Inf*ones(2,1);
x0 = rand(2,1);
theta1 = fmincon(@(t)CubicCostFunction(t, double(D_unique_bin), chi_med0_bin), x0, [], [], [], [], lb, ub, [], options);

if debug_flag
     chi_med0_pred_bin = CubicPred(theta1, D_unique_bin);
     fig100 = figure(100); 
     plot(D_unique_bin, chi_med0_bin, 'b*', D_unique_bin, chi_med0_pred_bin, 'r-')  % check measurement data  
     legend('observed', 'fitted')
     pause()
     close(fig100)
end

%% fitting the sigmoid function
[R2star_unique_bin, chi_med1_bin] = extract_training_data(R2starMap, chi_est, mask_LBV, 1);
nanindex = isnan(chi_med1_bin);
R2star_unique_bin(nanindex) = [];
chi_med1_bin(nanindex) = [];

R2star_thresh = 150; chi_med1_thresh = 1;
fit_index = (R2star_unique_bin <= R2star_thresh) & (chi_med1_bin < chi_med1_thresh);

A = zeros(4,4); A(1,1:2) = [1, -1];
b = zeros(4,1);
lb = [max(chi_med0_bin(:))/40, eps, 10, eps];     % make P_map ~ 40
ub = [Inf, Inf, 150, 50];           % 10 < s1 < 150; s2 < 50
x0 = rand(4,1);
theta2 = fmincon(@(t)SigmoidCostFunction(t, double(R2star_unique_bin(fit_index)), chi_med1_bin(fit_index)), x0, A, b, [], [], lb, ub, [], options);

if debug_flag
     chi_med1_pred_bin = SigmoidPred(theta2, R2star_unique_bin);
     fig101 = figure(101); 
     plot(R2star_unique_bin, chi_med1_bin, 'b*', R2star_unique_bin, chi_med1_pred_bin, 'r-')  % check measurement data  
     legend('observed', 'fitted')
     pause()
     close(fig101)
end

%% calculating the final P_map with theta1 and theta2
sigma1 = theta2(1); 
P_map = CubicPred(theta1, D_map).*(~maskErode) + SigmoidPred(theta2, R2starMap).*maskErode;
P_map = P_map./sigma1;      % normalized around 1

%% internal functions
function [x, y] = extract_training_data(x_map, y_map, mask, BinWidth)
    [N, edges, bin] = histcounts(x_map(mask), 'BinWidth', BinWidth);  % bin has the same size of x_map(mask)
    x = edges(1:end-1)+BinWidth/2;
    x = x(:);
    y0 = abs(y_map(mask));              % take abs                                
    y = zeros(length(N), 1);                              
    for i_bin = 1:length(N)
        y(i_bin) = median(y0(bin == i_bin));
    end        
end

function y = CubicPred(theta, X)
    y = theta(1).*(1+X./theta(2)).^-3;
end

function J = CubicCostFunction(theta, X_t, y_t)
    % theta has sigma0 and r0, with both sigma0 >= eps, r0 >= eps
    y_t = y_t(:);
    y_pred = CubicPred(theta, X_t(:));
    J = 0.5*(y_pred - y_t)'*(y_pred - y_t);
end

function y = SigmoidPred(theta, X_t)
    % theta has sigma1, sigma2, and s1, s2
    y = (theta(2) - theta(1))./(1 + exp(-(X_t - theta(3))./theta(4))) + theta(1);
end

function J = SigmoidCostFunction(theta, X_t, y_t)
    y_t = y_t(:);
    y_pred = SigmoidPred(theta, X_t(:));
    J = 0.5*(y_pred - y_t)'*(y_pred - y_t);
end

end