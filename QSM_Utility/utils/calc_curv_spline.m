function [ Kappa ] = calc_curv_spline( Lambda, regularization, consistency, isLinear, plotflag)
% Calculate the curvature of the L-curve by using spline differentiation
% This function outputs the curvature Kappa.
%
% Parameters:
% Lambda: vector with regularization weight values
% regularization: vector with regularization costs (use the appropiate function to evaluate such costs)
% consistency: vector with data fidelity costs (use the appropiate function to evaluate such costs)
% isLinear: true for a calculation in the linear domain and false to use the log domain.
%
% 
% Based on the code by Bilgic Berkin at http://martinos.org/~berkin/software.html
% Last modified by Carlos Milovic, 08.07.2020 
%
% cubic spline differentiation to find Kappa (largest curvature) 
% modifed by Xu Li, 2025-05-25

if nargin < 5
    plotflag = 0;
end

if isLinear == true
eta = regularization; 
rho = consistency;
else
eta = log(regularization.^2); 
rho = log(consistency.^2);
end

M = [0 3 0 0;0 0 2 0;0 0 0 1;0 0 0 0];

pp = spline(Lambda, eta);
ppd = pp;

ppd.coefs = ppd.coefs*M;
eta_del = ppval(ppd, Lambda);
ppd.coefs = ppd.coefs*M;
eta_del2 = ppval(ppd, Lambda);

pp = spline(Lambda, rho);
ppd = pp;

ppd.coefs = ppd.coefs*M;
rho_del = ppval(ppd, Lambda);
ppd.coefs = ppd.coefs*M;
rho_del2 = ppval(ppd, Lambda);


% Kappa = 2 * (rho_del2 .* eta_del - eta_del2 .* rho_del) ./ (rho_del.^2 + eta_del.^2).^1.5;
Kappa = abs(rho_del2 .* eta_del - eta_del2 .* rho_del) ./ (rho_del.^2 + eta_del.^2).^1.5;

index_opt = find(Kappa == max(Kappa));

disp(['Optimal lambda, consistency, regularization: ', num2str([Lambda(index_opt), consistency(index_opt), regularization(index_opt)])])

if plotflag 
    figure, subplot(1,2,1), plot(rho, eta, 'b-*')
    subplot(1,2,2), semilogx(Lambda, Kappa, 'r-*')
end

end

