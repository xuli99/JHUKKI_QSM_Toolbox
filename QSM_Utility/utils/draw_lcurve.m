function [ Kappa ] = draw_lcurve( Lambda, regularization, consistency, fig )
% Draw the L-curve both in the linear and log domains.
% Calculate the curvature in the log domain and plot it as function of Lambda.
% This function outputs the curvature Kappa.
%
% Parameters:
% Lambda: vector with regularization weight values
% regularization: vector with regularization costs (use the appropiate function to evaluate such costs)
% consistency: vector with data fidelity costs (use the appropiate function to evaluate such costs)
% fig: figure number to display the L-curve and the curvature.
%
% Last modified by Carlos Milovic, 08.07.2020 
% modified by Xu Li, 2021

% First plot the L-curve in the linear domain
figure(fig), subplot(1,3,1), plot((consistency), (regularization), 'marker', '*')
set(gca,'FontSize',24)
set(gcf,'Color','white')
xlabel('Fidelity cost')
ylabel('Regularization cost')
% Now plot the L-curve in the log domain
figure(fig), subplot(1,3,2), plot(log(consistency), log(regularization), 'marker', '*')
set(gca,'FontSize',24)
set(gcf,'Color','white')
xlabel('Log Fidelity cost')
ylabel('Log Regularization cost')

% Calculate the curvature in the log domain, using splines
[ Kappa ] = calc_curv_spline( Lambda, regularization, consistency, false )

index_opt = find(Kappa == max(Kappa));
disp(['Optimum lambda, consistency, regularization: ', num2str([Lambda(index_opt), consistency(index_opt), regularization(index_opt)])])

% Plot the curvature as function of the regularization weight
figure(fig), subplot(1,3,3), semilogx(Lambda, Kappa, 'marker', '*')
set(gca,'FontSize',24)
set(gcf,'Color','white')
xlabel('Regularization weight')
ylabel('Curvature')

end

