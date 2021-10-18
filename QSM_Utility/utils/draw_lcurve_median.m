function [ Kappa ] = draw_lcurve_median( Lambda, regularization, consistency, fig )
% Draw the L-curve in log domain.
% Curvature is calculated in the log domain. Then, it is filtered by a moving median.
% Both curvature estimations are plotted (raw = blue, filtered = red) as function of Lambda.
% This function outputs the filtered curvature Kappa.
%
% Parameters:
% Lambda: vector with regularization weight values
% regularization: vector with regularization costs (use the appropiate function to evaluate such costs)
% consistency: vector with data fidelity costs (use the appropiate function to evaluate such costs)
% fig: figure number to display the L-curve and the curvature.
%
% Last modified by Carlos Milovic, 08.07.2020 


% First plot the L-curve in the log domain
figure(fig), subplot(1,2,1), plot(log(consistency), log(regularization), 'marker', '*')
set(gca,'FontSize',24)
set(gcf,'Color','white')
xlabel('Log Fidelity cost')
ylabel('Log Regularization cost')

% Calculate the curvature in the log domain, using splines
[ Kappa ] = calc_curv_spline( Lambda, regularization, consistency, false )

figure(fig), subplot(1,2,2), semilogx(Lambda, Kappa,'b', 'marker', '*')

hold on;
Kappa = medfilt2( Kappa, [1 5] );

semilogx(Lambda, Kappa,'r', 'marker', '*')
hold off;
set(gca,'FontSize',24)
set(gcf,'Color','white')
xlabel('Regularization weight')
ylabel('Curvature')
legend('Original','Median filtered')
legend('Location','northeast')
end

