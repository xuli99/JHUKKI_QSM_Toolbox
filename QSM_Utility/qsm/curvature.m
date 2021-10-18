function Kappa = curvature(muTest, consistency, regularization, plotflag)
% cubic spline differentiation to find Kappa (largest curvature) 
if nargin < 3
    plotflag = 1;
end

if plotflag == 1
    figure(20), plot(consistency, regularization, 'marker', '*')
end

eta = log(regularization.^2);
rho = log(consistency.^2);

M = [0 3 0 0;0 0 2 0;0 0 0 1;0 0 0 0];

pp = spline(muTest, eta);
ppd = pp;

ppd.coefs = ppd.coefs*M;    % first derivative
eta_del = ppval(ppd, muTest);
ppd.coefs = ppd.coefs*M;    % second derivative
eta_del2 = ppval(ppd, muTest);

pp = spline(muTest, rho);
ppd = pp;

ppd.coefs = ppd.coefs*M;
rho_del = ppval(ppd, muTest);
ppd.coefs = ppd.coefs*M;
rho_del2 = ppval(ppd, muTest);

Kappa = 2 * (rho_del2 .* eta_del - eta_del2 .* rho_del) ./ (rho_del.^2 + eta_del.^2).^1.5;

if plotflag == 1
    figure(21), semilogx(muTest, Kappa, 'marker', '*');
end