function G=Gaussiankernel3d(D, sigma)   
% function G=Gaussiankernel3d(D, sigma)     
% D can be window size or size of matrix
% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu

r = round((D-1)/2);                     % radius

if length(D) == 1
    [XX,YY,ZZ]=meshgrid(-r:r,-r:r,-r:r);  
elseif length(D) == 3
    [XX,YY,ZZ]=meshgrid(-r(1):r(1),-r(2):r(2),-r(3):r(3));  
else
    error('input window size of vector of matrix size.')
end

Rsq = XX.^2 + YY.^2 + ZZ.^2;
G=-Rsq/(sigma.^2*2);
G=exp(G)/(sigma*sqrt(2*pi));
G=G./sum(G(:));                         % normalized
