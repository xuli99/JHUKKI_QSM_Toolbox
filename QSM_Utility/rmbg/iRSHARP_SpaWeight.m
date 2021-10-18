function [C,Grad]=iRSHARP_SpaWeight(data,Params, Wrapcount)
% [C,Grad]=iRSHARP_SpaWeight(data,Params, Wrapcount)
% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
% 
% Calculating the iRSHARP spatial weight

% Gaussian kernel
 Kwin = 5;      % window size, cubic
 khsize = floor((Kwin-1)/2);

 sigma = khsize;
 G = Gaussiankernel3d(Kwin, sigma); 

% Gradient
Grad = cgrad(data, Params.voxSize );    % phase gradient
Grad = sqrt(sum(Grad.^2, 4));           % sum of square
Grad = Grad./max(Grad(:));              % normalized gradient norm
data = data./max(data(:));              % normalized field shift
Wrapcount = Wrapcount./max(Wrapcount(:));    % normalized wrapcount

% combined weighting
MultiPhaseGrad =(1 + Grad).*(1 + abs(data)).*(1 + abs(Wrapcount));

% Smoothing with GaussianKernel
% convolution
N = size(data);

GKernel = zeros(N);
GKernel( floor(1+N(1)/2) - khsize : floor(1+N(1)/2) + khsize, ...
        floor(1+N(2)/2) - khsize : floor(1+N(2)/2) + khsize, ...
        floor(1+N(3)/2) - khsize : floor(1+N(3)/2) + khsize ) = G;

GKernel = fftn(fftshift(GKernel));
C = real(ifftn(fftn(MultiPhaseGrad).*GKernel));







