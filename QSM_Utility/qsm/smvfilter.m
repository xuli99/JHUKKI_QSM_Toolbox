function fimg = smvfilter(img, radius)
% spherical mean filter averaging on img, with
% 
%% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu

N = size(img);

% spherical kernel
kernelsize = 2*radius + 1;
[x, y, z] = meshgrid(-radius:radius, -radius:radius, -radius:radius);

kernel = sqrt((x.^2 + y.^2 + z.^2)) <= radius;
kernel = kernel./sum(kernel(:));

fftsize = (N + kernelsize - 1);

% averageing / convolution in real space, multiplication in k-space

temp = (fftn(img, fftsize)).*(fftn(kernel, fftsize));
temp = (ifftn((temp)));

fimg = temp((radius+1):(end-radius), (radius+1):(end-radius), (radius+1):(end-radius));