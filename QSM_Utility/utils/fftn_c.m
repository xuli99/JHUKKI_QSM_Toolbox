function [out] = fftn_c(in)
% function [out] = fftn_c(in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2d & 3d centered Fourier transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tatiana Latychevskaia and Hans-Werner Fink
% "Practical algorithms for simulation and reconstruction of digital in-line holograms",
% Appl. Optics 54, 2424 - 2434 (2015)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modifyed by Xu Li, 2020
%
% Author:           Xu Li, PhD
% Email:            xuli@mri.jhu.edu
%================================

ndim = ndims(in);
N = size(in);
in = cast(in, 'double');

if ndim == 3
    [x,y,z] = meshgrid(0:N(2)-1, 0:N(1)-1, 0:N(3)-1);   
    f1 = exp(1i*pi*(x+y+z));  % same size as in
elseif ndim == 2
    [x,y] = meshgrid(0:N(2)-1, 0:N(1)-1);   
    f1 = exp(1i*pi*(x+y));  % same size as in
end

out = f1.*fftn(f1.*in)*exp(-1i*sum(N)*pi/2);
 