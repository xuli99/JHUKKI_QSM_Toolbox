function res = D(image)

%
% res = D(image)
%
% image = a 3D image
%
% This function computes the finite difference transform of the image
%
% Related functions:
%       adjD , invD 
%
%

[sx,sy,sz] = size(image);

Dx = image([2:end,end],:,:) - image;
Dy = image(:,[2:end,end],:) - image;
Dz = image(:,:,[2:end,end]) - image;

%res = [sum(image(:))/sqrt(sx*sy); Dx(:);  Dy(:)]; 
res = cat(4,Dx,Dy,Dz);


