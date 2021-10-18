function maskErode = imerode3(mask, radius, flag3d)
% Faster with build-in!
% Updated 2014-12-15
% added 2D compatibility, 2017-04-11
% updated, X.L., 2017-08-04
% updated, added 3dflag

if nargin < 3
    flag3d = 0;
end

if (radius < 6) && (~flag3d)
    [xx,yy,zz] = ndgrid(-radius:radius);
    if size(mask, 3) > 1
        nhood = sqrt(xx.^2 + yy.^2 + zz.^2) <= radius;
    else
        nhood = sqrt(xx.^2 + yy.^2) <= radius;
    end
    maskErode = imerode(mask,nhood);
else
    % do 3D eroding in k space
    [Ny, Nx, Nz] = size(mask);

    sesize = radius*2 + 1;
    x = linspace(-radius, radius, sesize);
    y = linspace(-radius, radius, sesize);
    z = linspace(-radius, radius, sesize);

    [X_Grid, Y_Grid, Z_Grid] = meshgrid(x, y, z);  % mesh in space
    R_Grid = sqrt(X_Grid.^2 + Y_Grid.^2 + Z_Grid.^2);

    se = R_Grid <= radius;

    fftsize = [Ny+sesize-1, Nx+sesize-1, Nz+sesize-1];
    temp = fftn(mask, fftsize).*fftn(se, fftsize);
    temp = real(ifftn((temp)));

    maskErode = temp((radius+1):(end-radius), ...
                 (radius+1):(end-radius), (radius+1):(end-radius));
    maskErode = round(maskErode);

    maskErode = (maskErode == sum(se(:)));

    maskErode = maskErode > 0;
end