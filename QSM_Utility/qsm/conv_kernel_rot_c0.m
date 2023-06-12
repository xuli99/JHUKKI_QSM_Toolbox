function [C, KSq] = conv_kernel_rot_c0(Params, R, datatype, kerneltype)
% [C, KSq] = conv_kernel_rot_c0(Params, R, datatype, kerneltype)
% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
%
% Updated X.L., 2019-07
% Updated X.L., 2020-10-10

if nargin < 3
    datatype = 'double';
    kerneltype  = 0;
elseif nargin < 4
    kerneltype = 0;     % default kernel D
end

if isempty(datatype)
    datatype = 'double';
end

if ~isfield(Params, 'permuteflag')
    Params.permuteflag = 1;     % version before v2.8
end
    
warning off all

Nx = Params.sizeVol(1);
Ny = Params.sizeVol(2);
Nz = Params.sizeVol(3);

dkx = 1/Params.fov(1);
dky = 1/Params.fov(2);
dkz = 1/Params.fov(3);

%% convolution kernel 

kx = linspace(-Nx/2, Nx/2-1, Nx).*dkx;
ky = linspace(-Ny/2, Ny/2-1, Ny).*dky;
kz = linspace(-Nz/2, Nz/2-1, Nz).*dkz;

kx = cast(kx, datatype);
ky = cast(ky, datatype);
kz = cast(kz, datatype);

[KX_Grid, KY_Grid, KZ_Grid] = meshgrid(ky, kx, kz);
KSq = KX_Grid.^2 + KY_Grid.^2 + KZ_Grid.^2;

if isfield(Params, 'Tsominv') && isfield(Params, 'Tpom')
    
    H0 = [0, 0, 1]';                   % Lab space XYZ  
    extra = [0,-1,0;-1,0,0;0,0,1];     % 
    Hsub = -extra*Params.Tsominv*R'*Params.Tpom*H0; 
        
    R31 = Hsub(1);
    R32 = Hsub(2);
    R33 = Hsub(3);

elseif isfield(Params, 'sliceOri')
    
    switch Params.sliceOri
        case 1  % axial
            H0 = [0, 0, 1]';     % Lab space XYZ
        case 2
            H0 = [1, 0, 0]';    % need further testing
        case 3  % coronal 
            H0 = [0, 1, 0]';    % in Y directoin
        otherwise
            error('unknown slice orientation.')
    end

    if isfield(Params, 'datatype')
        if strcmpi(Params.datatype, '2dseq')
            H0 = [0, 0, 1]';     % for 2dseq, apply TAng on LPS coordinate
        end
    end

    Hsub = R'*H0;              
    
    R31 = Hsub(1);
    R32 = Hsub(2);
    R33 = Hsub(3);
    
else
    % default use axial H0 = [0, 0, 1], Hsub = R'*H0  
    R31 = R(3, 1);
    R32 = R(3, 2);
    R33 = R(3, 3);
end

KZp_Grid = R31.*KX_Grid + R32.*KY_Grid + R33.*KZ_Grid;  % i.e. ((R'*H0).*K), in subject frame

if kerneltype == 0              
    C = 1/3 - (KZp_Grid.^2)./(KSq);         % normal D2 kernel
    C(isnan(C)) = 0;  
elseif kerneltype == 1    
    C = 1/3.*KSq - KZp_Grid.^2;             % Laplacian DL kernel
end