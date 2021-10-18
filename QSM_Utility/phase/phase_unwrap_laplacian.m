function [unwrappedPhase, PhaseLaplacian]  = phase_unwrap_laplacian(dataPhase, Params, RefVox, zp, padsize)
% [unwrappedPhase, PhaseLaplacian]  = phase_unwrap_laplacian(dataPhase, Params, RefVox, zp, padsize)
%================================
% HEADER
%================================
%
% Name:             phase_unwrap_laplacian.m
% Author:           Xu Li, PhD
% Email:            xuli@mri.jhu.edu
%
%================================
% PURPOSE
% unwrap phase according to the laplacian method
% i.e. del2(theta) = cos(theta)*del2(sin(theta)) + sin(theta)*del2(cos(theta))
%  Updated by Xu Li, 2019-08-15, added tsvd threshold
%  Updated by Xu Li, 2020-10-01, updated with fftn_c
%  Updated by Xu Li, 2020-10-15, zp=2 matching STIsuite
% Updated 2021-06-27 X.L., cluster version

if ~isfield(Params, 'cluster')
    textWaitbar = 'Performing phase unwrapping';
    multiWaitbar(textWaitbar, 0, 'Color', 'b' );
end

warning off all

N = size(dataPhase); 
Ny = N(1);
Nx = N(2);
Nz = N(3);
VoxSize = Params.voxSize;  % VoxSize might need permute

sizeVol = Params.sizeVol;  % check wrt N
if sizeVol(1) ~= N(1)      
    VoxSize = VoxSize([2,1,3]);    % match dataPhase
end

if nargin < 3      
    zp = 2;         % code switch: 1 with discrete kernel, 2 with continuous k-space kernel
    RefVox = [floor(N(1)/2)+1, floor(N(2)/2)+1, floor(N(3)/2)+1]; 
    padsize = [12, 12, 12];
elseif nargin < 4
    zp = 2;    
    padsize = [12, 12, 12];
elseif nargin < 5
    padsize = [12, 12, 12];
end

if isempty(RefVox)
    RefVox = [floor(N(1)/2)+1, floor(N(2)/2)+1, floor(N(3)/2)+1]; 
end

RefFlag = 1;                % Do Ref in general
if RefVox == 0         
    RefFlag = 0;            % do not do Ref
end

%% padarray
% check if there are odd dimensions and do padding
% only slice number could be odd number
dimsOdd = mod(N, 2);
dataPhase = padarray(dataPhase, dimsOdd, 'replicate', 'post');  % padding to even dim first     
dataPhase = padarray(dataPhase, padsize, 0);

N = size(dataPhase);   % new dimentions after padding, now all even

%%
if zp == 1
    % using discrete kernel
    ksize = [3, 3, 3];               
    khsize = (ksize-1)/2;
    
    kernelcore = [];
    kernelcore(:,:,1) = [0 0 0; 0 1 0; 0 0 0];
    kernelcore(:,:,2) = [0 1 0; 1 -6 1; 0 1 0];
    kernelcore(:,:,3) = [0 0 0; 0 1 0; 0 0 0];

    Kernel = zeros(N);
    Kernel( 1+N(1)/2 - khsize(1) : 1+N(1)/2 + khsize(1), 1+N(2)/2 - khsize(2) : 1+N(2)/2 + khsize(2), ...
        1+N(3)/2 - khsize(3) : 1+N(3)/2 + khsize(3) ) = kernelcore;

    lap_opK = fftn_c(Kernel);          
    
else
    % using k-space continuous kernel
    FOV = VoxSize.*N;     
    %-------------------   calculate the k^2 in k space
    % Get k space cooridnate
    dky = 1/FOV(1);
    dkx = 1/FOV(2);
    dkz = 1/FOV(3);

    ky = linspace(-N(1)/2, N(1)/2-1, N(1)).*dky;
    kx = linspace(-N(2)/2, N(2)/2-1, N(2)).*dkx;
    kz = linspace(-N(3)/2, N(3)/2-1, N(3)).*dkz;
    
    [KX_Grid, KY_Grid, KZ_Grid] = meshgrid(kx, ky, kz);  % mesh in k space
    lap_opK = (KX_Grid.^2 + KY_Grid.^2 + KZ_Grid.^2);   % lap_opK
end

if ~isfield(Params, 'cluster')
    multiWaitbar(textWaitbar, 0.4, 'Color', 'b' );
end

% inverse laplacian kernel       
inv_lap_opK = zeros(N);
lap_opKthresh = 0; % 1e-3;       % 0 or threshold 1e-3                                                 % tsvd, X.L. 2019
inv_lap_opK(abs(lap_opK)>lap_opKthresh) = 1./lap_opK(abs(lap_opK)>lap_opKthresh);                   % inverse of kernel in k-space

if ~isfield(Params, 'cluster')
    multiWaitbar(textWaitbar, 0.6, 'Color', 'b' );
end

% laplacian based phase unwrapping 
PhaseSin = sin(double(dataPhase));
PhaseCos = cos(double(dataPhase));    
PhaseLaplacian  = PhaseCos.*ifftn_c(lap_opK.* fftn_c(PhaseSin)) - PhaseSin.*ifftn_c(lap_opK.* fftn_c(PhaseCos));  %  

clear PhaseSin PhaseCos
if ~isfield(Params, 'cluster')
    multiWaitbar(textWaitbar, 0.8, 'Color', 'b' );
end

unwrappedPhase  =  real(ifftn_c( inv_lap_opK.* fftn_c(PhaseLaplacian)));   %

PhaseLaplacian = real(PhaseLaplacian(padsize(1)+(1:Ny), padsize(2)+(1:Nx), padsize(3)+(1:Nz)));         % crop back
unwrappedPhase = unwrappedPhase(padsize(1)+(1:Ny), padsize(2)+(1:Nx), padsize(3)+(1:Nz));               % crop back

%% shift to a known constant (find the largest phase region and select the middle point for reference)
                       % use shift in most cases
if RefFlag == 1
    c = - unwrappedPhase(RefVox(1), RefVox(2), RefVox(3));
else
    c = 0;
end
unwrappedPhase = unwrappedPhase + c;        % phi' as in Schofield & Zhu, Optics Letters, 2003

if ~isfield(Params, 'cluster')
    % Update
    multiWaitbar(textWaitbar, 0.99, 'Color', 'b' );
    multiWaitbar(textWaitbar, 'Close');
end
