function mask_unrelyPhase = create_mask_unrelyPhase2(phasedata, mask, Params, thresh)
%% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
%
% This function creates unreliable phase mask for QSM
% Ref. Wang et. al., bioRxiv preprint doi: https://doi.org/10.1101/2021.06.28.450248

mask_unrelyPhase = zeros(size(phasedata(:,:,:,1)));  % 3D mask
N_mask = size(mask_unrelyPhase);

% making a spherical kernel
N_kernel = [5,5,5];
kernel_radius = 2;         % in mm
[Y,X,Z] = meshgrid(-floor(N_kernel(2)/2):ceil(N_kernel(2)/2-1),-floor(N_kernel(1)/2):ceil(N_kernel(1)/2-1),-floor(N_kernel(3)/2):ceil(N_kernel(3)/2-1));
X = X * Params.voxSize(1);
Y = Y * Params.voxSize(2);
Z = Z * Params.voxSize(3);
kernel = (X.^2 + Y.^2 + Z.^2) <= kernel_radius^2;
kernel = kernel / sum(kernel(:));

cmplx_data = exp(1i.*phasedata);
norm_factor = convn(mask, kernel, 'same').*mask;
phase_reliability = abs(convn(cmplx_data, kernel, 'same')).*mask;

% correct for bounday voxels
boundvox = (norm_factor < 0.999).*mask;
for kk = 1:N_mask(3)
    for jj = 1:N_mask(2)
        for ii = 1:N_mask(1)
            % update dataConvRho inside the mask
            if boundvox(ii,jj,kk) > 0
                 xmin = max([ii-kernel_radius, 1]); xmax = min([ii+kernel_radius, N_mask(1)]);
                 ymin = max([jj-kernel_radius, 1]); ymax = min([jj+kernel_radius, N_mask(2)]);
                 zmin = max([kk-kernel_radius, 1]); zmax = min([kk+kernel_radius, N_mask(3)]);
                 
                 tempdata = cmplx_data(xmin:xmax, ymin:ymax, zmin:zmax);
                 tempmask = mask(xmin:xmax, ymin:ymax, zmin:zmax);                                           
                 phase_reliability(ii, jj, kk) = abs(mean(tempdata(tempmask>0)));
            end
        end
    end
end

mask_unrelyPhase = (phase_reliability < thresh).*mask;