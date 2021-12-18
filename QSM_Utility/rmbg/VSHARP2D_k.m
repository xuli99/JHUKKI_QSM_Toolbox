function [phase_local, background, mask_eval] = VSHARP2D_k(GREPhaseSE, BrainMask, radiusArray, SMV_thres, Params, handles)
% [phase_local, background, mask_eval] = VSHARP2D_k(GREPhaseSE, BrainMask, radiusArray, SMV_thres, Params, handles)
% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
%
% Ref: Schweser et al. 2011, NIMG, Wu, et al., 2011, MRM, Bilgic, et. al NMRB 2017
%
% NOTE: k-space based VSHARP implementation, 2017-04
% radiusArray is max_radius:step_size_radius:min_radius, does not need to match resolution in the phase data
% Do k-space based V-SHARP for each echo, does not average
% 1. generate k-space based SMV kernals and masks for different radius
% 2. calculate SMV
% 3. do inverse with kernel with maximum radius
% 2D version based on 2D SMV
% updated 2019-07-09
% Updated 2021-06-27 X.L., cluster version
% Updated 2021-12-10 X.L., bug fix

warning off all

if radiusArray(1)<radiusArray(end)
    radiusArray = sort(radiusArray, 'descend');     % sort radiusArray in descending order
end

% Start the clock
textWaitbar = ['Performing 2D V-SHARP on ' num2str(length(Params.echoNums)) ' echoes'];
if ~isfield(Params, 'cluster')
    multiWaitbar(textWaitbar, 0, 'Color', 'b', 'CanCancel', 'On' );
else
    disp(textWaitbar);
end

% Cut only what we want
if max(Params.echoNums) < size(GREPhaseSE, 4)
    GREPhaseData = GREPhaseSE(:,:,:,Params.echoNums);                   % bug fix 2021
else
    GREPhaseData = GREPhaseSE;
end

phase_local = zeros(size(GREPhaseData), class(GREPhaseData));       % final output

% Create k-space kernel with different radius, 2D 
N = Params.sizeVol;                                
[Y,X] = meshgrid(-floor(N(2)/2):ceil(N(2)/2-1),-floor(N(1)/2):ceil(N(1)/2-1));

X = X * Params.voxSize(1);
Y = Y * Params.voxSize(2);

num_kernel = length(radiusArray);           % number of kernels
unrely_tol = 1e-3;                          % tol for unreliable boundary voxels
SMV_kernel = zeros([N(1:2), num_kernel]);    
mask_Sharp = zeros([N, num_kernel]);        % 3D, total mask is still 3D
SMV_inv_kernel = zeros(N(1:2));              % tsvd

for k = 1:num_kernel
    SMV = gen_SMVkernel_voxel_scaled( X, Y, radiusArray(k));     % kernel in k-space
    SMV_kernel(:,:,k) = SMV;
    if k == 1                               % with max_radius
        SMV_inv_kernel( abs(SMV) > SMV_thres ) = 1 ./ SMV( abs(SMV) > SMV_thres );
    end
end


for sliceii = 1:N(3)    
    mask_prev = zeros(N(1:2));      % 2D
    
    for k = 1:num_kernel
        % SMV = gen_SMVkernel_voxel_scaled( X, Y, radiusArray(k));     % kernel in k-space        
        mask_rely_iik = gen_SMVMask_new(SMV_kernel(:,:,k), BrainMask(:,:,sliceii), unrely_tol);       % reliable mask

        if sum(mask_rely_iik(:)) == 0
            continue
        end

        mask_Sharp(:,:,sliceii,k) = (mask_rely_iik-mask_prev);
        mask_prev = mask_rely_iik;

    end
    if ~isfield(Params, 'cluster')
        hasCanceled = multiWaitbar(textWaitbar, 0.5*(sliceii/N(3)));
        HandleStopReconstruction;    
    end
end
mask_eval = sum(mask_Sharp, 4) > 0;         %   final mask for evalualtion 


%% Do the loop for all echoes
procNechos = size(GREPhaseData,4);
for selectedEcho = 1:procNechos

    for sliceii = 1:N(3)
        
        phase_local(:,:,sliceii,selectedEcho) = fftn(GREPhaseData(:,:,sliceii,selectedEcho));
        phase_Sharp = zeros(N(1:2));

        % convolve with V-SHARP kernels
        for k = 1:num_kernel
            phase_Sharp = phase_Sharp +  mask_Sharp(:,:,sliceii,k).* ifftn(SMV_kernel(:,:,k) .*phase_local(:,:,sliceii,selectedEcho));
        end      

        % deconvolution
        phase_local(:,:,sliceii,selectedEcho) = mask_eval(:,:,sliceii).*ifftn(SMV_inv_kernel.*fftn(phase_Sharp));
        
        if ~isfield(Params, 'cluster')
            hasCanceled = multiWaitbar(textWaitbar, 0.5+0.5*(selectedEcho/procNechos)*(sliceii/N(3)));
            HandleStopReconstruction;
        end
    end
   
end

% Calculate background
background = (GREPhaseData - phase_local).*repmat(mask_eval, [1 1 1 size(GREPhaseData, 4)]);

% Average everything
phase_local = mean(phase_local,4);
background = mean(background,4);

%% subfunctions
function SMV = gen_SMVkernel_voxel_scaled( X, Y, smv_rad)
  
    smv = (X.^2 + Y.^2) <= smv_rad^2;
    smv = smv / sum(smv(:));                    % normalized 

    smv_kernel = zeros(size(X));
    smv_kernel(1+floor(end/2),1+floor(end/2)) = 1;   
    smv_kernel = smv_kernel - smv;              % delta - smv

    SMV = fftn(fftshift(smv_kernel));           % kernel in k-space 

end

function mask_rely = gen_SMVMask_new( SMV, chi_mask, unrely_tol)
  
    mask_unrely = ifftn(SMV .* fftn(chi_mask));             
    mask_unrely = abs(mask_unrely) > unrely_tol;            

    mask_rely = chi_mask .* (chi_mask - mask_unrely);       

end


end
