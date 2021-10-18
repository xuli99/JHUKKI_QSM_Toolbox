function [GREPhase, GREMag] = mcpc3Ds(GREMag, GREPhase, Params)
% [GREPhase, GREMag] = mcpc3Ds(GREMag, GREPhase, Params)
% multiple coil phase combination with MCPC-3D-S and ASPIRE
%% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
% 
% Multiple channel phase combination - using MCPC-3D-S or ASPIRE
% Ref. Eckstein & Robinson et al, MRM, 2018
%           Both MCPC-3D-S and ASPIRE needs multiple echoes, i.e. Params.nEchoes > 1
%           ASPIRE needs special acquisition with TE = n*deltaTE, n is interger 
% if data in 6D array, Philips/Siemens: ncol, nrow, nslice, necho, ndyanmics, ncoil
% if data in 8D array, Bruker format: dim1 dim2 dim3 Echo Slice Cycle Repetition Channel,
% output GREPhase, GREMag are 6D array
% Updated 2021-09-23

sizedim = size(GREPhase);
if ndims(GREPhase) == 8
    switch Params.VisuCoreDim
        case 3
            GREPhase = reshape(GREPhase, [sizedim(1:4),sizedim(7:8)]);
            GREMag = reshape(GREMag, [sizedim(1:4),sizedim(7:8)]);
        case 2
            GREPhase = reshape(GREPhase, [sizedim(1:2), sizedim(4:5),sizedim(7:8)]);
            GREPhase = permute(GREPhase, [1,2,4,3,5,6]);    % match 3D
            GREMag = reshape(GREMag, [sizedim(1:2),sizedim(4:5),sizedim(7:8)]);
            GREMag = permute(GREMag, [1,2,4,3,5,6]);    % match 3D
        otherwise
            error('Data is not 2D or 3D format!')
    end
end

sizedim = size(GREPhase); % updated
Nc = Params.coilNum;
Ne = Params.nEchoes;
Nd = Params.nDynamics;
TEs = Params.TEs;         % TEs
dTE = TEs(2) - TEs(1);    % deltaTE

textWaitbar = 'coil combination for phase with MCPC-3D-S/ASPIRE';
if ~isfield(Params, 'cluster')
    multiWaitbar(textWaitbar, 0);
else
    disp(textWaitbar);
end

%% main procedure
% select referecne voxel
if Ne < 2 
    error('MCPC-3D-S/ASPIRE needs to use multiple echoes.')
end   

DeltaThetaKJ = zeros([sizedim(1:3), 1, sizedim(5)]);  % collapse on echo and coil dim

% Estimate delta-theta_k_j for each dynamic, mag weighted sum
for IndDyn = 1:Nd
    DeltaThetaKJ(:,:,:,1,IndDyn) = ...
    angle(sum(GREMag(:,:,:,2,IndDyn,:).*exp(1i*GREPhase(:,:,:,2,IndDyn,:)).*...
        (GREMag(:,:,:,1,IndDyn,:).*exp(-1i*GREPhase(:,:,:,1,IndDyn,:))), 6));    
end

% Estimate phi_0_c, original phase for each coil
% check for ASPIRE eligibility
if mod(TEs(1), dTE) ~= 0
    % do unwrapping on DeltaThetaKJ
    for IndDyn = 1:Nd
        DeltaThetaKJ(:,:,:,1,IndDyn) = phase_unwrap_path_mex(DeltaThetaKJ(:,:,:,1,IndDyn));
    end
end
Phi0 = angle(exp(1i*GREPhase(:,:,:,1,:,:)).*exp(-1i*TEs(1)/dTE*DeltaThetaKJ));

if ~isfield(Params, 'cluster')
    multiWaitbar(textWaitbar, 0.5);
end

% Low-pass filtering/smoothing on Real/Imaginary
Ni = 1;
Nall = Nd*Nc;
for IndDyn = 1:Nd
   for IndChannel = 1:Nc
        temp = GREMag(:,:,:,1,IndDyn,IndChannel).*exp(1i*Phi0(:,:,:,1,IndDyn,IndChannel));
        
        Phi0(:,:,:,1,IndDyn,IndChannel) = angle(smooth3(real(temp), 'gaussian') + ...
                                1i*smooth3(imag(temp), 'gaussian'));     % smooth3 (default box3) or medfilt3

        if ~isfield(Params, 'cluster')
            multiWaitbar(textWaitbar, 0.5+(Ni/Nall)*0.4);
        end
        Ni = Ni + 1;
   end
end

% normalize to VRC and combine, use mag weighted sum
GREPhase = angle(sum(GREMag.*exp(1i.*GREPhase).*GREMag.*exp(-1i.*Phi0), 6));   % remove coil initial phase and do multiple coils combined

if ~isfield(Params, 'cluster')
    multiWaitbar(textWaitbar, 1);
    % Close waitbar
    multiWaitbar( 'CloseAll' );
else
    disp('Done.')
end

