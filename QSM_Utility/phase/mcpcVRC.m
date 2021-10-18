function [GREPhase, GREMag] = mcpcVRC(GREMag, GREPhase, Params)
% [GREPhase, GREMag] = mcpcVRC(GREMag, GREPhase, Params)
% multiple coil phase combination with VRC
%% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
% 
% Multiple channel phase combination - using virtual reference coil method
% (vrc method), using mcpcc to get the virtual reference first
% Ref. Robinson et al, NMRB, 2016
%        Parker et al, MRM, 2014
% updated 2017-03-15, X.L.
% updated 2019-09-20, X.L. fit JHUKKI_QSM_Toolbox use
% updated 2021-03-11, X.L. added option for 2D scans, added GREMag output
% (matching GREPhase dimensions)

% if data in 6D array, ncol, nrow, nslice, necho, ndyanmics, ncoil
% if data in 8D array, dim1 dim2 dim3 Echo Slice Cycle Repetition Channel
% Updated 2021-06-26, add update for cluster version

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

sizedim = size(GREPhase);
centVox = floor(0.5*(sizedim(1:3)));    
Nc = Params.coilNum;
Ne = Params.nEchoes;
Nd = Params.nDynamics;

textWaitbar = 'coil combination for phase with VRC';
if ~isfield(Params, 'cluster')
    multiWaitbar(textWaitbar, 0);
else
    disp(textWaitbar);
end

%% main procedure
% select referecne voxel
if Ne < 2 && Nd < 2   
    % if single echo and single dynamic, reference to a voxel with good SNR for every coil
    GREMagCoilMin = min(GREMag, [], 6);
    [~, Ind] = max(GREMagCoilMin(:));
    [centVox(1), centVox(2), centVox(3)] = ind2sub(sizedim(1:3), Ind);
end   

Ni = 1;
Nall = Nd*Ne*Nc;
% do global phase correction, reference to center voxel
for IndDyn = 1:Nd
    for IndEcho = 1:Ne
        for IndChannel = 1:Nc
               
            GREPhase(:,:,:, IndEcho, IndDyn, IndChannel) = ...
            angle(exp(1i.*GREPhase(:,:,:,IndEcho, IndDyn, IndChannel)).*...
                (exp(-1i.*GREPhase(centVox(1), centVox(2), centVox(3), IndEcho, IndDyn, IndChannel))));    
            
            if ~isfield(Params, 'cluster')
                multiWaitbar(textWaitbar, (Ni/Nall)*0.4);
            end
            Ni = Ni + 1;
        end
    end
end

% creaste Virtual reference coiil
VRC = sum(GREMag.*exp(1i.*GREPhase), 6)./sum(GREMag, 6);        % virtual reference coil, coil @ 6 dim
VRC(isnan(VRC)) = 0;

% find difference
PhaseDiff = angle(exp(1i.*GREPhase).*repmat(conj(VRC), [1, 1, 1, 1, 1, Nc])); % coil phase difference
if ~isfield(Params, 'cluster')
    multiWaitbar(textWaitbar, 0.5);
end

% Low-pass filtering/smoothing
Ni = 1;
for IndDyn = 1:Nd
    for IndEcho = 1:Ne
       for IndChannel = 1:Nc
            PhaseDiff(:,:,:,IndEcho,IndDyn,IndChannel) = smooth3(PhaseDiff(:,:,:,IndEcho,IndDyn,IndChannel));     % smooth3 (default box3) or medfilt3
            % PhaseDiff(:,:,:,IndEcho,IndDyn,IndChannel) = medfilt3(PhaseDiff(:,:,:,IndEcho,IndDyn,IndChannel));     % smooth3 or medfilt3
            if ~isfield(Params, 'cluster')
                multiWaitbar(textWaitbar, 0.5+(Ni/Nall)*0.4);
            end
            Ni = Ni + 1;
       end
    end
end

% normalize to VRC and combine 
GREPhase = angle(exp(1i.*GREPhase).*exp(-1i.*PhaseDiff));        % normalized to VRC    
GREPhase = angle(sum(GREMag.*exp(1i.*GREPhase), 6));             % multiple coils combined

if ~isfield(Params, 'cluster')
    multiWaitbar(textWaitbar, 1);
    % Close waitbar
    multiWaitbar( 'CloseAll' );
else
    disp('Done.')
end

