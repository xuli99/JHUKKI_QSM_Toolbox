function GREPhase = mcpcc(GREMag, GREPhase, Params)
% GREPhase = mcpcc(GREMag, GREPhase, Params)
% 
%% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
% 
% Multiple channel phase combination - constant
% updated 2017-03-15, X.L.

%% get virtual reference phase using mcpcc method first
centVox = floor(0.5*(Params.sizeVol));    
Nc = Params.nchannel;
Ne = Params.nEchoes;

if Params.nEchoes > 1   % multiecho data X*Y*Z*NE*Nc    
    % do global phaes correction, reference to center voxel    
    for IndEcho = 1:Ne
        for IndChannel = 1:Nc
            GREPhase(:,:,:,IndEcho,IndChannel) = ...
            angle(exp(1i.*GREPhase(:,:,:,IndEcho,IndChannel)).*...
                (exp(-1i.*GREPhase(centVox(1), centVox(2), centVox(3), IndEcho, IndChannel))));    
                
        end
    end
    
    VRC = sum(GREMag.*exp(1i.*GREPhase), 5)./sum(GREMag, 5);    % virtual reference coil
    VRC(isnan(VRC)) = 0;
               
else                    % single echo data
    % reference to a voxel with good SNR for every coil
    GREMagCoilMin = min(GREMag, [], 4);
    [~, Ind] = max(GREMagCoilMin(:));
    [centVox(1), centVox(2), centVox(3)] = ind2sub(Params.sizeVol, Ind);
    
    % do global phaes correction, reference to center voxel    
    for IndChannel = 1:Nc
        GREPhase(:,:,:,IndChannel) = ...
        angle(exp(1i.*GREPhase(:,:,:,IndChannel)).*(exp(-1i.*GREPhase(centVox(1), centVox(2), centVox(3), IndChannel))));
    end
    
    VRC = sum(GREMag.*exp(1i.*GREPhase), 4)./sum(GREMag, 4);    % virtual reference coil
    VRC(isnan(VRC)) = 0;
    
end

GREPhase = angle(VRC);