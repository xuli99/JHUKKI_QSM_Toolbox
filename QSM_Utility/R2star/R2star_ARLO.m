function R2starmap = R2star_ARLO(GREMag, TEs, BrainMask)
% R2starmap = R2star_ARLO(GREMag, TEs, BrainMask)
% 
%% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
% 
% calculat R2* maps from GRE Magnitude, using the ARLO method
% Ref: Pei et al, MRM 73:843 2015

classname = class(GREMag);
N = size(GREMag);

if N(end) ~= length(TEs)
    error('data size does not agree.');
end

if nargin < 3
    BrainMask =  GREMag(:,:,:,1) > 1;           % default BrainMask, amplitude > 1
end
BrainMask = BrainMask > 0;                      % to logical

%% 
TEs = TEs(:);                                   % in unit of sec
dTE = TEs(2) - TEs(1);                          % needs to be uniform in sec

% accumulating var
Si2 = zeros(N(1:3), classname);
SiDi = zeros(N(1:3), classname);
Di2 = zeros(N(1:3), classname);

%% Calculate components related to the AR model
for ii = 1:length(TEs)-2
    
    Si = (dTE/3*(GREMag(:,:,:,ii) + 4*GREMag(:,:,:,ii+1) + GREMag(:,:,:,ii+2)));
    Di = (GREMag(:,:,:,ii) - GREMag(:,:,:,ii+2));
    
    Si2 = Si2 + Si.^2;
    SiDi = SiDi + Si.*Di;
    Di2 = Di2 + Di.^2;

end

R2starmap = (dTE/3*Di2 + SiDi)./(Si2 + dTE/3*SiDi);

R2starmap(isnan(R2starmap)) = 0;
R2starmap(isinf(R2starmap)) = 0;

R2starmap = R2starmap.*BrainMask;
