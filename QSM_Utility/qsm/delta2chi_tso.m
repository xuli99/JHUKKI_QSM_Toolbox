function chi = delta2chi_tso(delta, D, thresh)
% chi = delta2chi_tso(delta, D, thresh)
% 
%% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
% QSM TKD method 

warning off all
%% threshold based susceptibility mapping (TSO or TKD)
delta_k = fftshift(fftn(delta));

IndIll = abs(D) < thresh;
Dinv = 1./D;
Dinv(IndIll) = sign(D(IndIll)).*(1/thresh);

chi = real(ifftn(ifftshift(delta_k.*Dinv)));

%% adding correcting factor for TKD
% ref: Scheweser, 2012, MRM, toward online QSM

psf = real(ifftn(ifftshift(D.*Dinv)));
c = psf(1);

chi = chi./c;
