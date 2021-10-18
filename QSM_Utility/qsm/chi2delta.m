function delta = chi2delta(chi, D)
% delta = chi2delta(chi, D)
% 
% susceptibiltiy to field change with dipoe kernel D

%% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu

%% convolution
chi_k = fftshift(fftn(chi));

delta = real(ifftn(ifftshift(chi_k.*D)));