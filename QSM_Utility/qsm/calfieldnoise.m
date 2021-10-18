% Noise standard deviation calculation
%   noise_level = calfieldnoise(iField, Mask)
% 
%   output
%   noise_level - the noise standard devivation of complex MR signal
%
%   input
%   iField - the complex MR image
%   Mask - the a region of interest that is not pure noise
%
%   Created by Tian Liu in 20013.07.24
%   Last modified by Tian Liu on 2013.07.24
%   modified by X.L. 2015-04-20 to fit the QSMToolbox dataformat and procedure

function noise_level = calfieldnoise(GREMag, GREPhaseRaw1, Mask)

expected_SNR = 40;
iMag = sqrt(sum(abs(GREMag).^2,4));
iField1 = GREMag(:,:,:,1).*exp(1i*GREPhaseRaw1).*(~Mask).*(iMag<max(iMag(:))/expected_SNR).*(iMag>0);
noise_level = std(real(iField1(iField1~=0))); 