function [slope,intercept] = find_slopeintercept_phasevstime_fast(dataPhase,TEs)
% [slope,intercept] = find_slopeintercept_phasevstime_fast(dataPhase,TEs)
%% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
% updated 2015-04, X.L.
% updated 2017-02, X.L.
% updated 2021, X.L.

dim = ndims(dataPhase);
siz = size(dataPhase);

sizVol = siz(1:(dim-1));      % assuming the last dimension is echoes

dataPhase = reshape(dataPhase, [prod(sizVol), siz(dim)]);     % reshape to NrXNe

N = length(TEs);                    % number of echoes
if N ~= siz(dim)
    error('echo number does not match');
end

if (abs(sum(diff(TEs(:), 2, 1))) < 100*eps) && (abs(mean(diff(TEs(:), 1, 1)) - TEs(1)) < 100*eps)
    A = [(1:N)', ones(size(TEs(:)))];           % NeX2
    scalingfactor = mean(diff(TEs(:), 1, 1));
else
    A = [TEs(:), ones(size(TEs(:)))];   
    scalingfactor = 1;
end

disp('fitting frequency ...')

ip = A\dataPhase';                  % dataPhase' is NeXNr; ip is 2XNr
slope = ip(1,:);
intercept = ip(2,:);

slope = slope./scalingfactor;

% reshape
slope = reshape(slope, sizVol);
intercept = reshape(intercept, sizVol);

disp('Done')


