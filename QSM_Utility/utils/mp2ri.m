function [rr,ii] = mp2ri( mm, pp )
%  Pass in the magnitude and phase and return
% the real and imaginary data

rr = mm .* cos( pp );
ii = mm .* sin( pp );

