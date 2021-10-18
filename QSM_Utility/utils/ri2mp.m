function [mm,pp] = ri2mp( rr, ii )
%  Pass in the real and imaginary data and return
%  magnitude and phase 

mm = abs( rr + 1i * ii );
pp = atan2( ii, rr );

