function [output] = polyfitn3d(input,mask,order,sc)

% ==========================
% HEADER
% ==========================
%
% Function name:    polyfitn3d.m
% Author:           Jonathan Farrell, Xu Li
% Created:          March 3 2010
% Last Modified:    2021-02-09
% Pupose:           fit a nth order polynomial in 3D to input in variable
%                   input, return fitted polynomial
%
% ==========================
% INPUT
% ==========================
%
% input:  3D volume
%  mask:  3D volume with 0 outside the head, 1 inside the head
% order:  integer to specify the maximum order of polynomial terms. 
%
% sc:   sc = 0, or sc = 1.  sc = 1 will result in rescaling and resampling
%
% ==========================
% OUTPUT
% ==========================
%
% output:     3D volume containing the fitted output
%

% ===========================
% CODE
% ===========================

N = size(input);

if (sc)
     % center and scale input
     temp = input;
     temp = temp(mask>0);
     mu1 = mean(temp(:));
     mu2 = std(temp(:));
     input = (input-mu1)./mu2;
end

[X,Y,Z] = meshgrid(1:N(2), 1:N(1), 1:N(3));
x = [X(:) Y(:) Z(:)];

p = polyfitn(x(mask(:)>0,:),input(mask(:)>0),order);
pp = polyvaln(p,x);
output = reshape(pp,N).*mask;

if (sc)
    output = output*mu2 + mu1;
    output = mask.*output;
end
