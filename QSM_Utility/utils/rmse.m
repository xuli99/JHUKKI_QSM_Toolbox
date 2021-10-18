function [ rmse ] = rmse( in, true, use_abs )
%RMSE Summary of this function goes here
%   Detailed explanation goes here
%   output percentage
%   updated by X.L. 202003

if nargin < 3
	use_abs = 0;
end

if (nargin == 3) && (use_abs == 1)
    rmse = 100 * norm(abs(in(:)) - abs(true(:))) / norm(abs(true(:)));
else
    rmse = 100 * norm(in(:) - true(:)) / norm(true(:));
end

end

