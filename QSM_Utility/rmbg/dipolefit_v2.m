function [Blocal, chi_source] = dipolefit_v2(data, mask, weight, Params, tol, itermax, centerswitch)
% [Blocal, chi_source] = dipolefit_v2(data, mask, weight, Params, tol, itermax, centerswitch)
%% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
%
% Purpose:
% dipole field fitting for removing the dipole field caused by large
% air-tissue interface
% 
% INPUT: data: data to fit, such as the wassr or GRE raw frequency map (unwrapped)
%         mask: mask define the ROI, such as the brain
%         weight: weighting matrix, generally use the modulus image
%         Params: imaging parameters
%         tol: tolerance for the iterative solver
% OUTPUT: Blocal: local field after subtracting the background
%         chi_source: fitted dipole sources outside of ROIs
% Other function used:
%         conv_kernel_rot_c0.m, cgsolve.m
%
% Ref: T. Liu et al. NMR Biomed 2011;24(9):1129-36
% 
% Updated 2013-04-01 Xu Li
% Updated 2018-01-04, X.L. adjusted to match the original PDF method
% Updated 2020-05-13, X.L. added itermax and centerswitch as input parameters


%% BODY
warning off all

if (nargin<5)
    tol = 0.1;
    itermax = 30;
    centerswitch = true;
elseif (nargin<6)
    itermax = 30;
    centerswitch = true;
end

if isempty(tol)
    tol = 0.1;
end

%% Mask, normalize
mask = mask > 0;        % change to logical
N = size(data);         % original size of volume

weight = weight.*mask;
temp = weight(mask);
weight = weight./mean(temp(:));     % normalize weighting

%% zero padding
if centerswitch
    d1 = max(max(mask,[],2),[],3);
    d1first = find(d1,1,'first');
    d1last = find(d1,1,'last');

    d2 = max(max(mask,[],1),[],3);
    d2first = find(d2,1,'first');
    d2last = find(d2,1,'last');

    d3 = max(max(mask,[],1),[],2);
    d3first = find(d3,1,'first');
    d3last = find(d3,1,'last');
else
    d1first = 1; d1last = N(1);
    d2first = 1; d2last = N(2);
    d3first = 1; d3last = N(3);
end
    
% check odd dimension
[d1first, d1last] = makeeven(d1first, d1last);
[d2first, d2last] = makeeven(d2first, d2last);
[d3first, d3last] = makeeven(d3first, d3last);

while d1last > N(1)
    d1last = d1last - 2;
end

while d2last > N(2)
    d2last = d2last - 2;
end

while d3last > N(3)
    d3last = d3last - 2;
end

%%% center piece
data=data(d1first:d1last,d2first:d2last, d3first:d3last);
weight=weight(d1first:d1last,d2first:d2last, d3first:d3last);
mask=mask(d1first:d1last,d2first:d2last, d3first:d3last);

padsize = [32, 32, 32];

%%% pad
data = padarray(data, padsize, 0,'both');
weight = padarray(weight, padsize, 0,'both');
mask = padarray(mask, padsize, 0,'both');

%%% dipole kernel
Params.sizeVol = size(data);
Params.fov = Params.sizeVol.*Params.voxSize;

D = conv_kernel_rot_c0(Params, Params.TAng);
D = fftshift(D);

%%%%%% start the PDF method %%%%%
W_std = weight;
W_var = weight.^2;

rhs_temp = real(ifftn(D.*fftn(W_var.*(data) )));
b = rhs_temp( mask(:) == 0);

% set erance level and maximum iteration allowed
E_nl = real(ifftn(D.*fftn(W_std)));

A = @(xx)(Ax(W_var,D,mask,xx) );
cg_tol = tol*norm(E_nl(mask(:) == 0))/norm(b(:));   % tolerance estimate

[chi0, res, num_iter] = cgsolve(A, b, cg_tol, itermax, 0);
fprintf('iterations: %f ; relres: %f\n', num_iter, res);

chi = zeros(size(D));
chi(mask(:) == 0) = chi0(1:end);
chi(mask(:) > 0) = 0;

% background dipole field
p_dipole = real(ifftn(D.*fftn(chi)));
p_final = (data-p_dipole).*mask;  

%% remove zero pad
Blocal = zeros(N);
Blocal(d1first:d1last,d2first:d2last, d3first:d3last) = ...
    p_final((padsize(1)+1):(d1last-d1first+padsize(1)+1), ...
                (padsize(2)+1):(d2last-d2first+padsize(2)+1), ...
                    (padsize(3)+1):(d3last-d3first+padsize(3)+1));

chi_source = zeros(N);
chi_source(d1first:d1last,d2first:d2last, d3first:d3last) = ...
    chi((padsize(1)+1):(d1last-d1first+padsize(1)+1), ...
                (padsize(2)+1):(d2last-d2first+padsize(2)+1), ...
                    (padsize(3)+1):(d3last-d3first+padsize(3)+1));

%% internal function
    function y = Ax(W,D,Mask,xx)
        x = zeros(size(D));
        x(Mask(:) == 0) = xx(1:end);
        x(Mask(:) == 1) = 0;            % source are outside of mask only, Mbar*x

        AtA = real(ifftn(D.* fftn(W.*real(ifftn(D.* fftn(x) )))));
        y = AtA( Mask(:) == 0);
    end

    function [first, last] = makeeven(first, last)
        if mod(last - first + 1, 2) > 0
            if first > 1
                first = first -1;
            else
                last = last + 1;
            end
        end
    end
end
