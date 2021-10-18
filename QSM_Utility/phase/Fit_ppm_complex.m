% Nonlinear fitting for Frequencey map estimation
%   [p1, dp1, relres, p0]=Fit_ppm_complex(M)
%    
%   output
%   p1 - field map, may need further spatial unwrapping
%   dp1 - a priori error estimate (N_std), noise estimation.
%   relres - relative residual
%   p0 - initla phase
%
%   input
%   M - a multi-echo and could be a multi-channel dataset
%       echo needs to be the 4th dimension
%       channel needs to be the 5th dimension
%
%   When using the code, please cite 
%   T. Liu et al. MRM 2013;69(2):467-76, small difference exists as compared to this paper
%   B. Kressler et al. IEEE TMI 2010;29(2):273-81
%   de Rochefort et al. MRM 2008;60(4):1003-1009
%
%   Adapted from a linear fitting created by Ludovic de Rochefort
%   Modified by Tian Liu on 2011.06.01
%   Last modified by Tian Liu on 2013.07.23
%   commented by X.L., 2017.12.11


function [p1, dp1, relres, p0]=Fit_ppm_complex(M)

if size(M,5)>1
% combine multiple coils together, assuming the coil is the fifth dimension
    M = sum(M.*conj( repmat(M(:,:,:,1,:),[1 1 1 size(M,4) 1])),5);  % sum(M.*conj(M1)) over coil, M1 is ref M@TE1
    M = sqrt(abs(M)).*exp(1i*angle(M));                             % Mag --> sqrt(Mi.*M1), Phase-->(Phi_i - Phi1)
end

% M = conj(M);            % conj, change sign of phase, change at later stage in our code
s0 = size(M);             % 4D if multi-ehcho, 3D if single-echo, original size
L_s0 = length(s0);        % dimension size, 4/3 for multi/single echo, or last dimension
nechos = size(M,L_s0);    % number of echoes Necho

M = reshape(M,[prod(s0(1:L_s0-1)),s0(L_s0)]); % reshape to NxNecho
s = size(M);                % new size

Y = angle(M(:,1:min(3,nechos)));  % phase of the first 3 echoes for initial estimate
c = ((Y(:,2)-Y(:,1)));            % phase difference between first 2 echoes, in radian, in the range of [-2*pi, 2*pi]
[~, ind] = min([abs(c-2*pi),abs(c),abs(c+2*pi)],[],2); % c,c+2pi, c-2pi find min abs
c(ind==1) = c(ind==1) - 2*pi;       % if c is in [pi, 2pi], change to [-pi, 0]
c(ind==3) = c(ind==3) + 2*pi;       % ajust c to -pi to pi, similar to time unwrapping Y(:,1:2) (for initial frequency estimator)

for n=1:min(2,nechos-1)
    cd=((Y(:,n+1)-Y(:,n)))-c;                   
    Y(cd<-pi,(n+1):end)=Y(cd<-pi,n+1:end)+2*pi; 
    Y(cd>pi,(n+1):end)=Y(cd>pi,n+1:end)-2*pi;
end
A = [1  0 ;1 1; 1 2 ];
ip = A\Y(:,1:3)';       % least square fitting
p0 = ip(1,:)';          % initial guess of phi0, Nx1
p1 = ip(2,:)';          % initial guess of normalized freq, (2*pi*dTE)

dp1 = p1;               % 
tol = norm(p1(:))*1e-4; % tol, 1e-4*norm(p1), instead of p0/p1
iter = 0;
max_iter = 30;          % instead of 10

% weigthed least square
% calculation of WA'*WA, note that norm(W)..^2 = norm(M).^2
v1=ones(1,nechos);
v2=(0:(nechos-1));
a11=sum(abs(M).^2.*(ones(s(1),1)*(v1.^2)),2);   % WtW, W.^2 is equivalent to M.^2
a12=sum(abs(M).^2.*(ones(s(1),1)*(v1.*v2)),2);  % 
a22=sum(abs(M).^2.*(ones(s(1),1)*(v2.^2)),2);   % 
% inversion using QR, similar as pinv
d=a11.*a22-a12.^2;  % determinant
ai11=a22./d;
ai12=-a12./d;
ai22=a11./d;

while ((norm(dp1)>tol) &&(iter<max_iter))
    iter = iter+1;
    W = abs(M).*exp(1i*(p0*v1 + p1*v2) );       % estimate of M

    % projection
    pr1=sum(conj(1i*W).*(ones(s(1),1)*v1).*(M-W),2);    % RHS, b1, At*conj(iW)*(M-W), sum(NxNe, 2)
    pr2=sum(conj(1i*W).*(ones(s(1),1)*v2).*(M-W),2);    % RHS, b2

    dp0=real(ai11.*pr1+ai12.*pr2);          % update on phi0
    dp1=real(ai12.*pr1+ai22.*pr2);          % update on freq
    dp1(isnan(dp1))=0;
    dp0(isnan(dp0))=0;
    
    %update
    p1 = p1+dp1;
    p0 = p0+dp0;
    

end

% error propagation
dp1=sqrt(ai22);
dp1(isnan(dp1)) = 0;
dp1(isinf(dp1)) = 0;

% relative residual
res = M - abs(M).*exp(1i*(p0*v1 + p1*v2) );
relres = sum(abs(res).^2,2)./sum(abs(M).^2,2);
relres(isnan(relres)) = 0;

p1(p1>pi)=mod(p1(p1>pi)+pi,2*pi)-pi;        % correct shift
p1(p1<-pi)=mod(p1(p1<-pi)+pi,2*pi)-pi;

p0=reshape(p0,s0(1:L_s0-1));
p1=reshape(p1,s0(1:L_s0-1));
dp1=reshape(dp1,s0(1:L_s0-1));
relres = reshape(relres,s0(1:L_s0-1));
    

