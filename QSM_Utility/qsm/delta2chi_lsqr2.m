function chi_res = delta2chi_lsqr2(deltaB, D, tol, maskErode)
% chi_res = delta2chi_lsqr2(deltaB, D, tol, maskErode)
%
% Input:
% deltaB: relative field shift
% D: convolution kernel, fftshifted
% tol: tolerance for the LSQR solver
% maskErode: brain mask
%
% Authors: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
% 
% Ref: Li et al, NIMG, 2015, 108:111-22. iLSQR method
% ---------------------------------------
% Added a weighting using laplacian of phase


D = ifftshift(D);
N = size(deltaB);

%% get weighting matrix using laplacian of phase

LapPhi_min = 60;            % percentile parameters that may be adjusted
LapPhi_max = 99.9;
   
lap_deltaB = laplacianK(deltaB);
lap_deltaB = lap_deltaB.*maskErode;

a = prctile(abs(lap_deltaB(maskErode>0)), [LapPhi_min, LapPhi_max]);
W0 = (a(2) - abs(lap_deltaB))./(a(2)-a(1));            
Wi = W0;
Wi(abs(lap_deltaB)<a(1)) = 1; 
Wi(abs(lap_deltaB)>a(2)) = 0;
Wi = Wi.*maskErode;

%% set up lsqr solver

maxit = 200;
b = real(ifftn(D.*fftn(Wi.*deltaB)));
b = b(:);

iter_counter = 0;

[chi_res,flag,relres,iter] = lsqr(@afun,b,tol,maxit);
disp(['lsqr flag: ', num2str(flag)]);
disp(['lsqr relres: ', num2str(relres)]);
disp(['lsqr iter: ', num2str(iter)]);

chi_res = reshape(chi_res, N);

function y = afun(x,transp_flag)
   if strcmp(transp_flag,'transp')       % y = A'*x
       x = reshape(x, N); 
       y = real(ifftn(D.*fftn(Wi.*real(ifftn(D.*fftn(x))))));
       y = y(:);             
   elseif strcmp(transp_flag,'notransp') % y = A*x;  
       x = reshape(x, N);
       y = real(ifftn(D.*fftn(Wi.*real(ifftn(D.*fftn(x))))));
       y = y(:);
       iter_counter = iter_counter + 1;
       if mod(iter_counter, 10) == 0
           disp(['iteration: ', num2str(iter_counter), '...']);
       end
   end
end


end
