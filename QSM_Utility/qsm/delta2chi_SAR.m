function [chi_res, chi_SA] = delta2chi_SAR(deltaB, D, chi_res0, maskErode, thresh_SAR, chi_ap)     
% function [chi_res, chi_SA] = delta2chi_SAR(deltaB, D, chi_res0, maskErode, thresh_SAR, chi_ap)  
% QSM with streaking artifact removal (SAR)
% Input:
% deltaB: relative field shift
% D: convolution kernel
% chi_res0: initial estimate 
% maskErode: brain mask
%
% Authors: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
% 
% Ref: Li et al, NIMG, 2015, 108:111-22. iLSQR method
% ---------------------------------------

N = size(deltaB);
VoxNum = numel(deltaB);

D = ifftshift(D);
Mic = abs(D) < thresh_SAR;              % ill conditioned k-space mask

% ---------------------------------------
LapPhi_min = 40;            % percentile parameters that may be adjusted
LapPhi_max = 60;
   
lap_deltaB = laplacianK(deltaB);
lap_deltaB = lap_deltaB.*maskErode;

a = prctile(abs(lap_deltaB(maskErode>0)), [LapPhi_min, LapPhi_max]);
W0 = (a(2) - abs(lap_deltaB))./(a(2)-a(1));            
Wi = W0;
Wi(abs(lap_deltaB)<a(1)) = 1;        % low gradient
Wi(abs(lap_deltaB)>a(2)) = 0;        % high gradient

% Calculate weighting from chi_ap (some SA inside)
Wtx = abs(mygradient(chi_ap, 1));   % spatial gradient
Wtxth = prctile(Wtx(maskErode > 0), [40, 60]);
Wtxf = (Wtxth(2) - Wtx)./(Wtxth(2) - Wtxth(1));
Wtxf(Wtx < Wtxth(1)) = 1;
Wtxf(Wtx > Wtxth(2)) = 0;
Wtxf = Wtxf.*Wi;

Wty = abs(mygradient(chi_ap, 2));
Wtyth = prctile(Wty(maskErode > 0), [40, 60]);
Wtyf = (Wtyth(2) - Wty)./(Wtyth(2) - Wtyth(1));
Wtyf(Wty < Wtyth(1)) = 1;
Wtyf(Wty > Wtyth(2)) = 0;  
Wtyf = Wtyf.*Wi;

Wtz = abs(mygradient(chi_ap, 3));
Wtzth= prctile(Wtz(maskErode > 0), [40, 60]);
Wtzf = (Wtzth(2) - Wtz)./(Wtzth(2) - Wtzth(1));
Wtzf(Wtz < Wtzth(1)) = 1;
Wtzf(Wtz > Wtzth(2)) = 0; 
Wtzf = Wtzf.*Wi;

% iterative solve for chi_SA using LSQR    
    % b
    tempx = mygradient(chi_res0, 1);
    tempy = mygradient(chi_res0, 2);
    tempz = mygradient(chi_res0, 3);
    b = cat(1, Wtxf(:).*tempx(:), Wtyf(:).*tempy(:), Wtzf(:).*tempz(:));   % if A is 3N x N

    tol = 1e-2;
    maxit = 200;    
    
    iter_counter = 0;
    disp('LSQR for estimating chi_SA ...')
    [chi_SA, flag, relres, iter, resvec] = lsqr(@afun,b(:),tol,maxit);
    disp(['flag; relres; iter = ', num2str([flag, relres, iter])]);
    disp('Done.')

    chi_SA = reshape(chi_SA, N);
        
    % final QSM 
    chi_res = (chi_res0 - chi_SA).*maskErode;
    
% ---------------------------------------        
    function y = afun(x,transp_flag)
       if strcmp(transp_flag,'transp')      % y = A'*x

           x1 = reshape(x(1:VoxNum), N);
           x2 = reshape(x(VoxNum+1:2*VoxNum), N);
           x3 = reshape(x(2*VoxNum+1:3*VoxNum), N);
           temp = mygradient(Wtxf.*x1, 1, 1) + mygradient(Wtyf.*x2, 2, 1)+ mygradient(Wtzf.*x3, 3, 1);

           y = real(ifftn(fftn(temp).*Mic));
           y = y(:);
           
       elseif strcmp(transp_flag,'notransp') % y = A*x;  

           x = reshape(x, N);
           temp = real(ifftn(fftn(x).*Mic));
           tempx = mygradient(temp, 1);
           tempy = mygradient(temp, 2);
           tempz = mygradient(temp, 3);
           
           y = cat(1, Wtxf(:).*tempx(:), Wtyf(:).*tempy(:), Wtzf(:).*tempz(:)); 
           y = y(:);
           iter_counter = iter_counter + 1;
           if mod(iter_counter, 10) == 0
               disp(['iteration: ', num2str(iter_counter), '...']);
           end
       end
    end

end
