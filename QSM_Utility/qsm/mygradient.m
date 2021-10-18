function [ GA ] = mygradient( A, dim, flag )
% Simple 3D gradient function using difference that satisfy the matrix
% transpose operation
% [ GA ] = mygradient( A, dim, flag )
% Authors: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
% ---------------------------------------

if nargin < 2
    dim = 1;
    flag = 0;   % 0: nontranspose, 1: transpose
elseif nargin < 3
    flag = 0;
end

[n1, n2, n3] = size(A);

if flag == 0        % nontranspose
    
    GA = diff(A, 1, dim);      % Y = diff(X,n,dim)
    switch dim
        case 1
            GA = cat(1, GA, zeros(1, n2, n3));        
        case 2
            GA = cat(2, GA, zeros(n1, 1, n3));        
        case 3
            GA = cat(3, GA, zeros(n1, n2, 1));
        otherwise
            disp('Only for 3D')
    end    
    
    
elseif flag == 1    % transpose

    GA = flipdim(diff(flipdim(A, dim), 1, dim), dim);      % Y = diff(X,n,dim)
    
    switch dim
        case 1
            GA = cat(1, -A(1, :, :), GA);
            GA(end,:,:) = A(end-1,:,:);

        case 2            
            GA = cat(2, -A(:, 1, :), GA);
            GA(:,end,:) = A(:,end-1,:);         
            
        case 3
            GA = cat(3, -A(:, :, 1), GA);
            GA(:,:,end) = A(:,:,end-1);
            
        otherwise
            disp('Only for 3D')
    end        
    
else
    disp('unknown flag')
end

end

