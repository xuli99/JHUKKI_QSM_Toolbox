function output  = laplacianK(input)
% do numerical laplacian in fourier space using real kernel

N = size(input);

if N(1) < 256
    padsize = [16,16,16];
else     % high res
    padsize = [64, 64, 64];
end

dimsOdd = mod(size(input), 2);
input = padarray(input, dimsOdd, 'replicate', 'post');  % padding to even dim first     
input = padarray(input, padsize, 'replicate');

Npad = size(input); 

ksize = [3, 3, 3];               
khsize = (ksize-1)/2;

kernelcore = [];
kernelcore(:,:,1) = [0 0 0; 0 1 0; 0 0 0];
kernelcore(:,:,2) = [0 1 0; 1 -6 1; 0 1 0];
kernelcore(:,:,3) = [0 0 0; 0 1 0; 0 0 0];

Kernel = zeros(Npad);
Kernel( 1+Npad(1)/2 - khsize(1) : 1+Npad(1)/2 + khsize(1), 1+Npad(2)/2 - khsize(2) : 1+Npad(2)/2 + khsize(2), ...
    1+Npad(3)/2 - khsize(3) : 1+Npad(3)/2 + khsize(3) ) = kernelcore;

% laplacian kernel in k space
lap_opK = fftn(fftshift(Kernel));                                     % make real delta function, only do fftn

% laplacian of the phase to get weighting matrix   
output = ifftn(fftn(input).*lap_opK);

output = output(padsize(1)+(1:N(1)), padsize(2)+(1:N(2)), padsize(3)+(1:N(3)));