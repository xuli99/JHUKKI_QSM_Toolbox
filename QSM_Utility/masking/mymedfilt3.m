function B = mymedfilt3(A, filterSize, padopt)
% my 3D median filter without mex
if nargin < 2
    filterSize = [3, 3, 3];
    padopt = 0;
elseif nargin < 3
    validateattributes(filterSize, ...
    {'numeric'}, ...
    {'real','positive','odd','integer','vector','numel',3}, ...
    mfilename, ...
    '[m n p]');
    padopt = 0;
end

radius = (filterSize - 1) / 2;
N = size(A);
Ap = padarray(A, radius, padopt);

B = zeros(N);
Index1 = (1:N(1)) + radius(1);
Index2 = (1:N(2)) + radius(2);
Index3 = (1:N(3)) + radius(3);

textWaitbar = 'generating unreliable phase mask ...';
multiWaitbar(textWaitbar, 0, 'CanCancel', 'On' );
% main loop
for ii = Index1
    for jj = Index2
        for kk = Index3            
            box = Ap(ii-radius(1):ii+radius(1), ...
                      jj - radius(2):jj+radius(2), ...
                      kk - radius(3):kk+radius(3));
            B(ii-radius(1), jj-radius(2), kk-radius(3)) = median(box(:));
        end
    end
    hasCanceled = multiWaitbar(textWaitbar, (ii/length(Index1)));
    HandleStopReconstruction;
end

