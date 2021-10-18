function h = mimage(A, range_min, range_max)

%% display montage image of 3D volume data

[ImageHeight, ImageWidth, ImageSlices, ImageVolume] = size(A);

if nargin < 1
    disp('incorrect input number');
elseif (nargin == 1) && (ImageVolume ~= 3)
    range_min = min(A(:)).*1;
    range_max = max(A(:)).*1;    
end
    
if ImageVolume == 1     % gray scale
    A_I = zeros(ImageHeight,ImageWidth, 1, ImageSlices);
    for i = 1:ImageSlices
        A_I(:,:,1,i) = A(:,:,i);
    end

    h = figure;
    warning off
    montage(A_I, 'DisplayRange', [range_min, range_max]);
    warning on
    
elseif ImageVolume == 3 % true color

    A_I = permute(A, [1,2,4,3]); 

    h = figure;
    warning off
    if size(A_I, 4) > 1     % compatibility with matlab R2018
        montage(A_I);
    else
        imshow(A_I);
    end
    warning on

elseif ImageSlices == 1 % dynamics

    h = figure;
    warning off
    montage(A, 'DisplayRange', [range_min, range_max]);
    warning on
    
else
    disp('unknow image type...')
end

