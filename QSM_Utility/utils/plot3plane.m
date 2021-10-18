function [h1, h2, h3] = plot3plane(data, sag_ind, cor_ind, axi_ind, cmin, cmax, lineflag, plotflag)      
%     

%% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
% 
% plot 3 planes of the 3D dataset
% lineflag: 0: no line, 1: with lines to mark the slice locations
% plotflag: 1: use imagesc: 2: using mimage

% updated 01/02/2013
% updated 05/05/2017 for new format for v2.8

dimnum = length(size(data));
[ImageHeight, ImageWidth, ImageSlices] = size(data);

if nargin < 2
    sag_ind = floor(ImageHeight/2) + 1;
    cor_ind = floor(ImageWidth/2) + 1;
    axi_ind = floor(ImagesSlices/2) + 1;
    cmin = min(data(:));
    cmax = max(data(:));
    lineflag = 1;
    plotflag = 1;
    
elseif nargin < 5
    cmin = min(data(:));
    cmax = max(data(:));
    lineflag = 1;
    plotflag = 1;
    
elseif nargin < 7
    lineflag = 1;
    plotflag = 1;
    
elseif nargin < 8
    plotflag = 1;
end

if (dimnum == 3)

    axi_slice = flipdim(permute(squeeze(data(:,:,axi_ind,:)), [2, 1, 3:length(size(data))]), 2) ;
    cor_slice = flipdim(permute(squeeze(data(:,cor_ind, :,:)), [2, 1, 3:length(size(data))]), 1);
    sag_slice = flipdim(permute(squeeze(data(sag_ind,:,:,:)), [2, 1, 3:length(size(data))]), 1);

elseif (dimnum == 4)
    
    axi_slice = (data(:,:,axi_ind,:));
    cor_slice = flipdim(data(:, cor_ind, :,:), 3);
    sag_slice = flipdim(permute((data(sag_ind,:,:,:)), [2, 1, 3:length(size(data))]), 3);
end


if plotflag == 1

    % plot axial plane
    h1 = figure; imagesc(axi_slice); caxis([cmin, cmax]); axis equal; axis tight; axis off; colorbar;
    colormap(gray); 
    if lineflag == 1
        hold on; plot(1:size(axi_slice,2), cor_ind.*ones(size(axi_slice,2), 1), 'g-');
        plot(sag_ind.*ones(size(axi_slice,1), 1), 1:size(axi_slice,1), 'b-');
    end
    set(gca, 'FontSize', 20)


    % plot coronal plane
    h2 = figure; imagesc(cor_slice); caxis([cmin, cmax]); axis equal; axis tight; axis off; colorbar;
    colormap(gray); 
    if lineflag == 1
        hold on; plot(1:size(cor_slice,2), (ImageSlices-axi_ind+1).*ones(size(cor_slice,2), 1), 'r-');
        plot(sag_ind.*ones(size(cor_slice,1), 1), 1:ImageSlices, 'b-');
    end
    set(gca, 'FontSize', 20)

    % plot saggital plane
    h3 = figure; imagesc(sag_slice); caxis([cmin, cmax]); axis equal; axis tight; axis off; colorbar;
    colormap(gray); 
    if lineflag == 1
        hold on; plot(1:size(sag_slice,2), (ImageSlices-axi_ind+1).*ones(size(sag_slice,2), 1), 'r-');
        plot(cor_ind.*ones(size(sag_slice,1), 1), 1:ImageSlices, 'g-');
    end
    set(gca, 'FontSize', 20)
    
else
    
    h1 = mimage(axi_slice, cmin, cmax);
    h2 = mimage(cor_slice, cmin, cmax);
    h3 = mimage(sag_slice, cmin, cmax);
    
end