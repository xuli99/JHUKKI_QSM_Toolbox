function [SMVdata, background] = SHARP_Adaptive_Multi2D(GREPhaseSE, mask, radiusMAX, thresh_reg, Params)
% [SMVdata, background] = SHARP_Adaptive_Multi2D(GREPhaseSE, mask, radiusMAX, thresh_reg, Params)
% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
%
% Ref: Schweser et al. 2011, NIMG, Wu, et al., 2011, MRM
%
% NOTE:
%   mask is the full brain mask without erosion, it is also the region of supoort (ROS)
%     where phase data is trustable
%     SHARP convolution is calculated in the ROS with adaptive radius
%   radiusMAX is the maximum radius in the unit of mm

% Modified to 2D version, 2017-04-03, 2D version does not use waitbar thus
% no need input handles
% Updated by Xu Li, 2019-07-03

warning off all

% Cut only what we want
GREPhaseData = GREPhaseSE(:,:,:,Params.echoNums);

Nx = Params.sizeVol(2);
Ny = Params.sizeVol(1);

dx = Params.voxSize(2);
dy = Params.voxSize(1);

radiusMIN = max([dx, dy]);
radiusSTEP = radiusMIN;

radiusArray = radiusMIN:radiusSTEP:radiusMAX;
radiusNum = int16(length(radiusArray));

if radiusNum < 1
    disp('Error in SHARP radius definition.');
    return;
end

%% making convolution kernels
%% convolution kernel with maximum radius, deconvolution kernel: deconvk
rhoArray = cell(radiusNum, 1);
radius_xArray = zeros(radiusNum, 1, 'int16');
radius_yArray = zeros(radiusNum, 1, 'int16');

for i = 1:radiusNum
    
    radius = radiusArray(i);
    
    radius_x = ceil(radius./dx);       % in the unit of pixel
    radius_y = ceil(radius./dy);
    
    Nx_kernel = 2*radius_x + 1;        % diameter, number of voxel
    Ny_kernel = 2*radius_y + 1;
    
    x = linspace(-radius_x*dx, radius_x*dx, Nx_kernel);
    y = linspace(-radius_y*dy, radius_y*dy, Ny_kernel);
    
    [X_Grid, Y_Grid] = meshgrid(x, y);  % mesh in space
    R_Grid = sqrt(X_Grid.^2 + Y_Grid.^2);
    
    rho = (R_Grid <= radius);
    rho = rho./sum(rho(:));                        % normalizing
    
    rhoArray{i} = rho;
    radius_xArray(i) = int16(radius_x);
    radius_yArray(i) = int16(radius_y);
    
    if i == radiusNum        
        % deconvolution kernel
        deconvk = rho;
        deconvk(radius_y+1, radius_x+1) = deconvk(radius_y+1, radius_x+1) - 1;
        deconvk = -deconvk;
    end
    
end

clear X_Grid Y_Grid R_Grid

%% convolution with rho & subtraction
% Do the loop for all echoes!
dataConvRho = zeros(size(GREPhaseData), class(GREPhaseData(1)));
for selectedEcho = 1:size(GREPhaseData,4)
    temp = zeros(Ny, Nx);
    temp(floor(Ny/2)-radius_yArray(1)+1:floor(Ny/2)+radius_yArray(1)+1, ...
        floor(Nx/2)-radius_xArray(1)+1:floor(Nx/2)+radius_xArray(1)+1) = rhoArray{1};
    dataConvRho(:,:,:,selectedEcho) = real(ifftn(fftn(GREPhaseData(:,:,:,selectedEcho)).*fftn(fftshift(temp))));

end

    for jj = 1:Nx
        for ii = 1:Ny
            % update dataConvRho inside the mask
            if mask(ii,jj) > 0
                
                % determine the right kernel size
                ir = radiusNum;
                while ir > 1
                    tempii = ii-radius_yArray(ir):ii+radius_yArray(ir);
                    tempjj = jj-radius_xArray(ir):jj+radius_xArray(ir);
                    
                    if (~isempty(find(tempii < 1, 1)) || ~isempty(find(tempii > Ny, 1)) || ...
                            ~isempty(find(tempjj < 1, 1)) || ~isempty(find(tempjj > Nx, 1)) )
                        % if tempii, tempjj or tempkk is invalide index
                        ir = ir - 1;
                        continue;
                    end
                    
                    tempmask = mask(tempii, tempjj);
                    if isempty(find(tempmask(rhoArray{ir}>0) < 10*eps, 1))    %  there is no zeros
                        % For each echo
                        for selectedEcho = 1:size(GREPhaseData,4)
                            % update dataConvRho with the right kernel
                            temp = GREPhaseData(tempii, tempjj, 1, selectedEcho).*rhoArray{ir};  % convolving with the selected kernel
                            dataConvRho(ii, jj, 1, selectedEcho) = sum(temp(:));
                        end
                        break;
                    else
                        ir = ir - 1;
                    end
                end

            end
        end
    end

% Repeat mask
data_interm = (GREPhaseData - dataConvRho).*repmat(mask, [1 1 1 size(GREPhaseData, 4)]);

clear temp tempmask tempii tempjj dataConvRho

%% deconvolution of (delta - rho) with masking
SMVdata = zeros(size(GREPhaseData), class(GREPhaseData(1)));

temp = zeros(Ny, Nx);
temp(floor(Ny/2)-radius_yArray(end)+1:floor(Ny/2)+radius_yArray(end)+1, ...
    floor(Nx/2)-radius_xArray(end)+1:floor(Nx/2)+radius_xArray(end)+1) = deconvk;


Hk = fftn(fftshift(temp));
mask_reg = abs(Hk) < thresh_reg;   % changed

% For each echo
for selectedEcho = 1:size(GREPhaseData,4)
        
    temp = fftn(data_interm(:,:,:, selectedEcho))./Hk;    % tsvd
    temp(mask_reg) = 0;
    SMVdata(:,:,:,selectedEcho) = real(ifftn(temp));
    
    SMVdata(:,:,:,selectedEcho) = temp(1:Ny, 1:Nx);
    
    %% Applying the mask
    SMVdata(:,:,:,selectedEcho) = SMVdata(:,:,:,selectedEcho).*mask;
    
end

% Calculate background
background = (GREPhaseData - SMVdata).*repmat(mask, [1 1 1 size(GREPhaseData, 4)]);

end
