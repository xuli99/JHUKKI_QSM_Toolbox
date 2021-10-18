function [SMVdata, background] = SHARP_Adaptive_Multi(GREPhaseSE, mask, radiusMAX, thresh_reg, Params, handles)
% [SMVdata, background] = SHARP_Adaptive(data, mask, radiusMAX,thresh_reg,Params)
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
%   Params.NSR: estimated NSR used in the Wiener filter, may need to be tested out
%   The smaller the NSR, the bigger SNR, the less smoothing effect

% Updated by Xu Li, 2014-02-28
% Updated by Jiri, 2014-02-20 - Made it MULTI, so all echoes at once
% Updated by Jiri, 2014-03-01 - Added handles to allow for cancel button on the waitbar
% Updated by Xu Li, 2019-07-09
% Updated 2019-09-07: Added LapPhaseCorrection, X.L.
% Updated 2021-06-28 X.L., cluster version

warning off all

% Start the clock
textWaitbar = ['Performing V-SHARP on ' num2str(length(Params.echoNums)) ' echoes'];
if ~isfield(Params, 'cluster')
    multiWaitbar(textWaitbar, 0, 'CanCancel', 'On' );
else
    disp(textWaitbar);
end

% Cut only what we want
GREPhaseData = GREPhaseSE(:,:,:,Params.echoNums);
GRETEs = Params.TEs(Params.echoNums);

% Huge data?
if isa(GREPhaseData(1), 'single')
    dataTypeFlag = 1;
else
    dataTypeFlag = 2;
end

Nx = Params.sizeVol(2);
Ny = Params.sizeVol(1);
Nz = Params.sizeVol(3);

dx = Params.voxSize(2);
dy = Params.voxSize(1);
dz = Params.voxSize(3);

radiusMIN = max([dx, dy, dz]);
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
radius_zArray = zeros(radiusNum, 1, 'int16');

for i = 1:radiusNum
    
    radius = radiusArray(i);
    
    radius_x = ceil(radius./dx);       % in the unit of pixel
    radius_y = ceil(radius./dy);
    radius_z = ceil(radius./dz);
    
    Nx_kernel = 2*radius_x + 1;        
    Ny_kernel = 2*radius_y + 1;
    Nz_kernel = 2*radius_z + 1;
    
    x = linspace(-radius_x*dx, radius_x*dx, Nx_kernel);
    y = linspace(-radius_y*dy, radius_y*dy, Ny_kernel);
    z = linspace(-radius_z*dz, radius_z*dz, Nz_kernel);
    
    [X_Grid, Y_Grid, Z_Grid] = meshgrid(x, y, z);  % mesh in space
    R_Grid = sqrt(X_Grid.^2 + Y_Grid.^2 + Z_Grid.^2);
    
    rho = (R_Grid <= radius);
    rho = rho./sum(rho(:));                        % normalizing
    
    rhoArray{i} = rho;
    radius_xArray(i) = int16(radius_x);
    radius_yArray(i) = int16(radius_y);
    radius_zArray(i) = int16(radius_z);
    
    if i == radiusNum
        % deconvolution kernel
        deconvk = rho;
        deconvk(radius_y+1, radius_x+1, radius_z+1) = deconvk(radius_y+1, radius_x+1, radius_z+1) - 1;
        deconvk = -deconvk;
    end
    
end

clear X_Grid Y_Grid Z_Grid R_Grid

%% convolution with rho & subtraction
% Do the loop for all echoes!
dataConvRho = zeros(size(GREPhaseData));
for selectedEcho = 1:size(GREPhaseData,4)
    temp = zeros(Ny, Nx, Nz);
    temp(floor(Ny/2)-radius_yArray(1)+1:floor(Ny/2)+radius_yArray(1)+1, ...
        floor(Nx/2)-radius_xArray(1)+1:floor(Nx/2)+radius_xArray(1)+1, ...
        floor(Nz/2)-radius_zArray(1)+1:floor(Nz/2)+radius_zArray(1)+1) = rhoArray{1};
    dataConvRho(:,:,:,selectedEcho) = real(ifftn(fftn(GREPhaseData(:,:,:,selectedEcho)).*fftn(fftshift(temp))));
end

% Count datapoints for the waitbar
totPoints = sum(tovec(mask > 0));
pointCounter = 0;

% 3D convolution in real space, adaptive inside the brainmask
for kk = 1:Nz
    for jj = 1:Nx
        for ii = 1:Ny
            % update dataConvRho inside the mask
            if mask(ii,jj,kk) > 0
                
                % determine the right kernel size
                ir = radiusNum;
                while ir > 1
                    tempii = ii-radius_yArray(ir):ii+radius_yArray(ir);
                    tempjj = jj-radius_xArray(ir):jj+radius_xArray(ir);
                    tempkk = kk-radius_zArray(ir):kk+radius_zArray(ir);
                    
                    if (~isempty(find(tempii < 1, 1)) || ~isempty(find(tempii > Ny, 1)) || ...
                            ~isempty(find(tempjj < 1, 1)) || ~isempty(find(tempjj > Nx, 1)) || ...
                            ~isempty(find(tempkk < 1, 1)) || ~isempty(find(tempkk > Nz, 1)) )
                        % if tempii, tempjj or tempkk is invalide index
                        ir = ir - 1;
                        continue;
                    end
                    
                    tempmask = mask(tempii, tempjj, tempkk);
                    if isempty(find(tempmask(rhoArray{ir}>0) < 10*eps, 1))    %  there is no zeros
                        % For each echo
                        for selectedEcho = 1:size(GREPhaseData,4)
                            % update dataConvRho with the right kernel
                            temp = GREPhaseData(tempii, tempjj, tempkk, selectedEcho).*rhoArray{ir};  % convolving with the selected kernel
                            dataConvRho(ii, jj, kk, selectedEcho) = sum(temp(:));
                        end
                        break;
                    else
                        ir = ir - 1;
                    end
                end
                % Increase counter
                pointCounter = pointCounter+1;
                % Only show waitbar every 1000 iterations
                if(mod(pointCounter, 1000) == 0)
                    if ~isfield(Params, 'cluster')
                        hasCanceled = multiWaitbar(textWaitbar, (pointCounter/totPoints));
                        HandleStopReconstruction;
                    else
                        disp([num2str(100*(pointCounter/totPoints)), '% Done.'])
                    end
                end
            end
        end
    end
end

% Repeat mask
data_interm = (GREPhaseData - dataConvRho).*repmat(mask, [1 1 1 size(GREPhaseData, 4)]);

clear temp tempmask tempii tempjj tempkk dataConvRho

%% deconvolution of (delta - rho) with masking
if dataTypeFlag == 1
    deconvk = single(deconvk);
end

% Repload
SMVdata = zeros(size(GREPhaseData));

temp = zeros(Ny, Nx, Nz);
temp(floor(Ny/2)-radius_yArray(end)+1:floor(Ny/2)+radius_yArray(end)+1, ...
    floor(Nx/2)-radius_xArray(end)+1:floor(Nx/2)+radius_xArray(end)+1, ...
    floor(Nz/2)-radius_zArray(end)+1:floor(Nz/2)+radius_zArray(end)+1) = deconvk;

Hk = fftn(fftshift(temp));
mask_reg = abs(Hk) < thresh_reg;   

% For each echo
for selectedEcho = 1:size(GREPhaseData,4)
        
    temp = fftn(data_interm(:,:,:, selectedEcho))./Hk;    % tsvd
    temp(mask_reg) = 0;    
    SMVdata(:,:,:,selectedEcho) = real(ifftn(temp));

    %% Applying the mask
    SMVdata(:,:,:,selectedEcho) = SMVdata(:,:,:,selectedEcho).*mask;
    
    % LapPhaseCorrection, 2019-09-07, xl
    if Params.LapPhaseCorrection == 1
        temp = angle(exp(1i*SMVdata(:,:,:,selectedEcho).*(2*pi*GRETEs(selectedEcho))));
        temp = phase_unwrap_laplacian(temp, Params, 0, 2);     % no Ref
        SMVdata(:,:,:,selectedEcho) = temp./(2*pi*GRETEs(selectedEcho)).*mask;
    end    
    
end

% Calculate background
background = (GREPhaseData - SMVdata).*repmat(mask, [1 1 1 size(GREPhaseData, 4)]);

% Average everything
SMVdata = mean(SMVdata,4);
background = mean(background,4);
end
