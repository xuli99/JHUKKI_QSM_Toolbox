function [SMVdata, background] = iRSHARPv1(phaseMap, wrap_raw, mask, Params, handles)
% ========================================================================================= 
% function [SMVdata, background] = iRSHARPv1(phaseMap, wrap_raw, mask, Params, handles)
% 
% Inputs:
%   phaseMap        - unwrapped phase map
%   wrap_raw        - wrapped phase map
%   mask            - Binary mask
%   Params          - Parameters of the acquired GRE, including:
%                     matrix size -Params.sizeVol
%                     resolution  -Params.voxSize
%                     FOV         -Params.fov
%                     echo time   -Params.TEs
%                     main field  -Params.B0
%                     gyro        -Params.gamma 
%   important parameters:
%   Params.SHARPradius      - Maximum radius
%   Params.thresh_tsvd      - Truncated SVD threshold
%   Params.iRSHARP_C        - The adjusting constant
% outputs:
%   SMVdata       - Local field
%   background    - Background field
%
% Originally Written by Jinsheng
% Created:  2017.08
% Modified: 2018.08
% When using the code, please cite
% J.Fang, et al, MRM 2018;
% J.Fang, et al, JMR 2017;
% Schweser et al. 2011, NIMG,
% Wu, et al., 2011, MRM
% -------------------------------------------------------------------------
% Acknowledgement
% Dr. Xu Li, Prof. Peter van Zijl, from JHU, USA;
% =========================================================================================== 
%
% modifed to fit the QSM_Toolbox format, X.L., 2019
% phaseMap input could be multi-echo 4D dataset, do it echo by echo
% updated 2019-07-09
% Updated 2021-06-27 X.L., cluster version

warning off all

% Start the clock
textWaitbar = ['Performing iRSHARP on ' num2str(length(Params.echoNums)) ' echoes'];
if ~isfield(Params, 'cluster')
    multiWaitbar(textWaitbar, 0, 'CanCancel', 'On' );
else
    disp(textWaitbar);
end

% Cut only what we want
GREPhaseData = phaseMap(:,:,:,Params.echoNums);
wrap_raw = wrap_raw(:,:,:,Params.echoNums);
Necho = size(GREPhaseData, 4);

Nx = Params.sizeVol(2);
Ny = Params.sizeVol(1);
Nz = Params.sizeVol(3);

dx = Params.voxSize(2);
dy = Params.voxSize(1);
dz = Params.voxSize(3);

radiusMAX = Params.SHARPradius;
radiusMIN = max([dx, dy, dz]);
radiusSTEP = radiusMIN;

radiusArray = radiusMIN:radiusSTEP:radiusMAX;
radiusNum = int16(length(radiusArray));

if radiusNum < 1
    disp('Error in SHARP radius definition.');
    return;
end

%% making SMV kernels
%% convolution kernel with maximum radius is used to generate deconvolution kernel: deconvk
SMVArray = cell(radiusNum, 1);      % SMV kernels

radius_xArray = zeros(radiusNum, 1, 'int16');
radius_yArray = zeros(radiusNum, 1, 'int16');
radius_zArray = zeros(radiusNum, 1, 'int16');

for i = 1:radiusNum
    
    radius = radiusArray(i);
    
    radius_x = ceil(radius./dx);       % in the unit of pixel
    radius_y = ceil(radius./dy);
    radius_z = ceil(radius./dz);
    
    Nx_kernel = 2*radius_x + 1;        % diameter, number of voxel
    Ny_kernel = 2*radius_y + 1;
    Nz_kernel = 2*radius_z + 1;
    
    x = linspace(-radius_x*dx, radius_x*dx, Nx_kernel);
    y = linspace(-radius_y*dy, radius_y*dy, Ny_kernel);
    z = linspace(-radius_z*dz, radius_z*dz, Nz_kernel);
    
    [X_Grid, Y_Grid, Z_Grid] = meshgrid(x, y, z);  % mesh in space
    R_Grid = sqrt(X_Grid.^2 + Y_Grid.^2 + Z_Grid.^2);
    
    rho = (R_Grid <= radius);
    rho = rho./sum(rho(:));                        % normalizing
   
    SMVArray{i} = rho;

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

%% SGK kernels depends on sigma, thus has to be calculated for each voxel

%% convolution with SMV (smallest) & subtraction
% Do the loop for all echoes!
dataConvRho = zeros(size(GREPhaseData));
SMVkernel = zeros(Ny, Nx, Nz);
SMVkernel(floor(Ny/2)-radius_yArray(1)+1:floor(Ny/2)+radius_yArray(1)+1, ...
    floor(Nx/2)-radius_xArray(1)+1:floor(Nx/2)+radius_xArray(1)+1, ...
    floor(Nz/2)-radius_zArray(1)+1:floor(Nz/2)+radius_zArray(1)+1) = SMVArray{1};
SMVkernel = fftn(fftshift(SMVkernel));

for selectedEcho = 1:Necho
    dataConvRho(:,:,:,selectedEcho) = real(ifftn(fftn(GREPhaseData(:,:,:,selectedEcho)).*SMVkernel));
end

% Count datapoints for the waitbar
totPoints = sum(tovec(mask > 0));

for selectedEcho = 1:Necho
    
    %% calculate spatial weighting for each echo
    wrapcount= abs((GREPhaseData(:,:,:,selectedEcho) - wrap_raw(:,:,:,selectedEcho))./(2*pi)).*mask;    
    SpaMat = iRSHARP_SpaWeight(GREPhaseData(:,:,:,selectedEcho),Params,wrapcount).*mask;  % Spatial weighting matrix
    
    SpaMat_max=max(SpaMat(:));
    SpaMat_mean=mean(SpaMat(mask>0));
    SpaMat_thresh=SpaMat_mean + (SpaMat_max-SpaMat_mean)/(SpaMat_max+SpaMat_mean);        % threshold  
        
    pointCounter = 0;
    
    % do image based convolution basd on spatial weighting SpaMat
    % 3D convolution in real space, adaptive inside the brainmask
    for kk = 1:Nz
        for jj = 1:Nx
            for ii = 1:Ny
                
                % update dataConvRho inside the mask
                if mask(ii,jj,kk) > 0
                    % determine the right kernel size
                    ir = radiusNum;
                    SpaWeight = SpaMat(ii,jj,kk);
                    
                    while ir > 0
                        
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
                        
                        if isempty(find(tempmask(SMVArray{ir}>0) < 10*eps, 1))    %  there is no zeros
                          
                            if (SpaWeight<=SpaMat_thresh)   
                                % update dataConvRho with the right SMV kernel
                                temp = GREPhaseData(tempii, tempjj, tempkk, selectedEcho).*SMVArray{ir};
                                dataConvRho(ii, jj, kk, selectedEcho) = sum(temp(:));
                            else
                                % calculate and use SGK kernel
                                sigma=1/(SpaWeight*Params.iRSHARP_C);
                                SGkernel = SMVArray{ir}.*Gaussiankernel3d(size(SMVArray{ir}), sigma);
                                SGkernel = SGkernel./(sum(SGkernel(:)));
                                temp = GREPhaseData(tempii, tempjj, tempkk, selectedEcho).*SGkernel;
                            end   

                             dataConvRho(ii, jj, kk, selectedEcho) = sum(temp(:)); 
                             break;
                        else
                            ir = ir - 1;
                        end
                        
                    end % while end

                    % Increase counter
                    pointCounter = pointCounter+1;
                    % Only show waitbar every 1000 iterations
                    if (mod(pointCounter, 1000) == 0)
                        if ~isfield(Params, 'cluster')
                            hasCanceled = multiWaitbar(textWaitbar, (pointCounter + (selectedEcho-1)*totPoints)/(totPoints*Necho));
                            HandleStopReconstruction;
                        end
                    end

                end
            end
        end
    end
    
end % echo end

% Repeat mask
data_interm = (GREPhaseData - dataConvRho).*repmat(mask, [1 1 1 size(GREPhaseData, 4)]);

clear temp tempmask tempii tempjj tempkk dataConvRho

%% deconvolution of (delta - rho) with masking
% Repload
SMVdata = zeros(size(GREPhaseData));

temp = zeros(Ny, Nx, Nz);
temp(floor(Ny/2)-radius_yArray(end)+1:floor(Ny/2)+radius_yArray(end)+1, ...
    floor(Nx/2)-radius_xArray(end)+1:floor(Nx/2)+radius_xArray(end)+1, ...
    floor(Nz/2)-radius_zArray(end)+1:floor(Nz/2)+radius_zArray(end)+1) = deconvk;

Hk = fftn(fftshift(temp));
mask_reg = abs(Hk) < Params.thresh_tsvd;   % changed

% For each echo
for selectedEcho = 1:size(GREPhaseData,4)
        
    temp = fftn(data_interm(:,:,:, selectedEcho))./Hk;    % tsvd
    temp(mask_reg) = 0;    
    SMVdata(:,:,:,selectedEcho) = real(ifftn(temp));

    %% Applying the mask
    SMVdata(:,:,:,selectedEcho) = SMVdata(:,:,:,selectedEcho).*mask;
    
end

% Calculate background
background = (GREPhaseData - SMVdata).*repmat(mask, [1 1 1 size(GREPhaseData, 4)]);

% Average everything
SMVdata = mean(SMVdata,4);
background = mean(background,4);

end

