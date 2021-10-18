%% Author: Jiri van Bergen and Xu Li
% Affiliation: Radiology @ JHU - Kirby Center
% Contact via xuli@mri.jhu.edu

% Updated by Xu Li, 2019-07-09

%% Show that image!

% Enable button
set(handles.ButtonAdjustContrast, 'Enable', 'On');

% Check dimensions
if size(handles.CurrentImage,5) > 1
    imageSlice = handles.CurrentImage(:,:,handles.slice2disp, min(handles.Params.echoNums), 1); % show the first dynamic
elseif size(handles.CurrentImage,4) > 1
    imageSlice = handles.CurrentImage(:,:,handles.slice2disp, min(handles.Params.echoNums));
else
    imageSlice = handles.CurrentImage(:,:,handles.slice2disp);
end

% New image? Then new IMSHOW
if(exist('isNewImage','var') && isNewImage)
    handles.theImage = imshow(imageSlice, [min(imageSlice(:))  max(imageSlice(:))], 'Parent', handles.ImageContainer);
    isNewImage = 0;
else
    % Overwrite the image data
    set(handles.theImage,'CData', imageSlice, 'parent', handles.ImageContainer);
end

colorbar;
colormap(gray)
set(handles.ImageContainer, 'XTick', [], 'YTick', [])

% Save changes
guidata(hObject, handles);