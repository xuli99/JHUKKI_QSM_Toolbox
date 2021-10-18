function handles = LoadImage(hObject, handles, imageData, imageLegend)
    %% Author: Jiri van Bergen
    % Affiliation: Radiology @ JHU - Kirby Center
    % Contact via xuli@mri.jhu.edu
    % updated 2019-06-20, X.L.
    
    handles.CurrentImage    = imageData;
    guidata(hObject, handles);
    
    % Update legend
    set(handles.TextLegend, 'String', imageLegend)
    
    % Display
    isNewImage = 1;     % flag used in ShowImage
    ShowImage;
end