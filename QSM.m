function varargout = QSM(varargin)
% QSM Toolbox
%      QSM Toolbox
%      Authors: Jiri van Bergen and Xu Li
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help QSM

% Last Modified by GUIDE v2.5 11-Oct-2023 12:22:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @QSM_OpeningFcn, ...
    'gui_OutputFcn',  @QSM_OutputFcn, ...
    'gui_LayoutFcn',  [], ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before QSM is made visible.
function QSM_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for QSM
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Add modules
FilePath = fileparts(mfilename('fullpath'));
addpath(FilePath);
addpath(genpath(fullfile(FilePath, 'QSM_Modules')));
addpath(genpath(fullfile(FilePath, 'QSM_NIFTI')));
addpath(genpath(fullfile(FilePath, 'QSM_Utility')));

% --- Outputs from this function are returned to the command line.
function varargout = QSM_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in ButtonSelectFile.
function ButtonSelectFile_Callback(hObject, eventdata, handles)
% Load constants
Constants;
% Open files
OpenFiles;

% --- Executes on button press in CheckSaveData.
function CheckSaveData_Callback(hObject, eventdata, handles)
% Will output be saved?
handles.Params.saveOutput = get(hObject,'Value');
% Save
guidata(hObject, handles);

% --- Executes on slider movement.
function ImageSlider_Callback(hObject, eventdata, handles)
% hObject    handle to ImageSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.slice2disp = round(get(hObject,'Value'));
ShowImage;

% --- Executes during object creation, after setting all properties.
function ImageSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in ButtonStart.
function ButtonStart_Callback(hObject, eventdata, handles)

%-----------------   Get all GUI input variables here first and save to Params
Params      = handles.Params;
% Get variables if there is not user interaction

% 1. Get variable for phase unwrapping
Params.UnwrappingMethod   = get(handles.VarUnwrappingMethod, 'Value');

% Updated by Xu Li, 2019-09-07, Added LapPhaseCorrection
% Add in corrections while using Path based unwrapping
% Such correction may affect quantifications, thus switched back
switch Params.UnwrappingMethodsDict{Params.UnwrappingMethod}
    case {'Path', 'NonlinearFit + Path'}
        Params.LapPhaseCorrection = 0;    % 1: add Correction / switched back
    case 'Laplacian'
        Params.LapPhaseCorrection = 0;
    otherwise
        error('Unknown unwrapping method.')        
end

Params.TemplateEcho = str2double(get(handles.VarTemplateEcho, 'String'));

% 2. Get variable for BrainMask
Params.SaveEcho = eval(get(handles.VarMaskEchoes, 'String'));
Params.FSLFolder = get(handles.VarFSL,'String');
Params.FSLThreshold = get(handles.VarFSLThres,'String');   % string send to bet command
Params.ErodeRadius  = str2double(get(handles.VarRadiusDisk,'String'));  

% 3. Get variable for background removal
Params.BgRemoval    = get(handles.VarBgRemoval,'Value');
Params.SHARPradius  = str2double(get(handles.VarSHARPradius,'String'));

% 4. Get variables for QSM
Params.QSMSolver    = get(handles.VarQSMSolver,'Value');
Params.B0           = str2double(get(handles.VarB0,'String'));

% if no click event
Params.phase2DprocFlag = get(handles.checkbox2Dproc,'Value');   
Params.EchoAvg = get(handles.checkbox_EchoAvg,'Value');         
Params.R2starFlag = get(handles.checkboxR2star,'Value');
Params.AutoRefFlag = get(handles.checkboxAutoRef,'Value');

% Update handles
handles.Params = Params;
guidata(hObject, handles);

% -----------------------------------------------
% LETS DO THIS THING!!
PerformUnwrapping;
CreateBrainMask;
RemoveBackground;
CalculateQSM;   % after done will enable MultiProcess


% --- Executes on key press with focus on VarB0 and none of its controls.
function VarB0_KeyPressFcn(hObject, eventdata, handles)
% Act on right key
if (strcmpi(eventdata.Key, '7'))
    % Find within the limits of TE
    handles.Params.echoNums = find(handles.Params.TEs >= handles.TELowerLimit_7T & handles.Params.TEs <= handles.TEUpperLimit_7T);
    set(handles.TextEchoInfo, 'String',  sprintf('Using echoes  %d - %d  (TE %0.3g - %0.3g ms)', min(handles.Params.echoNums) , max(handles.Params.echoNums), handles.Params.TEs(min(handles.Params.echoNums))*1000, handles.Params.TEs(max(handles.Params.echoNums))*1000));
elseif (strcmpi(eventdata.Key, '3'))
    % Find within the limits of TE
    handles.Params.echoNums = find(handles.Params.TEs >= handles.TELowerLimit_3T & handles.Params.TEs <= handles.TEUpperLimit_3T);
    set(handles.TextEchoInfo, 'String',  sprintf('Using echoes  %d - %d  (TE %0.3g - %0.3g ms)', min(handles.Params.echoNums) , max(handles.Params.echoNums), handles.Params.TEs(min(handles.Params.echoNums))*1000, handles.Params.TEs(max(handles.Params.echoNums))*1000));
elseif (strcmpi(eventdata.Key, '1.5'))
    % Find within the limits of TE
    handles.Params.echoNums = find(handles.Params.TEs >= handles.TELowerLimit_1p5T & handles.Params.TEs <= handles.TEUpperLimit_1p5T);
    set(handles.TextEchoInfo, 'String',  sprintf('Using echoes  %d - %d  (TE %0.3g - %0.3g ms)', min(handles.Params.echoNums) , max(handles.Params.echoNums), handles.Params.TEs(min(handles.Params.echoNums))*1000, handles.Params.TEs(max(handles.Params.echoNums))*1000));
elseif (strcmpi(eventdata.Key, '11.7'))
    % Find within the limits of TE
    handles.Params.echoNums = find(handles.Params.TEs >= handles.TELowerLimit_11p7T & handles.Params.TEs <= handles.TEUpperLimit_11p7T);
    set(handles.TextEchoInfo, 'String',  sprintf('Using echoes  %d - %d  (TE %0.3g - %0.3g ms)', min(handles.Params.echoNums) , max(handles.Params.echoNums), handles.Params.TEs(min(handles.Params.echoNums))*1000, handles.Params.TEs(max(handles.Params.echoNums))*1000));
end
% Fill masking echo
handles.Params.SaveEcho = handles.Params.echoNums(1);
set(handles.VarMaskEchoes,'String',sprintf('[%d]',handles.Params.echoNums(1)));
% Save
guidata(hObject, handles);

% for displaying image slices after all the calculations ...
% --- Executes on button press in ButtonShowRawPhase.
function ButtonShowRawPhase_Callback(hObject, eventdata, handles)
% Set raw phase
LoadImage(hObject, handles, handles.GREPhaseRaw, 'Phase (rad)');

% --- Executes on button press in ButtonShowPhase.
function ButtonShowPhase_Callback(hObject, eventdata, handles)
% Set phase
LoadImage(hObject, handles, handles.GREPhase, 'Phase (rad)');

% --- Executes on button press in ButtonShowMask.
function ButtonShowMask_Callback(hObject, eventdata, handles)
% Set Mask - TIMES TEN TO STOP IMSHOW FROM COMPLAINING
LoadImage(hObject, handles, handles.GREPhase(:,:,:,handles.Params.echoNums(1)).*handles.maskErode, 'Phase (rad)');

% --- Executes on button press in ButtonShowBg.
function ButtonShowBg_Callback(hObject, eventdata, handles)
% Set freqMap
LoadImage(hObject, handles, handles.freqMap, 'Frequency (Hz)');

% --- Executes on button press in ButtonShowQSM.
function ButtonShowQSM_Callback(hObject, eventdata, handles)
% Set QSM
LoadImage(hObject, handles, handles.chi_res, 'Susceptibility (ppm)');

% --- Executes on button press in ButtonAdjustContrast.
function ButtonAdjustContrast_Callback(hObject, eventdata, handles)
% Start that contrast thing
% set(handles.ImageContainer,'CLim', handles.ImageCLim);
imcontrast(handles.ImageContainer);

% --- Executes on button press in ButtonAddDataset.
function ButtonAddDataset_Callback(hObject, eventdata, handles)
AddDatasets;

% --- Executes on button press in ButtonStartDatasets.
function ButtonStartDatasets_Callback(hObject, eventdata, handles)
StartMultiProcess;

% --- Executes on button press in ButtonEditEchoes.
function ButtonEditEchoes_Callback(hObject, eventdata, handles)
% Popup for asking which echoes
if handles.Params.nEchoes > 1
    newEchoes = inputdlg('What echoes to use?', 'What echoes to use? (start:step:end)', [1,30], { sprintf('%d:%d:%d', min(handles.Params.echoNums) , handles.Params.echoStep, max(handles.Params.echoNums)) });
    if ~isempty(newEchoes)
        newEchoes = eval(newEchoes{:});    
        % Save
        handles.Params.echoNums = newEchoes;
        if length(newEchoes) > 1
            handles.Params.echoStep = newEchoes(2) - newEchoes(1);  
        else
            handles.Params.echoStep = 0;
        end
    end
else
    newEchoes = inputdlg('What echoes to use?', 'What echoes to use? (start:end)', [1,30], { sprintf('%d:%d', min(handles.Params.echoNums), max(handles.Params.echoNums)) });
    newEchoes = eval(newEchoes{:});    
    % Save
    handles.Params.echoNums = newEchoes;
end

set(handles.TextEchoInfo, 'String',  sprintf('Using echoes  %d - %d  (TE %0.3g - %0.3g ms)', min(handles.Params.echoNums) , max(handles.Params.echoNums), handles.Params.TEs(min(handles.Params.echoNums))*1000 , handles.Params.TEs(max(handles.Params.echoNums))*1000));
guidata(hObject, handles);

% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jFrame = get(hObject, 'JavaFrame');
try
  jProx = jFrame.fFigureClient.getWindow;
catch
  % jProx = jFrame.fHG1Client.getWindow;  % [EDITED] Fallback % obsolete  
end
% Now try the setting
try
    jProx.setMinimumSize(java.awt.Dimension(1250, 900));
catch
    % Nothing
end

% --- Executes on selection change in VarUnwrappingMethod.
function VarUnwrappingMethod_Callback(hObject, eventdata, handles)
% hObject    handle to VarUnwrappingMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns VarUnwrappingMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from VarUnwrappingMethod
handles.Params.UnwrappingMethod   = get(hObject,'Value');   % 1: Lap 2: Nonlinear + path 3: Path
% Save
guidata(hObject, handles);

contents = cellstr(get(hObject,'String'));
% set method combinations
% if Path --> show TemplateEcho option
if matches(contents{get(hObject,'Value')}, 'Path')
    set([handles.Tag_TemplateEcho_text], 'Visible', 'on')
    set([handles.VarTemplateEcho], 'Visible', 'on')
else
    set([handles.Tag_TemplateEcho_text], 'Visible', 'off')
    set([handles.VarTemplateEcho], 'Visible', 'off')
end

% % --- Executes during object creation, after setting all properties.
% function VarUnwrappingMethod_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to VarUnwrappingMethod (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: popupmenu controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end


% --- Executes on button press in ButtonLoadDataList.
function ButtonLoadDataList_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonLoadDataList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LoadDataList;

% --- Executes on button press in checkboxR2star.
function checkboxR2star_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxR2star (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkboxR2star

% if SFCR/nSFCR show option of AutoRef. 
QSMSolverDict   = cellstr(get(handles.VarQSMSolver, 'String'));
QSMSolver       = get(handles.VarQSMSolver,'Value');
R2starFlag      = get(hObject, 'Value');

if R2starFlag == 1 && contains(QSMSolverDict(QSMSolver), 'SFCR')
    set([handles.checkboxAutoRef], 'Visible', 'On')
else
    set([handles.checkboxAutoRef], 'Visible', 'Off')
    set([handles.checkboxAutoRef], 'Value', 0)
end

% Save
guidata(hObject, handles);

% --- Executes on selection change in VarQSMSolver.
function VarQSMSolver_Callback(hObject, eventdata, handles)
% hObject    handle to VarQSMSolver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns VarQSMSolver contents as cell array
%        contents{get(hObject,'Value')} returns selected item from VarQSMSolver
contents = cellstr(get(hObject,'String'));
% set method combinations
% if ~SFCR --> no AutoRef
% if SFCR + R2* --> AutoRef available
if ~contains(contents{get(hObject,'Value')}, 'SFCR')
    set([handles.checkboxAutoRef], 'Visible', 'Off')
    set([handles.checkboxAutoRef], 'Value', 0)
else
    R2starFlag = get(handles.checkboxR2star, 'Value');
    if R2starFlag == 1
        set([handles.checkboxAutoRef], 'Visible', 'On')
    else
        set([handles.checkboxAutoRef], 'Value', 0)
    end
end

function VarTemplateEcho_Callback(hObject, eventdata, handles)
% hObject    handle to VarTemplateEcho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VarTemplateEcho as text
%        str2double(get(hObject,'String')) returns contents of VarTemplateEcho as a double

% Check Values
TemplateEcho = str2double(get(hObject,'String'));
if TemplateEcho > handles.Params.nEchoes || TemplateEcho < 0
    errordlg('Invalide Template Echo.', 'Error in Entered Values');
end
