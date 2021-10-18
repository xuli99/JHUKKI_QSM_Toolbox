function f = showInfoWindow(msg,wtitle)

%
% SHOWINFOWINDOW creates a small dialog window which displays a user-defined message. It
% is quite like MATLAB's msgbox but without an OK button. Although msgbox could be used
% for the simple display of a message that informs the user about a process being
% performed (and is not expected to last long, else a waitbar/timebar could be used but at
% a cost of running time), the presence of the OK button that closes the window on press
% could be annoying. SHOWINFOWINDOW resolves this problem by removing the button.
% Additionally, if you wish to write a script or create a GUI where you do not wish the
% user to accidenatlly hit the close button in the upper right corner of the window, you
% can remove the comment on line 82 of the code concerning the 'CloseRequestFcn' property
% of the dialog. However, if you do that you must assign a handle to showinfowindow so you
% can change its 'CloseRequestFcn' property to 'closereq' when your process is done or
% else you might end up with a non-closing window! If that happens, one way to resolve
% this could be to use the findobj function:
%
% non_closing_window_handle = findobj('CloseRequestFcn','');
% set(non_closing_window_handle,'CloseRequestFcn','closereq')
% close(non_closing_window_handle)
%
% Syntax:
%
% showinfowindow
% h = showinfowindow
% showinfowindow(msg)
% h = showinfowindow(msg)
% showinfowindow(msg,wtitle)
% h = showinfowindow(msg,wtitle)
%
% showinfowindow and h = showinfowindow without any input arguments create an example.
% 
% showinfowindow(msg) and h = showinfowindow(msg) display the user-defined message msg
% which should be a character string or a cell array of strings.
%
% showinfowindow(msg,wtitle) and h = showinfowindow(msg,wtitle) display the user-defined 
% message msg which should be a character string or a cell array of strings. This time the
% window has the title wtitle
%
% Although this script was created for something trivial, I found it useful especially in
% my GUIs where I wanted to briefly display some information in a very simple way and
% especially when I wanted to avoid someone causing such windows to close when running a
% process within a GUI, probably resulting in errors.
%
% Examples
%
% h1 = showinfowindow('Running - Please wait...');
% 
% h2 = showinfowindow({'This is a long message for a long process.',...
%                      'It is displayed in two lines.'},'Long message');
%
%==========================================================
%
% Author        : Panagiotis Moulos (pmoulos@eie.gr)
% First created : May 7, 2007
% Last modified : June 11, 2007
%

% Check input arguments
if nargin<1
    msg='Example of showinfowindow.m';
    wtitle='Example';
elseif nargin<2
    wtitle='Info';
end
if ~iscell(msg)
    msg={msg};
end

% Set default sizes
DefFigPos=get(0,'DefaultFigurePosition');
MsgOff=7;
FigWidth=125;
FigHeight=50;
DefFigPos(3:4)=[FigWidth FigHeight];
MsgTxtWidth=FigWidth-2*MsgOff;
MsgTxtHeight=FigHeight-2*MsgOff;

% Initialize dialog window
f=dialog('Name',wtitle,'Units','points','WindowStyle','normal','Toolbar','none',...
         'DockControls','off','MenuBar','none','Resize','off','ToolBar','none',...
         'NumberTitle','off');%,'CloseRequestFcn','');

% Initialize message
msgPos=[MsgOff MsgOff MsgTxtWidth MsgTxtHeight];
msgH=uicontrol(f,'Style','text','Units','points','Position',msgPos,'String',' ',...
               'Tag','MessageBox','HorizontalAlignment','left','FontSize',8);
[WrapString,NewMsgTxtPos]=textwrap(msgH,msg,75);
set(msgH,'String',WrapString)
delete(msgH);

% Fix final message positions
MsgTxtWidth=max(MsgTxtWidth,NewMsgTxtPos(3));
MsgTxtHeight=min(MsgTxtHeight,NewMsgTxtPos(4));
MsgTxtXOffset=MsgOff;
MsgTxtYOffset=MsgOff;
FigWidth=MsgTxtWidth+2*MsgOff;
FigHeight=MsgTxtYOffset+MsgTxtHeight+MsgOff;

DefFigPos(3:4)=[FigWidth FigHeight];
set(f,'Position',DefFigPos);

% Create the message
AxesHandle=axes('Parent',f,'Position',[0 0 1 1],'Visible','off');
txtPos=[MsgTxtXOffset MsgTxtYOffset 0];
text('Parent',AxesHandle,'Units','points','HorizontalAlignment','left',...
     'VerticalAlignment','bottom','Position',txtPos,'String',WrapString,...
     'FontSize',8,'Tag','MessageBox');
 
% Move the window to the center of the screen
set(f,'Units','pixels')
screensize=get(0,'screensize');                       
winsize=get(f,'Position');
winwidth=winsize(3);
winheight=winsize(4);
screenwidth=screensize(3);                           
screenheight=screensize(4);                           
winpos=[0.5*(screenwidth-winwidth),0.5*(screenheight-winheight),winwidth,winheight];                          
set(f,'Position',winpos);

% Give priority to displaying this message
drawnow
