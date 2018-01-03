function varargout = DTIvol(varargin)
% DTIVOL M-file for DTIvol.fig
%      DTIVOL, by itself, creates a new DTIVOL or raises the existing
%      singleton*.
%
%      H = DTIVOL returns the handle to a new DTIVOL or the handle to
%      the existing singleton*.
%
%      DTIVOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DTIVOL.M with the given input arguments.
%
%      DTIVOL('Property','Value',...) creates a new DTIVOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DTIvol_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DTIvol_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DTIvol

% Last Modified by GUIDE v2.5 24-Apr-2004 23:35:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DTIvol_OpeningFcn, ...
                   'gui_OutputFcn',  @DTIvol_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before DTIvol is made visible.
function DTIvol_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DTIvol (see VARARGIN)

% Choose default command line output for DTIvol
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DTIvol wait for user response (see UIRESUME)
% uiwait(handles.figure1);

assignin('base','doFAmap',0);
assignin('base','dotADC',0);
assignin('base','dolambda1',0);
assignin('base','dolambda2',0);
assignin('base','dolambda3',0);
assignin('base','docm',0);
assignin('base','noise',60);
extension='IMA';
assignin('base','extension',extension);

% --- Outputs from this function are returned to the command line.
function varargout = DTIvol_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in tagFA.
function tagFA_Callback(hObject, eventdata, handles)
% hObject    handle to tagFA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tagFA
doFAmap=get(hObject,'Value');
assignin('base','doFAmap',doFAmap);

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over tagFA.
function tagFA_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to tagFA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton1.
function pushbutton1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in tagtADC.
function tagtADC_Callback(hObject, eventdata, handles)
% hObject    handle to tagtADC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tagtADC
dotADC=get(hObject,'Value');
assignin('base','dotADC',dotADC);
% --- Executes on button press in tagColormap.
function tagColormap_Callback(hObject, eventdata, handles)
% hObject    handle to tagColormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tagColormap
docm=get(hObject,'Value');
assignin('base','docm',docm);

% --- Executes on button press in tage1.
function tage1_Callback(hObject, eventdata, handles)
% hObject    handle to tage1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tage1
dolambda1=get(hObject,'Value');
assignin('base','dolambda1',dolambda1);
% --- Executes on button press in tage2.
function tage2_Callback(hObject, eventdata, handles)
% hObject    handle to tage2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tage2
dolambda2=get(hObject,'Value');
assignin('base','dolambda2',dolambda2)
% --- Executes on button press in tage3.
function tage3_Callback(hObject, eventdata, handles)
% hObject    handle to tage3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tage3
%dolambda3='y';
dolambda3=get(hObject,'Value');
assignin('base','dolambda3',dolambda3);


% --- Executes during object creation, after setting all properties.
function tagnoise_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tagnoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function tagnoise_Callback(hObject, eventdata, handles)
% hObject    handle to tagnoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tagnoise as text
%        str2double(get(hObject,'String')) returns contents of tagnoise as a double
noise=str2double(get(hObject,'String'));


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


DTIguicode



% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over tagnoise.
function tagnoise_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to tagnoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%noise=str2double(get(hObject,'String'))


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over checkbox7.
function checkbox7_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over noise0.
function noise0_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to noise0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in noise0.
function noise0_Callback(hObject, eventdata, handles)
% hObject    handle to noise0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noise0
noisestate=get(hObject,'Value')
if (noisestate==1)
    noise=0
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over tagtADC.
function tagtADC_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to tagtADC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


dotADC=get(hObject,'Value')
assignin('base','dotADC',dotADC)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over frame1.
function frame1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to frame1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
noise=str2double(get(hObject,'String'));
assignin('base','noise',noise);


% --- Executes during object creation, after setting all properties.
function extension_CreateFcn(hObject, eventdata, handles)
% hObject    handle to extension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function extension_Callback(hObject, eventdata, handles)
% hObject    handle to extension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of extension as text
%        str2double(get(hObject,'String')) returns contents of extension as a double

extension=get(hObject,'String');
assignin('base','extension',extension);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over extension.
function extension_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to extension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


