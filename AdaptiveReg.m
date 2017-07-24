function varargout = AdaptiveReg(varargin)
% ADAPTIVEREG MATLAB code for AdaptiveReg.fig
%      ADAPTIVEREG, by itself, creates a new ADAPTIVEREG or raises the existing
%      singleton*.
%
%      H = ADAPTIVEREG returns the handle to a new ADAPTIVEREG or the handle to
%      the existing singleton*.
%
%      ADAPTIVEREG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ADAPTIVEREG.M with the given input arguments.
%
%      ADAPTIVEREG('Property','Value',...) creates a new ADAPTIVEREG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AdaptiveReg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AdaptiveReg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AdaptiveReg

% Last Modified by GUIDE v2.5 18-Jul-2017 17:56:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AdaptiveReg_OpeningFcn, ...
                   'gui_OutputFcn',  @AdaptiveReg_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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


% UIWAIT makes AdaptiveReg wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AdaptiveReg_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes just before AdaptiveReg is made visible.
function AdaptiveReg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AdaptiveReg (see VARARGIN)

% Choose default command line output for AdaptiveReg
handles.output = hObject;

if exist('Adaptive.mat','file')
    load('Adaptive','SelectionCriterion','SelectionTol','thresh',...
        'MaxOutlierIter','MinIter','AvgNVal','TargetErrorNorm',...
        'MaxOrder','StoppingCriterion');
    set(handles.edit_BoundariesTolerance,'String',sprintf('%d',SelectionTol));
    set(handles.edit_ZeroSupuriousMode,'String',sprintf('%d',thresh));
    set(handles.edit_MaxOutlierIt,'String',sprintf('%d',MaxOutlierIter));
    set(handles.MinIt,'String',sprintf('%d',MinIter));
    set(handles.edit_AvgNval,'String',sprintf('%d',AvgNVal));
    set(handles.edit_ToleranceForConverg,'String',sprintf('%d',TargetErrorNorm));
    set(handles.edit_MaxOrderDandB,'String',sprintf('%d',MaxOrder));   
    set(handles.popupmenu_SelectionCrit,'Value',SelectionCriterion);
    set(handles.popupmenu_StoppingCrit,'Value',StoppingCriterion);

end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

SelectionCriterion = get(handles.popupmenu_SelectionCrit,'Value');
SelectionTol = str2double(get(handles.edit_BoundariesTolerance,'String'));
thresh = str2double(get(handles.edit_ZeroSupuriousMode,'String'));
MaxOutlierIter = str2double(get(handles.edit_MaxOutlierIt,'String'));
MinIter = str2double(get(handles.MinIt,'String'));
AvgNVal = str2double(get(handles.edit_AvgNval,'String'));
TargetErrorNorm = str2double(get(handles.edit_ToleranceForConverg,'String'));
MaxOrder = str2double(get(handles.edit_MaxOrderDandB,'String'));
StoppingCriterion = get(handles.popupmenu_StoppingCrit,'Value');

if MinIter < AvgNVal
    errordlg('Minimum number of iterations must be larger than the number of iterations to compute the average energy','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% Saving the workspace to the local |mat| file. If requested by the
% user, the GUI information is also _appended_ to the file specified by the
% user.
load('HeatStructDef','DirName','FileName');
% if requested by the user, saves to file
if ~isempty(DirName)
    save(fullfile(DirName,FileName),'-append',...
        'SelectionCriterion','SelectionTol','thresh',...
        'MaxOutlierIter','MinIter','AvgNVal','TargetErrorNorm',...
        'MaxOrder','StoppingCriterion');
end

save('Adaptive','SelectionCriterion','SelectionTol','thresh',...
        'MaxOutlierIter','MinIter','AvgNVal','TargetErrorNorm',...
        'MaxOrder','StoppingCriterion');
close(handles.figure1);

CheckHeatReg;


% --- Executes on button press in previous.
function previous_Callback(hObject, eventdata, handles)
% hObject    handle to previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.figure1);
HeatRegBC2;




% --- Executes on selection change in popupmenu_SelectionCrit.
function popupmenu_SelectionCrit_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_SelectionCrit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_SelectionCrit contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_SelectionCrit


% --- Executes during object creation, after setting all properties.
function popupmenu_SelectionCrit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_SelectionCrit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in popupmenu_StoppingCrit.
function popupmenu_StoppingCrit_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_StoppingCrit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_StoppingCrit contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_StoppingCrit


% --- Executes during object creation, after setting all properties.
function popupmenu_StoppingCrit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_StoppingCrit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_BoundariesTolerance_Callback(hObject, eventdata, handles)
% hObject    handle to edit_BoundariesTolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_BoundariesTolerance as text
%        str2double(get(hObject,'String')) returns contents of edit_BoundariesTolerance as a double
a = str2double(get(hObject,'String'));

if isnan(a) || ~isreal(a) || a > 1 || a <= 0
    set(hObject,'String','0.99');    
    errordlg('You must enter a value between 0 and 1','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_BoundariesTolerance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_BoundariesTolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_ZeroSupuriousMode_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ZeroSupuriousMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ZeroSupuriousMode as text
%        str2double(get(hObject,'String')) returns contents of edit_ZeroSupuriousMode as a double
a = str2double(get(hObject,'String'));

if isnan(a) || ~isreal(a) || a < 0 || a > 0.1
    set(hObject,'String','1.e-12');
    errordlg('You must enter a real value between 0.0 and 0.1','Invalid input','modal');
    uicontrol(hObject);
    return;
end


% --- Executes during object creation, after setting all properties.
function edit_ZeroSupuriousMode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ZeroSupuriousMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_MaxOutlierIt_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MaxOutlierIt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MaxOutlierIt as text
%        str2double(get(hObject,'String')) returns contents of edit_MaxOutlierIt as a double
a = str2double(get(hObject,'String'));

if isnan(a) || ~isreal(a) || logical(abs(round(a)-a)<eps)==0 || isequal(a,0) || a < 0
    set(hObject,'String','6');
    errordlg('You must enter an integer positive value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_MaxOutlierIt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MaxOutlierIt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MinIt_Callback(hObject, eventdata, handles)
% hObject    handle to MinIt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinIt as text
%        str2double(get(hObject,'String')) returns contents of MinIt as a double

a = str2double(get(hObject,'String'));

if isnan(a) || ~isreal(a) || logical(abs(round(a)-a)<eps)==0 || a < 1
    set(hObject,'String','5');
    errordlg('You must enter an integer positive value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function MinIt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinIt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_AvgNval_Callback(hObject, eventdata, handles)
% hObject    handle to edit_AvgNval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_AvgNval as text
%        str2double(get(hObject,'String')) returns contents of edit_AvgNval as a double
a = str2double(get(hObject,'String'));

if isnan(a) || ~isreal(a) || logical(abs(round(a)-a)<eps)==0 || a < 1
    set(hObject,'String','3');    
    errordlg('You must enter an integer positive value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_AvgNval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_AvgNval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ToleranceForConverg_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ToleranceForConverg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ToleranceForConverg as text
%        str2double(get(hObject,'String')) returns contents of edit_ToleranceForConverg as a double
a = str2double(get(hObject,'String'));

if isnan(a) || ~isreal(a) || a < 0 || a > 1
    set(hObject,'String','1.e-4');
    errordlg('You must enter a real value between 0.0 and 1.0','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_ToleranceForConverg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ToleranceForConverg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_MaxOrderDandB_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MaxOrderDandB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MaxOrderDandB as text
%        str2double(get(hObject,'String')) returns contents of edit_MaxOrderDandB as a double
a = str2double(get(hObject,'String'));

if isnan(a) || ~isreal(a) || logical(abs(round(a)-a)<eps)==0 || a < 1
    set(hObject,'String','10');    
    errordlg('You must enter an integer positive value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_MaxOrderDandB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MaxOrderDandB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
