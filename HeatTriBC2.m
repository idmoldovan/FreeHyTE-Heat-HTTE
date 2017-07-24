function varargout = HeatTriBC2(varargin)
% HEATTRIBC2 MATLAB code for HeatTriBC2.fig
%      HEATTRIBC2, by itself, creates a new HEATTRIBC2 or raises the existing
%      singleton*.
%
%      H = HEATTRIBC2 returns the handle to a new HEATTRIBC2 or the handle to
%      the existing singleton*.
%
%      HEATTRIBC2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HEATTRIBC2.M with the given input arguments.
%
%      HEATTRIBC2('Property','Value',...) creates a new HEATTRIBC2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HeatTriBC2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HeatTriBC2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to next (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HeatTriBC2

% Last Modified by GUIDE v2.5 10-Aug-2016 12:06:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @HeatTriBC2_OpeningFcn, ...
    'gui_OutputFcn',  @HeatTriBC2_OutputFcn, ...
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


% UIWAIT makes HeatTriBC2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = HeatTriBC2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% OPENING FUNCTION
% * Executes just before |HeatTriBC2| is made visible;
% * Reads the mesh data from the |mat| file of the previous GUI and
% constructs the lists of Dirichlet and Neumann boundaries;
% * If the |mat| file of the current GUI exists (meaning that the 
% previous GUI was not changed), it loads the boundary information,
% otherwise it sets all boundaries to Dirichlet.

function HeatTriBC2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HeatTriBC2 (see VARARGIN)

%%
% Creates the |handles| structure
% Choose default command line output for |HeatTriBC2|
handles.output = hObject;

%%
% Loads the exterior edges and their types from the previous GUI
load('HeatTriBC1','edgesArray','data','DorN');

%%
% Creates two arrays listing the exterior Dirichlet and Neumann edges
% (|edgesDirichlet| and |edgesNeumann|)
for i=1:length(edgesArray)
    if strcmp(data{i,2},DorN{1}) % if the edge type is Dirichlet
        handles.edgesDirichlet(i,1)=edgesArray(i); % list with the Dirichlet edges
    else
        handles.edgesNeumann(i,1)=edgesArray(i);  % list with the Neumann edges
    end
end

%%
% *Dirichlet boundary conditions*

%%
% Generates the |dataTemp| matrix to store the id of the exterior Dirichlet 
% edges and their enforced boundary conditions. The |dataTemp| matrix is 
% created from from scratch or imported from the |HeatTriBC2| file, if 
% such file exists (meaning that the previous GUI was left unchanged).


if isfield(handles,'edgesDirichlet')==1
    handles.edgesDirichlet(handles.edgesDirichlet==0)=[];
    
    %%
    % Writes all |edgesDirichlet| into the |listboxDir|
    set(handles.listboxDir,'string',handles.edgesDirichlet);
    
    %%
    % Creates the |dataTemp| matrix.
    handles.dataTemp=cell(length(handles.edgesDirichlet),2);
    
    %%
    % If there exists a local |mat| file, it loads it to fill in
    % the fields with data taken from the previous next...
    if exist('./HeatTriBC2.mat','file') % reuse the previous data
        load('HeatTriBC2','dataTemp');
        handles.dataTemp=dataTemp;
    %%
    % ... otherwise it just sets all Dirichlet boundary conditions to zero    
    else
        for i=1:length(handles.edgesDirichlet)
            handles.dataTemp{i,1}=handles.edgesDirichlet(i);
            handles.dataTemp{i,2}='0';
        end
    end
    
    %%
    % Creates the table where the Dirichlet boundary conditions are listed, 
    % along with their types
    column1={'Boundary ID','Temperature'};
    uitable('units','Normalized','Position',[0.482,0.16,0.11,0.45],...
        'Data',handles.dataTemp,'ColumnName',column1,'RowName',[]);
end


%%
% *Neumann boundary conditions*

%%
% Generates the |dataHeat| matrix to store the id of the exterior Neumann 
% edges and their enforced boundary conditions. The |dataHeat| matrix is 
% created from from scratch or imported from the |HeatTriBC2| file, if such 
% file exists (meaning that the previous GUI was left unchanged).

if isfield(handles,'edgesNeumann')==1

    handles.edgesNeumann(handles.edgesNeumann==0)=[];
    
    %%
    % Writes all |edgesNeumann| into the |listboxNeu|
    set(handles.listboxNeu,'string',handles.edgesNeumann);
    
    %%
    % Creates the |dataHeat| matrix.
    handles.dataHeat=cell(length(handles.edgesNeumann),2);
    
    
    %%
    % If there exists a local |mat| file, it loads it to fill in
    % the fields with data taken from the previous next...
    if exist('./HeatTriBC2.mat','file')  % reuse the previous data
        load('HeatTriBC2','dataHeat');
        handles.dataHeat=dataHeat;
    %%
    % ... otherwise it just sets all Neumann boundary conditions to NaN
    else
        for i=1:length(handles.edgesNeumann)
            handles.dataHeat{i,1}=handles.edgesNeumann(i);
            handles.dataHeat{i,2}='0';
        end
    end
    
    %%
    % Creates the table where the Neumann boundary conditions are listed, 
    % along with their types
    column2={'Boundary ID','Heat Flux'};
    uitable('units','Normalized','Position',[0.85,0.16,0.11,0.45],...
        'Data',handles.dataHeat,'ColumnName',column2,'RowName',[]);
end

%%
% Generates the code for drawing the mesh, along with the mesh information
% buttons

%load the mesh data
load('HeatTriBC1','nodes','edges_nodes','edges_loops','loops_nodes',...
    'loops_edges');

nel = size(loops_nodes,1);        % Total Number of Elements in the Mesh
nnel = size(loops_nodes,2);           % Number of nodes per Element
nnode = size(nodes,1);      % Total Number of Nodes in the Mesh
nedge = size(edges_nodes,1); % Total Number of Edges in the Mesh

% For drawing purposes
limxmin = min(nodes(:,1));
limxmax = max(nodes(:,1));
limymin =  min(nodes(:,2));
limymax =  max(nodes(:,2));

%
% Plotting the Finite Element Mesh
% Initialization of the required matrices
X = zeros(nnel,nel) ;
Y = zeros(nnel,nel) ;
% Extract X,Y coordinates for the (iel)-th element
for iel = 1:nel
    X(:,iel) = nodes(loops_nodes(iel,:),1) ;
    Y(:,iel) = nodes(loops_nodes(iel,:),2) ;
end

patch(X,Y,'w')
axis([limxmin-0.01*abs(limxmin) limxmax+0.01*abs(limxmax) limymin-0.01*abs(limymin) limymax+0.01*abs(limymax)]);
axis equal;
axis off ;

% To display Node Numbers % Element Numbers
axpos = getpixelposition(handles.axes3); % position & dimension of the axes object
% Define button's weight and height
bweight = 60;
bheight = 20;
pos = [((axpos(1)+axpos(3))/2) (axpos(2)-1.5*bheight) bweight bheight]; % align the second button with the center of the axes obj limit
ShowNodes = uicontrol('style','toggle','string','Nodes',....
    'position',[(pos(1)-2*bweight) pos(2) pos(3) pos(4)],'background',...
    'white','units','normalized');

ShowEdges = uicontrol('style','toggle','string','Edges',....
    'position',[pos(1) pos(2) pos(3) pos(4)],'background','white',...
    'units','normalized');

ShowElements = uicontrol('style','toggle','string','Elements',....
    'position',[(pos(1)+2*bweight) pos(2) pos(3) pos(4)],'background',...
    'white','units','normalized');

set(ShowNodes,'callback',...
    {@SHOWNODES,ShowEdges,ShowElements,nodes,edges_nodes,loops_nodes,X,Y,nnode,nedge,nel});
set(ShowEdges,'callback',...
    {@SHOWEDGES,ShowNodes,ShowElements,nodes,edges_nodes,loops_nodes,X,Y,nnode,nedge,nel});
set(ShowElements,'callback',....
    {@SHOWELEMENTS,ShowNodes,ShowEdges,nodes,edges_nodes,loops_nodes,X,Y,nnode,nedge,nel});

%%
% Updates the handles structure
guidata(hObject, handles);



%% NEXT FUNCTION
% * Executes on button press in |next|;
% * Reads the boundary condition data, stores it in the local |mat| file
% and starts the next GUI;
% * If the user asked for the model to be saved, it saves the information
% regarding this GUI in the specified file and folder.

function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Recovering the user-defined data 
if isfield(handles,'edgesDirichlet')
    edgesDirichlet=handles.edgesDirichlet;
    dataTemp=handles.dataTemp;
end
if isfield(handles,'edgesNeumann')
    edgesNeumann=handles.edgesNeumann;
    dataHeat=handles.dataHeat;
end


%% 
% Saving the local data to save files
% If HeatTriBC2.mat does not exist, it creates it
if ~exist('./HeatTriBC2.mat','file')
    dummy = 0;
    save('HeatTriBC2','dummy');
end

%%
% If requested by the user, appends the local data to the save file
load('HeatStructDef','DirName','FileName');

if isfield(handles,'edgesDirichlet')
    if ~isempty(DirName)
        save(fullfile(DirName,FileName),'-append',...
            'edgesDirichlet','dataTemp');
    end
    save('HeatTriBC2','-append','edgesDirichlet','dataTemp');
end
if isfield(handles,'edgesNeumann')
    if ~isempty(DirName)
        save(fullfile(DirName,FileName),'-append',...
            'edgesNeumann','dataHeat');
    end
    save('HeatTriBC2','-append','edgesNeumann','dataHeat');
end  

%%
% Closes everything and launches whatever comes next
close(handles.figure1) % closing the GUI window

% Trying to close Visualize if it was opened
try
    close('Visualize');
catch
end

load('HeatStructDef','RunOption'); % checking if the run is adaptive

if RunOption == 2
    AdaptiveTri;
else
    CheckHeatTri;
end


%% ASSIGN FLUX FUNCTION
% * Executes on button press in |assignHeatflux|;
% * Fills in the enforced flux table for the boundaries selected in
% |listboxNeu| with the string defined in |editHeat|.

function assignHeatflux_Callback(hObject, eventdata, handles)
% hObject    handle to assignHeatflux (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Gets the selected boundaries in |listboxNeu| 
itemlist = get(handles.listboxNeu,'Value'); % list of selected items in the listbox
nitems = length(itemlist); % number of selected items in the listbox

%%
% Writes the flux definition in the column of |dataHeat|.
for ii = 1:nitems
    crtitem = itemlist(ii);
    handles.dataHeat{crtitem,2}=get(handles.editHeat,'string');
end

%%
% Updates the handles structure
guidata(hObject,handles);

%%
% Redraws the table where the Neumann boundary conditions are listed
column2={'Boundary ID','Heat Flux'};
uitable('units','Normalized','Position',[0.85,0.16,0.11,0.45],...
    'Data',handles.dataHeat,'ColumnName',column2,'RowName',[]);

%% ASSIGN TEMPERATURE FUNCTION
% * Executes on button press in |assignTemp|;
% * Fills in the enforced temperature table for the boundaries selected in
% |listboxDir| with the string defined in |editTemp|.

function assignTemp_Callback(hObject, eventdata, handles)
% hObject    handle to assignTemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Gets the selected boundaries in |listboxDir|
itemlist = get(handles.listboxDir,'Value'); % list of selected items in the listbox
nitems = length(itemlist); % number of selected items in the listbox

%%
% Writes the temperature definition in the column of |dataTemp|.
for ii = 1:nitems
    crtitem = itemlist(ii);
    handles.dataTemp{crtitem,2}=get(handles.editTemp,'string');
end

%%
% Updates the handles structure
guidata(hObject,handles);

%%
% Redraws the table where the Dirichlet boundary conditions are listed
column1={'Boundary ID','Temperature'};
uitable('units','Normalized','Position',[0.482,0.16,0.11,0.45],...
    'Data',handles.dataTemp,'ColumnName',column1,'RowName',[]);


%% RESET NEUMANN/DIRICHLET BOUNDARY CONDITION FUNCTIONS

%%
% *Reset Neumann*
%
% * Executes on button press in |resetNeumann|;
% * Substitutes all previous definitions in |dataHeat| by zero;
% * Redraws the force boundary condition table.

function resetNeumann_Callback(hObject, eventdata, handles)
% hObject    handle to resetNeumann (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Substitutes all previous definitions in |dataNeu| by zero
for i=1:length(handles.edgesNeumann)
    handles.dataHeat{i,2} = '0';
end

%%
% Updates the handles structure
guidata(hObject,handles);

%%
% Redraws the table where the Neumann boundary conditions are listed
column2={'Boundary ID','Heat Flux'};
uitable('units','Normalized','Position',[0.85,0.16,0.11,0.45],...
    'Data',handles.dataHeat,'ColumnName',column2,'RowName',[]);



%%
% *Reset Dirichlet*
%
% * Executes on button press in |resetDirichlet|;
% * Substitutes all previous definitions in |dataTemp| by zero;
% * Redraws the force boundary condition table.

function resetDirichlet_Callback(hObject, eventdata, handles)
% hObject    handle to resetDirichlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Substitutes all previous definitions in |dataTemp| by zero
for i=1:length(handles.edgesDirichlet)
    handles.dataTemp{i,2} = '0';
end

%%
% Updates the handles structure
guidata(hObject,handles);

%%
% Redraws the table where the Dirichlet boundary conditions are listed
column1={'Boundary ID','Temperature'};
uitable('units','Normalized','Position',[0.482,0.16,0.11,0.45],...
    'Data',handles.dataTemp,'ColumnName',column1,'RowName',[]);

%% PREVIOUS FUNCTION
% * Executes on button press in |previous|;
% * Just closes the current GUI and launches the previous one. All changes
% made in the current GUI are lost.

function previous_Callback(hObject, eventdata, handles)
% hObject    handle to previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(HeatTriBC2);

% Trying to close Visualize if it was opened
try
    close('Visualize');
catch
end

HeatTriBC1;


%% ENLARGE FUNCTION
% * Executes on button press in |pushbutton_enlarge|;
% * Just launches the Visualize GUI.

% --- Executes on button press in pushbutton_enlarge.
function pushbutton_enlarge_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_enlarge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Visualize;





%% GUI ENTITIES GENERATION CODE
% * Automatically generated code for the buttons and menus;
% * Some (not-so-sound) checks are performed on the data inserted by the
% user in |editDisp| and |editNeu| just to make sure they are
% mathematically legible. 

% --- Executes on selection change in listboxNeu.
function listboxNeu_Callback(hObject, eventdata, handles)
% hObject    handle to listboxNeu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxNeu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxNeu


% --- Executes during object creation, after setting all properties.
function listboxNeu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxNeu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editHeat_Callback(hObject, eventdata, handles)
% hObject    handle to editHeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editHeat as text
%        str2double(get(hObject,'String')) returns contents of editHeat as a double

[HeatFlux, status] = str2num(get(hObject,'string'));
if any(isnan(HeatFlux)) || ~status  % if the input is something else than
                                     % a vector of reals
    set(hObject,'String','');
    errordlg('All fluxes must have real value','Invalid input','modal');
    uicontrol(hObject);
    return;
end


% --- Executes during object creation, after setting all properties.
function editHeat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editHeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in listboxDir.
function listboxDir_Callback(hObject, eventdata, handles)
% hObject    handle to listboxDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxDir contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxDir


% --- Executes during object creation, after setting all properties.
function listboxDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editTemp_Callback(hObject, eventdata, handles)
% hObject    handle to editTemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTemp as text
%        str2double(get(hObject,'String')) returns contents of editTemp as a double
[T, status] = str2num(get(hObject,'string'));
if any(isnan(T)) || ~status  % if the input is something else than
                                     % a vector of reals
    set(hObject,'String','');
    errordlg('All temperatures must have real value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function editTemp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
