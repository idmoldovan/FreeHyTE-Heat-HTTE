function varargout = HeatStructDef(varargin)
% HEATSTRUCTDEF MATLAB code for HeatStructDef.fig
%      HEATSTRUCTDEF, by itself, creates a new HEATSTRUCTDEF or raises the existing
%      singleton*.
%
%      H = HEATSTRUCTDEF returns the handle to a new HEATSTRUCTDEF or the handle to
%      the existing singleton*.
%
%      HEATSTRUCTDEF('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HEATSTRUCTDEF.M with the given input arguments.
%
%      HEATSTRUCTDEF('Property','Value',...) creates a new HEATSTRUCTDEF or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HeatStructDef_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HeatStructDef_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HeatStructDef

% Last Modified by GUIDE v2.5 07-Oct-2016 14:13:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HeatStructDef_OpeningFcn, ...
                   'gui_OutputFcn',  @HeatStructDef_OutputFcn, ...
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


% UIWAIT makes HeatStructDef wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = HeatStructDef_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% OPENING FUNCTION
% * Executes just before HeatStructDef is made visible;
% * Reads the previous data in the local |mat| file and fills in the
% fields.
function HeatStructDef_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HeatStructDef (see VARARGIN)

%%
% If no varargin was defined, it is set to zero and saved to handles
if isempty(varargin)
    handles.varargin = 0;
else
    handles.varargin = varargin{1};
end

%%
% Setting warnings to off. These warnings are caused by missing files and
% variables before the local |mat| files are written for the first time and
% by the possibility that the problem is purely Neumann or Dirichlet. The
% warnings are re-activated after the successful execution, at the end of
% the |main.m| function.
warning('off','MATLAB:DELETE:FileNotFound');
warning('off','MATLAB:load:variableNotFound');


%%
% Creates the |handles| structure
% Choose default command line output for HeatStructDef
handles.output = hObject;

%% 
% Getting the current working folder (in R2012 may be different from the
% folder you started the app from!)
handles.WorkingFolder = pwd;

%%
% If there exists a local |mat| file, it loads its key elements to fill in
% the fields with data taken from the previous run.
if exist('./HeatStructDef.mat','file')
    load('HeatStructDef','L','B','Nx','Ny','EdgesOrder','LoopsOrder',...
        'NumberGaussPoints','RunOption','k','Q','MeshOption');
    
    %%
    % In rare situations, no values are stored for L, B, Nx, Ny. If this is
    % the case, it reads NaN and thus ChangesQ always results 1. To avoid
    % this, when NaN is red, it is substituted by zero.
    if isnan(L)
        L=0; B=0; Nx=0; Ny=0;
    end
 
    set(handles.edit_DimX,'String',sprintf('%d',L));
    set(handles.edit_DimY,'String',sprintf('%d',B));
    set(handles.edit_NLoopX,'String',sprintf('%d',Nx));
    set(handles.edit_NLoopY,'String',sprintf('%d',Ny));
    set(handles.edit_OrderEdge,'String',sprintf('%d',EdgesOrder));
    set(handles.edit_OrderLoop,'String',sprintf('%d',LoopsOrder));
    set(handles.edit_NGP,'String',sprintf('%d',NumberGaussPoints));
    set(handles.edit_Conductivity,'String',sprintf('%g',k));
    set(handles.edit_InternalHeat,'String',sprintf('%g',Q));
    set(handles.popupmenu_adaptive,'Value',RunOption);
    set(handles.popupmenu_mesh,'Value',MeshOption);
    
    
    %%
    % If |MeshOption = 2|, that is, the mesh generation is automatic, it
    % makes no sense to edit the fields associated to the regular mesh
    % generator. They become inactive.
    if MeshOption == 1
        set(handles.edit_DimX, 'enable', 'on');
        set(handles.edit_DimY, 'enable', 'on');
        set(handles.edit_NLoopX, 'enable', 'on');
        set(handles.edit_NLoopY, 'enable', 'on');
    else
        set(handles.edit_DimX, 'enable', 'off');
        set(handles.edit_DimY, 'enable', 'off');
        set(handles.edit_NLoopX, 'enable', 'off');
        set(handles.edit_NLoopY, 'enable', 'off');
    end
    
end


% Update handles structure
guidata(hObject, handles);



%% NEXT BUTTON
% * Executes on button press in |pushbutton_next|;
% * It recovers all data provided by the user in the GUI fields;
% * It computes the |edgesArray| list (consisting of the exterior edges of the structure);
% * Iit checks for relevant (i.e. mesh) changes as compared to the previous |mat| file. If
%   such changes exist, it deletes the |mat| file corresponding to the next
%   GUIs to force their definition from scratch. If no mesh changes were
%   detected, the next GUI will load its previous version;
% * If the user asked for the model to be saved, it saves the information
%   regarding this GUI in the specified file and folder.
function pushbutton_next_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Reading the information provided by the user in the GUI fields.
EdgesOrder = str2double(get(handles.edit_OrderEdge,'String'));
LoopsOrder = str2double(get(handles.edit_OrderLoop,'String'));
NumberGaussPoints = str2double(get(handles.edit_NGP,'String'));
RunOption = get(handles.popupmenu_adaptive,'Value');
MeshOption = get(handles.popupmenu_mesh,'Value');
k = str2double(get(handles.edit_Conductivity,'String'));
Q = str2double(get(handles.edit_InternalHeat,'String'));


%%
% If the user asked for the model to be saved, it stores the path and the
% file name.
if isfield(handles,'DirName')
    DirName = handles.DirName;
    FileName = handles.FileName;
else
    DirName = '';
    FileName = '';
end

%%
% *Procedure for the regular rectangular mesh*
if MeshOption == 1
    
    %%
    % Reading mesh information
    L = str2double(get(handles.edit_DimX,'String'));
    B = str2double(get(handles.edit_DimY,'String'));
    Nx = str2double(get(handles.edit_NLoopX,'String'));
    Ny = str2double(get(handles.edit_NLoopY,'String'));
    
    %%
    % Computing the |edgesArray| list
    L1 = linspace(1,Ny,Ny);
    L2 = linspace((Nx*Ny)+1,((Nx*Ny)+1)+(Ny-1),Ny);
    L3 = linspace(((Nx*Ny)+1)+Ny,((Nx*Ny)+1)+Ny+((Ny+1)*(Nx-1)),Nx);
    L4 = linspace(((Nx*Ny)+1)+2*Ny,((Nx*Ny)+1)+Ny+((Ny+1)*(Nx-1))+Ny,Nx);
    edgesArray = sort(cat(1,L1',L2',L3',L4'));
    
    %%
    % |p| and |t| variables are allocated dummy values to avoid errors
    % related to the comparison between |mat| files corresponding to
    % regular (rectangular) and triangular meshes.
    p = 0; t = 0;
    
    %%
    % Check if critical changes were made in the current GUI session. Note
    % that critical changes only refer to the mesh, not necessarily to the
    % geometry of the structure.
    if exist('./HeatStructDef.mat','file')
        save('ToDetectChanges','MeshOption','Nx','Ny','p','t'); % temporary file to detect critical changes
        % Fields whose change triggers the change of the mesh
        CriticalFields = {'MeshOption','Nx','Ny','p','t'};
        % Checks the old |HeatStructDef| and the new |ToDetectChanges| files
        % for changes in the critical fields. ChangesQ = 1 if there are
        % changes, 0 otherwise.
        ChangesQ = any(~isfield(comp_struct(load('HeatStructDef.mat'),...
            load('ToDetectChanges.mat')),CriticalFields));
        % deleting the auxiliary file
        delete('ToDetectChanges.mat');
    else
        % If there is no |HeatStructDef| file in the first place, it sets
        % ChangesQ to 1.
        ChangesQ = 1;
    end
    
    %%
    % Saving the workspace to the local |mat| file. If requested by the
    % user, the GUI information is also saved to the file specified by the
    % user.
    save('HeatStructDef','L','B','Nx','Ny','EdgesOrder','LoopsOrder',...
        'NumberGaussPoints','k','Q','p','t','edgesArray','DirName','FileName',...
        'RunOption','MeshOption','ChangesQ');
    
    % If the folder where it looks for pre-load files is different from the
    % current folder, it creates a copy of HeatStructDef in the former
    if  ~strcmp(pwd,handles.WorkingFolder)
        save(fullfile(handles.WorkingFolder,'HeatStructDef'),'L','B',...
            'Nx','Ny','EdgesOrder','LoopsOrder',...
            'NumberGaussPoints','k','Q','p','t','edgesArray',...
            'DirName','FileName','RunOption','MeshOption','ChangesQ');
    end
    
    % Saves to the save file
    if ~isempty(DirName)
        save(fullfile(DirName,FileName),...
            'L','B','Nx','Ny','EdgesOrder','LoopsOrder',...
            'NumberGaussPoints','k','Q','k','Q','p','t','edgesArray',...
            'RunOption','MeshOption','ChangesQ');
    end
    
    %%
    % Closing the GUI
    close(handles.figure1);
    
    %%
    % If there are relevant changes, it physically deletes the |mat| files
    % associated to the following GUIs.
    if ChangesQ
        delete('HeatRegBC1.mat','HeatRegBC2.mat');
    end
    
    HeatRegBC1;
    
    %%
    % *Procedure for the automatic triangular mesh*
else
    
    %%
    % Reading mesh information. This data is useless for the automatic mesh
    % generation. It is only stored to fill in the corresponding fields in
    % a future run.
    L = str2double(get(handles.edit_DimX,'String'));
    B = str2double(get(handles.edit_DimY,'String'));
    Nx = str2double(get(handles.edit_NLoopX,'String'));
    Ny = str2double(get(handles.edit_NLoopY,'String'));
    %%
    % |edgesArray| variable is allocated a dummy value to avoid errors
    % related to the comparison between |mat| files corresponding to
    % regular (rectangular) and triangular meshes.
    edgesArray = 0; 
    %%
    % Launching |pdetool| to define the mesh
    pdewindow = pdeinit; 
    %% 
    % Stopping the execution until the |pdetool| window is closed 
    waitfor(pdewindow);
    
    %%
    % This is a basic check to confirm that the user saved the mesh 
    % information. A new mesh need not be created if the corresponding
    % information already exist in |HeatStructDef.mat|, so the program checks
    % if a new mesh was defined or if |HeatStructDef.mat| exists. Of course,
    % the check fails if |HeatStructDef.mat| exists, but does not contain
    % relevant mesh information (the program exists with an error).
    BaseVars = evalin('base','whos'); % collects variable info from base
    % if the mesh variables are not found in base and an old definition
    % does not exist in HeatStructDef.mat
    if (~ismember('p',[BaseVars(:).name]) || ~ismember('t',[BaseVars(:).name])) && ...
            (~exist('./HeatStructDef.mat','file'))
        errordlg('No mesh information was exported or variable names were changed. Please press "Next" again and export the mesh info as p, e and t.','Invalid input','modal');
        uicontrol(hObject);
        return;
    end
    
    %%
    % if the mesh variables are found in base, it loads them...
    if ismember('p',[BaseVars(:).name]) && ismember('t',[BaseVars(:).name])
        p = evalin('base','p');
        t = evalin('base','t');
    %%
    % ... otherwise, it loads them from the HeatStructDef.mat
    else
        load('HeatStructDef','p','t');
    end
    
    %%
    % Check if critical changes were made in the current GUI session. Note
    % that critical changes only refer to the mesh, not necessarily to the
    % geometry of the structure.
    if exist('./HeatStructDef.mat','file')
        save('ToDetectChanges','MeshOption','Nx','Ny','p','t'); % temporary file to detect critical changes
        % Fields whose change triggers the change of the mesh
        CriticalFields = {'MeshOption','Nx','Ny','p','t'};
        % Checks the old |HeatStructDef| and the new |ToDetectChanges| files
        % for changes in the critical fields. ChangesQ = 1 if there are
        % changes, 0 otherwise.
        ChangesQ = any(~isfield(comp_struct(load('HeatStructDef.mat'),...
            load('ToDetectChanges.mat')),CriticalFields));
        % deleting the auxiliary file
        delete('ToDetectChanges.mat');
    else
        % If there is no |HeatStructDef| file in the first place, it sets
        % ChangesQ to 1.
        ChangesQ = 1;
    end
    
    %%
    % Saving the workspace to the local |mat| file. If requested by the
    % user, the GUI information is also saved to the file specified by the
    % user.
    save('HeatStructDef','L','B','Nx','Ny','EdgesOrder','LoopsOrder',...
        'NumberGaussPoints','k','Q','p','t','edgesArray','DirName','FileName',...
        'RunOption','MeshOption','ChangesQ');
    
    % If the folder where it looks for pre-load files is different from the
    % current folder, it creates a copy of HeatStructDef in the former
    if  ~strcmp(pwd,handles.WorkingFolder)
        save(fullfile(handles.WorkingFolder,'HeatStructDef'),'L','B',...
            'Nx','Ny','EdgesOrder','LoopsOrder',...
            'NumberGaussPoints','k','Q','p','t','edgesArray',...
            'DirName','FileName','RunOption','MeshOption','ChangesQ');
    end
    
    % Saves to the save file
    if ~isempty(DirName)
        save(fullfile(DirName,FileName),...
            'L','B','Nx','Ny','EdgesOrder','LoopsOrder',...
            'NumberGaussPoints','k','Q','k','Q','p','t','edgesArray',...
            'RunOption','MeshOption','ChangesQ');
    end
    
    
    %%
    % Closing the GUI
    close(handles.figure1);
    
    %%
    % If there are relevant changes, it physically deletes the |mat| files
    % associated to the following GUIs.
    if ChangesQ
        delete('HeatTriBC1.mat','HeatTriBC2.mat');
    end
    
    HeatTriBC1;
    
end


%% RESET BUTTON
% * Executes on button press in |pushbutton_reset|;
% * It deletes all entries of the edit fields and resets the pop-up menus;
% * It also deletes the save information so that the save button action is
% reversed.
function pushbutton_reset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.edit_DimX,'String','');
set(handles.edit_DimY,'String','');
set(handles.edit_NLoopX,'String','');
set(handles.edit_NLoopY,'String','');
set(handles.edit_OrderEdge,'String','');
set(handles.edit_OrderLoop,'String','');
set(handles.edit_NGP,'String','');
set(handles.edit_Conductivity,'String','');
set(handles.edit_InternalHeat,'String','');
set(handles.edit_Path,'String','');
set(handles.popupmenu_adaptive,'Value',1);
set(handles.popupmenu_mesh,'Value',1);


%%
% Reseting |edit_Path|, the pop-up menus, and the properties of the mesh
% definition fields
set(handles.edit_Path,'String',...
    'THE MODEL WILL NOT BE SAVED!  Please press Save if you wish to save the model.',...
    'BackgroundColor',[1.0 0.694 0.392]);

set(handles.edit_DimX, 'enable', 'on');
set(handles.edit_DimY, 'enable', 'on');
set(handles.edit_NLoopX, 'enable', 'on');
set(handles.edit_NLoopY, 'enable', 'on');

handles.FileName = '';
handles.DirName = '';

%%
% Updating the |handles| structure
guidata(hObject, handles);


%% SAVE BUTTON
% * Executes on button press in |Save|;
% * Reads the path and filename provided by the user to save the model;
% * Generates the |FileName| and |DirName| members in the |handles|
% structure to store this information and to make it available to other
% functions.
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% 
% Reads the |handles| structure
handles=guidata(hObject);

%%
% Gets the path and filename to save the model
[FileName,DirName] = uiputfile('*.mat','Save as');

if FileName
    %%
    % Registers the path to the save file in the |edit_Path| box
    set(handles.edit_Path,'string',fullfile(DirName,FileName),'BackgroundColor',[0 1 0]);
    
    %%
    % Generates the |FileName| and |DirName| members in the |handles|
    % structure to store this information and to make it available to other
    % functions.
    handles.FileName = FileName;
    handles.DirName = DirName;
    
    %%
    % If user cancelled the saving process and no save file was given
else
    handles.FileName = '';
    handles.DirName = '';
    
end

%%
% Updating the |handles| structure
guidata(hObject, handles);



%% LOAD BUTTON
% * Executes on button press in |Load|;
% * Loads all relevant variables in the file appointed by the user into the
% workspace. If the definition of the model was not completed in the
% previous session, some of the variables cannot be loaded (the
% corresponding warning is suppressed);
% * Stores the required variables into the HeatStructDef |mat| file.
% This data must always be available;
% * Verifies if the data generated by the other three GUIs are
% available. If they are, stores them into the local |mat| files, otherwise
% _deletes_ the |mat| files, to force their reinitialization;
% * Deletes the |FileName| and |DirName| variables to avoid overwriting of the
% loaded file;
% * Refreshes the interface, filling in the fields with the newly loaded
% values.
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%%
% Reads the |handles| structure
handles=guidata(hObject);

%%
% Loads the |mat| file indicated by the user
[FileName,DirName] = uigetfile('*.mat','File to load');

%%
% Deletes the |FileName| and |DirName| variables to avoid overwriting;
handles.FileName = '';
handles.DirName = '';

%%
% Loading all relevant variables into the current workspace
load(fullfile(DirName,FileName),'L','B','Nx','Ny','EdgesOrder','LoopsOrder',...
    'NumberGaussPoints','k','Q','p','t','MeshOption','RunOption',...
    'ChangesQ','edgesArray','data','dataTemp','dataHeat','edgesDirichlet',...
    'edgesNeumann');

% If RunOption is adaptive, loads the adaptive data
if RunOption == 2
    load(fullfile(DirName,FileName),'SelectionCriterion','SelectionTol',...
        'thresh','MaxOutlierIter','MinIter','AvgNVal',...
        'TargetErrorNorm','MaxOrder','StoppingCriterion');
end

%%
% Saving the local |HeatStructDef.mat| file
save('HeatStructDef','L','B','Nx','Ny','EdgesOrder','LoopsOrder',...
    'NumberGaussPoints','k','Q','p','t','RunOption','MeshOption',...
    'edgesArray','ChangesQ');

% If the folder where it looks for pre-load files is different from the
% current folder, it creates a copy of StructDef in the former
if  ~strcmp(pwd,handles.WorkingFolder)
    save(fullfile(handles.WorkingFolder,'HeatStructDef'),'L','B',...
        'Nx','Ny','EdgesOrder','LoopsOrder','NumberGaussPoints','k',...
        'Q','p','t','RunOption','MeshOption','edgesArray','ChangesQ');
end

%%
% If 'p' exists (is not zero), updates 'p' and 't' in the base workspace
if any(any(p~=0))
    assignin('base','p',p);
    assignin('base','t',t);
end

%%
% Depending on the kind of mesh that is selected ...
if MeshOption == 1 % regular rectangular mesh
    %%
    % ... checks if the variable created by the second GUI exists and if it
    % doesn't, it deletes the local |mat| file to force reinitialization...
    if ~exist('data','var')
        delete('HeatRegBC1.mat');
        %%
        % ... or it stores (and overwrites) the |mat| file if the variable
        % exists.
    else
        save('HeatRegBC1','data');
    end
    
    %%
    % ... checks if the variables created by the third GUI exist and if they
    % don't, it deletes the local |mat| file to force reinitialization...
    if ~exist('dataTemp','var') && ~exist('dataHeat','var')
        delete('HeatRegBC2.mat');
        %%
        % ... or it stores (and overwrites) the |mat| file if the variables
        % exist.
    else
        delete('HeatRegBC2.mat');
        dummy = 0;
        save('HeatRegBC2','dummy');
        if exist('dataTemp','var')
            save('HeatRegBC2','-append','edgesDirichlet','dataTemp');
        end
        if exist('dataHeat','var')
            save('HeatRegBC2','-append','edgesNeumann','dataHeat');
        end
    end
    
    %%
    % Same operation for the irregular mesh file.
else
    if ~exist('data','var')
        delete('HeatTriBC1.mat');
    else
        save('HeatTriBC1','data');
    end
    
    if ~exist('dataTemp','var') && ~exist('dataHeat','var')
        delete('HeatTriBC2.mat');
    else
        delete('HeatTriBC2.mat');
        dummy = 0;
        save('HeatTriBC2','dummy');
        if exist('dataTemp','var')
            save('HeatTriBC2','-append','edgesDirichlet','dataTemp');
        end
        if exist('dataHeat','var')
            save('HeatTriBC2','-append','edgesNeumann','dataHeat');
        end
    end
end

%%
% If RunOption is adaptive, checks if a variable created by the fourth GUI
% exists and if it doesn't, it deletes the local |mat| file to force
% reinitialization...
if RunOption == 2
    if ~exist('SelectionCriterion','var')
        delete('Adaptive.mat');
        %%
        % ... or it stores (and overwrites) the |mat| file if the variable
        % exists.
    else
        save('Adaptive','SelectionCriterion',...
            'SelectionTol','thresh','MaxOutlierIter','MinIter','AvgNVal',...
            'TargetErrorNorm','MaxOrder','StoppingCriterion');
    end
end

%%
% Refreshes the interface, writing the loaded values into its fields...
set(handles.edit_DimX,'String',sprintf('%d',L));
set(handles.edit_DimY,'String',sprintf('%d',B));
set(handles.edit_NLoopX,'String',sprintf('%d',Nx));
set(handles.edit_NLoopY,'String',sprintf('%d',Ny));
set(handles.edit_OrderEdge,'String',sprintf('%d',EdgesOrder));
set(handles.edit_OrderLoop,'String',sprintf('%d',LoopsOrder));
set(handles.edit_NGP,'String',sprintf('%d',NumberGaussPoints));
set(handles.edit_Conductivity,'String',sprintf('%g',k));
set(handles.edit_InternalHeat,'String',sprintf('%g',Q));
set(handles.popupmenu_adaptive,'Value',RunOption);
set(handles.popupmenu_mesh,'Value',MeshOption);
set(handles.edit_Path,'string',...
    'THE MODEL WILL NOT BE SAVED!  Please press Save if you wish to save the model.',...
    'BackgroundColor',[1.0 0.694 0.392]);

%%
% ... and activates or inactivates the regular mesh fields according to the
% type of mesh in the loaded model.
if MeshOption == 1
    set(handles.edit_DimX, 'enable', 'on');
    set(handles.edit_DimY, 'enable', 'on');
    set(handles.edit_NLoopX, 'enable', 'on');
    set(handles.edit_NLoopY, 'enable', 'on');
else
    set(handles.edit_DimX, 'enable', 'off');
    set(handles.edit_DimY, 'enable', 'off');
    set(handles.edit_NLoopX, 'enable', 'off');
    set(handles.edit_NLoopY, 'enable', 'off');
end

%%
% Updates the |handles| structure
guidata(hObject, handles);



%% CLEAR BUTTON
% * Executes on button press in |Clear|;
% * Deletes the path and name of the save file and reinitializes the
% |edit_Path| field.
function Clear_Callback(hObject, eventdata, handles)
% hObject    handle to Clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Reads the |handles| structure
handles=guidata(hObject);

%%
% Deletes the |FileName| and |DirName| variables
handles.FileName = '';
handles.DirName = '';
set(handles.edit_Path,'string',...
    'THE MODEL WILL NOT BE SAVED!  Please press Save if you wish to save the model.',...
    'BackgroundColor',[1.0 0.694 0.392]);

%%
% Updates the |handles| structure
guidata(hObject, handles);


%% POP-UP MENUS
% *Adaptive pop-up menu*
% * Controls whether a single-run or an adaptive analysis is required
% * Executes on selection change in popupmenu_adaptive.
function popupmenu_adaptive_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_adaptive (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_adaptive contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_adaptive
% Save the handles structure.
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_adaptive_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_adaptive (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_mesh.
function popupmenu_mesh_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_mesh contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_mesh
if get(handles.popupmenu_mesh,'Value') == 1
    set(handles.edit_DimX, 'enable', 'on');
    set(handles.edit_DimY, 'enable', 'on');
    set(handles.edit_NLoopX, 'enable', 'on');
    set(handles.edit_NLoopY, 'enable', 'on');
else
   set(handles.edit_DimX, 'enable', 'off');
   set(handles.edit_DimY, 'enable', 'off');
   set(handles.edit_NLoopX, 'enable', 'off');
   set(handles.edit_NLoopY, 'enable', 'off');
end

% --- Executes during object creation, after setting all properties.
function popupmenu_mesh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





%% EDIT FIELDS
% Field box editing functions. Most Callbacks check for the validity of the
% data.


function edit_DimX_Callback(hObject, eventdata, handles)
% hObject    handle to edit_DimX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_DimX as text
%        str2double(get(hObject,'String')) returns contents of edit_DimX as a double
L = str2double(get(hObject,'String'));
if isnan(L) || ~isreal(L) || isequal(L,0) || L < 0
    set(hObject,'String','');
    errordlg('You must enter a real positive value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_DimX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_DimX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_DimY_Callback(hObject, eventdata, handles)
% hObject    handle to edit_DimY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_DimY as text
%        str2double(get(hObject,'String')) returns contents of edit_DimY as a double
B = str2double(get(hObject,'String'));
if isnan(B) || ~isreal(B) ||  isequal(B,0) || B < 0
    set(hObject,'String','');
    errordlg('You must enter a real positive value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_DimY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_DimY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_NLoopX_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NLoopX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NLoopX as text
%        str2double(get(hObject,'String')) returns contents of edit_NLoopX as a double
Nx = str2double(get(hObject,'String'));
if isnan(Nx) || ~isreal(Nx) || logical(abs(round(Nx)-Nx)<eps)==0 || isequal(Nx,0) || Nx < 0
    set(hObject,'String','');
    errordlg('You must enter a positive integer value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_NLoopX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NLoopX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_NLoopY_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NLoopY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NLoopY as text
%        str2double(get(hObject,'String')) returns contents of edit_NLoopY as a double
Ny = str2double(get(hObject,'String'));
if isnan(Ny) || ~isreal(Ny) || logical(abs(round(Ny)-Ny)<eps)==0 || isequal(Ny,0) || Ny < 0
    set(hObject,'String','');
    errordlg('You must enter a positive integer value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_NLoopY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NLoopY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_OrderEdge_Callback(hObject, eventdata, handles)
% hObject    handle to edit_OrderEdge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_OrderEdge as text
%        str2double(get(hObject,'String')) returns contents of edit_OrderEdge as a double
EdgesOrder = str2double(get(handles.edit_OrderEdge,'String'));
if isnan(EdgesOrder) || ~isreal(EdgesOrder) || ...
        logical(abs(round(EdgesOrder)-EdgesOrder)<eps)==0 || EdgesOrder < 0
    set(hObject,'String','');
    errordlg('You must enter either 0 or a natural number','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_OrderEdge_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_OrderEdge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_OrderLoop_Callback(hObject, eventdata, handles)
% hObject    handle to edit_OrderLoop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_OrderLoop as text
%        str2double(get(hObject,'String')) returns contents of edit_OrderLoop as a double
LoopsOrder = str2double(get(handles.edit_OrderLoop,'String'));
if isnan(LoopsOrder) || ~isreal(LoopsOrder) || ...
        logical(abs(round(LoopsOrder)-LoopsOrder)<eps)==0 || LoopsOrder < 0
    set(hObject,'String','');
    errordlg('You must enter either 0 or a natural number','Invalid input','modal');
    uicontrol(hObject);
    return;
end


% --- Executes during object creation, after setting all properties.
function edit_OrderLoop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_OrderLoop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_NGP_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NGP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NGP as text
%        str2double(get(hObject,'String')) returns contents of edit_NGP as a double
NumberGaussPoints = str2double(get(handles.edit_NGP,'String'));
if isnan(NumberGaussPoints) || ~isreal(NumberGaussPoints) || ...
        logical(abs(round(NumberGaussPoints)-NumberGaussPoints)<eps)==0 ||...
        isequal(NumberGaussPoints,0) || NumberGaussPoints < 0
    set(hObject,'String','');
    errordlg('You must enter a positive integer value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_NGP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NGP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_InternalHeat_Callback(hObject, eventdata, handles)
% hObject    handle to edit_InternalHeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_InternalHeat as text
%        str2double(get(hObject,'String')) returns contents of edit_InternalHeat as a double
Q = str2double(get(handles.edit_InternalHeat,'String'));
if isnan(Q) || ~isreal(Q)
    set(hObject,'String','');
    errordlg('You must enter real numeric value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_InternalHeat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_InternalHeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Conductivity_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Conductivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Conductivity as text
%        str2double(get(hObject,'String')) returns contents of edit_Conductivity as a double
K = str2double(get(handles.edit_Conductivity,'String'));
if isnan(K) || ~isreal(K) || isequal(K,0) || K < 0
    set(hObject,'String','');
    errordlg('You must enter a real positive numeric value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_Conductivity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Conductivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_Path_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Path as text
%        str2double(get(hObject,'String')) returns contents of edit_Path as a double

% --- Executes during object creation, after setting all properties.
function edit_Path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





%% TEXT FIELDS
% Text fields' editing functions. The respective |ButtonDownFcn| defines
% the help of that field, accessed through right clicking.


% --- Executes during object creation, after setting all properties.
function text_DimX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_DimX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text_DimY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_DimY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text_NLoopX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_NLoopX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text_NLoopY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_NLoopY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text_OrderEdge_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_OrderEdge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text_OrderLoop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_OrderLoop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text_NGP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_NGP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text_Conductivity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_Conductivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text_InternalHeat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_InternalHeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
