function [NGP, Nodes, Edges, Loops, BConds] = InputProcReg
% INPUTPROCREG is the input processing function for regular meshes.
%
% INPUTPROCREG is called by MAINREG. It reads the input data inserted in 
% the GUIs and organizes it in the data structures Edges, Loops and BConds.
%
%
% BIBLIOGRAPHY
% 1. FreeHyTE Page - https://sites.google.com/site/ionutdmoldovan/freehyte
% 2. Moldovan ID, Cismasiu I - FreeHyTE: theoretical bases and developer’s 
% manual, https://drive.google.com/file/d/0BxuR3pKS2hNHTzB2N2Q4cXZKcGc/view
% 3. FreeHyTE Heat HTTE User's Manual - 
%    https://drive.google.com/drive/folders/0BxuR3pKS2hNHaFhiSjZHOE9TMzg
% 4. Geraldes, MRCB,  Elementos finitos híbridos-Trefftz adaptativos para
%   problemas de condução de calor. MSc Thesis, Universidade Nova de Lisboa,
%   2016 (in Portuguese).
%
%
% INPUT DATA (read from the *.mat files created by the GUIs):
% * from HeatStructDef: k thermal conductivity coefficient  and Q internal 
% generated heat, NumberGaussPoints (the nuber of Gauss-Legendre points for the side
% integration), EdgesOrder and LoopsOrder (the orders of the approximation
% bases on the essential edges of the mesh and in the finite elements);
% * from HeatTriBC1: nodes, edges_nodes, edges_loops, loops_nodes and 
% loops_edges. These variables contain geometrical and topological 
% information regarding the mesh. For a full description of these
% variables, please consult Section 5.3 of reference [2];
% * from HeatTriBC2: edgesDirichlet, edgesNeumann, dataTemp, dataHeat;
% * edgesDirichlet (edgesNeumann): list with the Dirichlet (Neumann) edges;
% * dataTemp (dataHeat): cell array with three columns and as many lines as
% Dirichlet (Neumann) edges. For each edge, it stores the edge's index, 
% and the Dirichlet (Neumann) boundary conditions in the normal and 
% tangential directions. The boundary conditions are stored as strings. For
% further details regarding the definition of the boundary conditions,
% please refer to Section 4.6 of reference [3].
%
%
% OUTPUT DATA (to function MAINREG):
% * NGP is the number of Gauss points for the line integration;
% * Nodes is a (NNODE x 2) matrix, where NNODE is the number of nodes in 
% mesh. It stores the coordinates of each node;
% * Edges, Loops and BConds are data structures storing information
% on the edges, finite elements (loops) and boundary conditions,
% respectively. They are documented in Section 5.3 of reference [2];


%% MESH DATA
% Creation of the Nodes, Edges and Loops data structures

% Loading mesh information from STRUCTDEF GUI
load('HeatRegBC1','nodes','edges_nodes','edges_loops','loops_nodes',...
     'loops_edges');

% Definition of the mesh-related data Nodes, and data structures Edges and
% Loops. For a full description of these data structure, please refer to
% Section 5.3 of reference [2].
Nodes = nodes;
Edges=struct('nini',edges_nodes(:,1),'nfin',edges_nodes(:,2),...
    'parametric',createLine(Nodes(edges_nodes(:,1),:),...
    Nodes(edges_nodes(:,2),:)),...
    'lleft',edges_loops(:,1),'lright',edges_loops(:,2),...
    'type',char(zeros(length(edges_nodes(:,1)),1)),'order',...
    zeros(length(edges_nodes(:,1)),1));
Edges.type(:) = 'D';        %all edges are predefined as Dirichlet
Edges.order(:) = NaN;       %all degrees are predefined as NaN

Loops = struct('nodes', loops_nodes,'edges',loops_edges,...
    'center',zeros(length(loops_nodes(:,1)),2),...
    'area',zeros(length(loops_nodes(:,1)),1),...
    'order',zeros(length(loops_nodes(:,1)),1),...
    'material',zeros(length(loops_nodes(:,1)),2));

% Computation of the barycenter and area of each finite element:     
% It uses the POLYGONCENTROID function by David Legland.
for i =1:length(loops_nodes(:,1)) 
    [Loops.center(i,:),Loops.area(i)] = ...
        polygonCentroid(Nodes(loops_nodes(i,:),:));
    Loops.area(i)=abs(Loops.area(i)); 
end                                   

%% MATERIAL DATA
% Loads the material data and stores them in the Loop structure. The
% material data is assumed uniform for all elements. However, users may
% overwrite the data loaded from the GUI to define elements with distinct
% material properties. An example is given below.

% Loading the material data.
load('HeatStructDef','k','Q');

% Allocate the material characteristics defined in the GUI to all elements
Loops.material(:,1) = k;
Loops.material(:,2) = Q;
% ... or overwrite the material characteristics manually
% Loops.material(5,1) = 1;      <------- Examples
% Loops.material(5,2) = 2;      <------- Examples


%% RUN CONTROL DATA
% Loading algorithmic, refinement and edge data
load('HeatStructDef','NumberGaussPoints','EdgesOrder','LoopsOrder');
load('HeatRegBC2','edgesDirichlet','edgesNeumann','dataTemp','dataHeat');
NGP = NumberGaussPoints;   %Number of Gauss integration points per interval

%% EDGE TYPE DATA
% Registration of the Neumann edges, using edgesNeumann vector from the GUI
% It is recalled that all edges were predefined as Dirichlet.
% Users may overwrite the data loaded from the GUI to change the boundary
% types. An example is given below.
if exist('edgesNeumann')
    for i=1:length(edgesNeumann)
        Edges.type(edgesNeumann(i),1) = 'N';
    end
end

%% EDGE REFINEMENT DATA
% Allocate the refinement order defined in the GUI to all Dirichlet
% boundaries ...
Edges.order(Edges.type=='D') = EdgesOrder;
% ... or overwrite orders manually, if you need to have different orders
% for different essential boundaries
% Edges.order(:) = 4;      % <------- Examples
% Edges.order(3:4) = 0;    % <------- Examples
% Edges.order(8:3:11) = 0; % <------- Examples

%% ELEMENT DATA
% Allocate the refinement order defined in the GUI to all elements
Loops.order(:) = LoopsOrder;
% ... or overwrite orders manually, if you need to have different orders
% for different elements
% Loops.order(2:2:4) = 10; %  <------- Examples
  

%% Temperature and flux data 
% Boundary conditions can be described by polynomials of any order. The 
% definition of a boundary condition is made by specifying its values in 
% as many equally spaced points along the boundary as needed to define its
% polynomial variation. For further details regarding the definition of the 
% boundary conditions, please refer to Section 4.6 of reference [3].
%
% BConds data structure collects information regarding the boundary 
% conditions. Its members are cell arrays with as many lines as the 
% external boundaries of the structure. The values of the fluxes (or
% temperatures) enforced on the boundaries are stored in the Neumann (or
% Dirichlet) fields of the structure. NaN is stored in the Dirichlet field
% of a Neumann boundary, and vice-versa.

% Initialization of the BConds structure
BConds = struct('Neumann',{cell(length(edges_nodes),1)},'Dirichlet',...
    {cell(length(edges_nodes),1)});
BConds.Neumann(:) = {NaN};
BConds.Dirichlet(:) = {NaN};

% Dirichlet boundary conditions are imported from the GUI and stored in the
% Dirichlet field of the structure, in the normal and tangential
% directions. Users may overwrite the data loaded from the GUI to change 
% the boundary conditions. 
if exist('dataTemp')
    for i=1:size(dataTemp,1)
        BConds.Dirichlet{dataTemp{i,1},1}=str2num(dataTemp{i,2});
    end
end

% Neumann boundary conditions are imported from the GUI and stored in the
% Neumann field of the BConds structure, in the normal and tangential
% directions. Users may overwrite the data loaded from the GUI to change 
% the boundary conditions. 
if exist('dataHeat')
    for i=1:size(dataHeat,1)
        BConds.Neumann{dataHeat{i,1},1}=str2num(dataHeat{i,2});
    end
end

% NoDiv is the number of points for plotting the colormaps of the
% solution, in each Cartesian direction. It is predefined here as the
% number of Gauss points used for the side integration, but it may be
% arbitrarily defined to some other value.
NoDiv = NGP;

%% ******************* Start Checks *******************************
% The consistency of the boundary definition is checked here.
for ii=1:1:length(Edges.type)
    if (Edges.type(ii)=='D' && any(~isnan(BConds.Neumann{ii})))
        % if both of these conditions are met..
        % if Dirichlet boundary and any of the Neumann BC is not NaN
        error('local:consistencyChk',...
            'Edge %d is defined as Dirichlet, but has Neumann boundary conditions. \n',...
            ii);
    elseif (Edges.type(ii)=='D' && isnan(Edges.order(ii)))
        % if Neumann boundary and any of the Dirichlet BC is not NaN
        error('local:consistencyChk',...
            'Edge %d is defined as Dirichlet, but has no order of approximation. \n',...
            ii);
    elseif (Edges.type(ii)=='N' && any(~isnan(BConds.Dirichlet{ii})))
        % if Neumann boundary and any of the Dirichlet BC is not NaN
        error('local:consistencyChk',...
            'Edge %d is defined as Neumann, but has Dirichlet boundary conditions. \n',...
            ii);
    elseif (Edges.type(ii)=='N' && ~isnan(Edges.order(ii)))
        % if Neumann boundary and the order is not NaN
        error('local:consistencyChk',...
            'Edge %d is defined as Neumann, but has non-zero order of approximation. \n',...
            ii);
    end
end

end
