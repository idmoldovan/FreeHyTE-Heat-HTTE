function [NGP, Nodes, Edges, Loops, BConds] = InputProcReg
%--------------------------------------------------------------------------
% This is the input processing routine for the regular mesh, heat conduction
% problems. It reads the input data as defined in the GUI. It should be
% fully automatic. Only change stuff if you want to overwite manually the
% definitions from the GUI.

% Loads mesh information from previous interfaces 
load('HeatRegBC1','nodes','edges_nodes','edges_loops','loops_nodes',...
     'loops_edges');

%Definition of the data structures
Nodes = nodes;
Edges=struct('nini',edges_nodes(:,1),'nfin',edges_nodes(:,2),...
    'parametric',createLine(Nodes(edges_nodes(:,1),:),...
    Nodes(edges_nodes(:,2),:)),...
    'lleft',edges_loops(:,1),'lright',edges_loops(:,2),...
    'type',char(zeros(length(edges_nodes(:,1)),1)),'order',...
    zeros(length(edges_nodes(:,1)),1));
Edges.type(:) = 'D';        %all edges are predefined as Dirichlet
Edges.order(:) = NaN;       %all degrees are predefined as NaN

Loops = struct('nodes', loops_nodes, 'edges',loops_edges,'center',...
    zeros(length(loops_nodes(:,1)),2),'area',...
    zeros(length(loops_nodes(:,1)),1),'order',...
    zeros(length(loops_nodes(:,1)),1),'material',...
    zeros(length(loops_nodes(:,1)),2));

for i =1:length(loops_nodes(:,1)) % computing the area and centroid of each element
    [Loops.center(i,:),Loops.area(i)] = ...
        polygonCentroid(Nodes(loops_nodes(i,:),:));
    Loops.area(i)=abs(Loops.area(i)); % for the area to be correct, the nodes
    % must be listed in clockwise or counter-clockwise order!
end                                   

%% -------------------------------------------------------------------
% -------------------------USER-DEFINED AREA--------------------------

% MATERIAL DATA

load('HeatStructDef','k','Q');

% Allocate the material characteristics defined in the GUI to all elements
Loops.material(:,1) = k;
Loops.material(:,2) = Q;
% ... or overwrite the material characteristics manually
% Loops.material(5,1) = 1;      <------- Examples
% Loops.material(5,2) = 2;      <------- Examples


% RUN CONTROL DATA
load('HeatStructDef','NumberGaussPoints','EdgesOrder','LoopsOrder');
load('HeatRegBC2','edgesDirichlet','edgesNeumann','dataTemp','dataHeat');
NGP = NumberGaussPoints;   %Number of Gauss integration points per interval

% Registration of the Neumann edges, using edgesNeumann vector from the GUI
if exist('edgesNeumann')
    for i=1:length(edgesNeumann)
        Edges.type(edgesNeumann(i),1) = 'N';
    end
end

% EDGE REFINEMENT DATA

% Allocate the refinement order defined in the GUI to all Dirichlet
% boundaries ...
Edges.order(Edges.type=='D') = EdgesOrder;
% ... or overwrite orders manually, if you need to have different orders
% for different essential boundaries
% Edges.order(:) = 4;
% Edges.order(3:4) = 0;
% Edges.order(8:3:11) = 0; % <------- Examples

% ELEMENT DATA

% Allocate the refinement order defined in the GUI to all elements
Loops.order(:) = LoopsOrder;
% ... or overwrite orders manually, if you need to have different orders
% for different elements
% Loops.order(2:2:4) = 10; %  <------- Examples
  

% LOAD & DISPLACEMENT DATA (evenly pts for polynomial interpolation)
% initialization...
BConds = struct('Neumann',{cell(length(edges_nodes),1)},'Dirichlet',...
    {cell(length(edges_nodes),1)});
BConds.Neumann(:) = {NaN};
BConds.Dirichlet(:) = {NaN};

% Dirichlet boundary conditions, as imported from the GUI
if exist('dataTemp')
    for i=1:size(dataTemp,1)
        BConds.Dirichlet{dataTemp{i,1},1}=str2num(dataTemp{i,2});
    end
end

% Neumann boundary conditions
if exist('dataHeat')
    for i=1:size(dataHeat,1)
        BConds.Neumann{dataHeat{i,1},1}=str2num(dataHeat{i,2});
    end
end

% GRID POINTS FOR THE POST-PROCESSING. OPTION FOR STORAGE IN GAUSS POINTS
NoDiv = NGP;

% --------------------- END OF USER-DEFINED AREA -------------------
% ------------------------------------------------------------------

%% ******************* Checks *******************************
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
%% ******************** End of Checks +++++++++++++++++++++++++++++++++

end
