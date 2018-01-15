function List = StoreIterationInfo(Loops,Edges,List,...
    iteration,Energy0)
% STOREITERATIONINFO saves relevant information on the current iteration
% for use mainly in post-processing
%
% StoreIterationInfo is called by MAIN*** 
% 
% Input:
%  Loops, Edges, and List structures, current iteration counter and 
%  initial solution energy Energy0
% Output/Returns to MAIN***:
%  List structure updated with: 
%  * List.EdgesOrderIt - stores the orders of the edges at every iteration;
%  * List.LoopsOrderIt - stores the orders of the loops at every iteration; 
%  * List.GDL_It - stores the total no of DOFs at every iteration;
%  * List.EnergyVariationIt - solution's energy variation, relative to the
%  previous iteration;
%  * List.ErrorEdgeVariationIt - boundary temperature balance error
%  relative to the first iteration.
% Calls:
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

%% RefinedLoopsIt
% If no loop refinement is required, it registers 'None'
if isempty(List.RefinedLoopsIt{iteration,:})
    List.RefinedLoopsIt{iteration,:}='None';
end

%% EdgesOrderIt, LoopsOrderIt, GDL_It
% Stores the relevant information for the next iteration
List.EdgesOrderIt(:,iteration)=Edges.order(:); % stores the edge's orders
List.LoopsOrderIt(:,iteration)=Loops.order(:); % Make a record per iteration
List.GDL_It(iteration)=sum(Edges.dim(:))+sum(Loops.dim(:));

%% EnergyVariationIt
% Computes the energy variation relative to the previous iteration
if iteration == 1
    List.EnergyVariationIt(iteration)=...
        abs(List.EnergyIt(iteration)-Energy0)/abs(Energy0);
elseif iteration ~= 1
    List.EnergyVariationIt(iteration)=abs(List.EnergyIt(iteration)-...
        List.EnergyIt(iteration-1))/abs(List.EnergyIt(iteration-1));
end

%% ErrorEdgeVariationIt
% Computes the error variation relative to the first iteration
if iteration == 1
    List.ErrorEdgeVariationIt(iteration)=1;
elseif iteration ~= 1
    List.ErrorEdgeVariationIt(iteration)=...
        abs(List.ErrorEdgeNormIt(iteration))/abs(List.ErrorEdgeNormIt(1));
end

end
