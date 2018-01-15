function List = SelectLoopsToRefine(Loops,Edges,List, iteration)
% SELECTLOOPSTOREFINE identifies the elements (loops) to be refined, namely
% 1. loops with boundaries selected for refinement that would cause
% spurious modes;
% 2. loops with non-positive kinematic indeterminacy numbers; and 
% 3. loops with boundaries of larger refinement than the loops themselves.
%
% SELECTLOOPSTOREFINE is called by MAIN*** 
% Input data:
%  the Loops, Egdes, and List strucutres, and the current iteration
%  counter.
% Output/Returns to MAIN***:
%  the updated List structure with the loops to be refined 
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
% The domain refinement is always triggered by boundary refinements, in
% order to avoid the over-constrainment of the domain bases and the
% instability of the solving system that could be caused by the addition of
% the boundary modes. 
% To prevent over-constrainment of the domain bases, the finite elements
% are refined whenever the addition of boundary modes would render them
% statically or kinematically (over-)determinate. 
% To avoid the local over-constrainment of the domain bases, their orders 
% are maintained larger than the orders of all neighbouring essential
% boundaries. 
% To handle the instability of the solving system, whenever the refinement
% of a boundary basis causes the emergence of singular value outliers
% new modes are added to the finite elements adjacent to that boundary
% until the outliers vanish. 
%
% Further details on the refinement selection process can be found in
% reference [2] (Section 6.3.2), and reference [4] (Sections 5.2 to 5.4,
% and 6.3). 


%% Refinement criterion # 1
% Identify loops adjacent to the boundaries with spurious modes selected 
% for refinement
if ~isempty(List.SpurEdgeToRefine)
    for i=1:length(List.SpurEdgeToRefine)
        index = List.SpurEdgeToRefine(i);
        List.LoopsToRefine = cat(2,List.LoopsToRefine,...
            Edges.lleft(index),Edges.lright(index));
    end
end

%% Refinement criterion # 2
% Identify loops with negative or null kinematic indeterminacy number
for i=1:length(Loops.area)
    
    % the kinematic indeterminacy number is equal to the difference between
    % the number of functions present in the domain basis of an element and
    % the total number of functions present in the bases of its essential
    % boundaries.
    List.BetaIt(iteration,i)=Loops.dim(i)-(sum(Edges.dim(Loops.edges(i,:))));
    
    % if the Beta of the current loop is not positive, it is listed for
    % refinement
    if List.BetaIt(iteration,i)<=0
        List.LoopsToRefine = cat(2,List.LoopsToRefine,i);
    end
end

%% Refinement criterion # 3
% Identify loops with adjacent essential boundaries of order larger
% than the order of the loop
for i=1:length(Loops.area)
    % compares the order of the loop with the maximum order of the
    % adjacent boundaries. It should omit (by default) the NaN values 
    % associated to Neumann boundaries.
    if Loops.order(i) <= max(Edges.order(Loops.edges(i,:)))
        List.LoopsToRefine = cat(2,List.LoopsToRefine,i);
    end
    
end

%% Constructing the final list with loops to be refined
if ~isempty(List.LoopsToRefine)
    % eliminate the duplicates
    List.LoopsToRefine = unique(List.LoopsToRefine); 
    % eliminate the zero entries, caused by the exterior edges
    List.LoopsToRefine = List.LoopsToRefine(List.LoopsToRefine~=0);
    List.RefinedLoopsIt{iteration,:}=List.LoopsToRefine;
end


end
