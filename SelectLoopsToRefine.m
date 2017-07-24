function List = SelectLoopsToRefine(Loops,Edges,List, iteration)
% Identifies the elements to be refined, namely loops with boundaries
% selected for refinement that would cause spurious modes, loops with
% non-positive kinematic indeterminacy numbers and loops with boundaries of
% larger refinement than the loops themselves

%%
% Identify loops adjacent to the boundaries with spurious modes selected for refinement
if ~isempty(List.SpurEdgeToRefine)
    for i=1:length(List.SpurEdgeToRefine)
        index = List.SpurEdgeToRefine(i);
        List.LoopsToRefine = cat(2,List.LoopsToRefine,Edges.lleft(index),Edges.lright(index));
    end
end

%%
% Identify loops with negative or null kinematic indeterminacy
for i=1:length(Loops.area)
    List.BetaIt(iteration,i)=Loops.dim(i)-(sum(Edges.dim(Loops.edges(i,:))));
    
    % if the Beta of the current loop is not positive, it is listed for
    % refinement
    if List.BetaIt(iteration,i)<=0
        List.LoopsToRefine = cat(2,List.LoopsToRefine,i);
    end
end

%%
% Identify loops with adjacent essential boundaries of order larger
% than the order of the loop
for i=1:length(Loops.area)
    % compares the order of the loop with the maximum order of the
    % adjacent boundaries. It omits the NaN values associated to
    % Neumann boundaries.
    if Loops.order(i) <= nanmax(Edges.order(Loops.edges(i,:)))
        List.LoopsToRefine = cat(2,List.LoopsToRefine,i);
    end
    
end

%%
% Constructing the final list with loops to be refined
if ~isempty(List.LoopsToRefine)
    List.LoopsToRefine = unique(List.LoopsToRefine); % eliminate the duplicates
    % eliminate the zero entries, caused by the exterior edges
    List.LoopsToRefine = List.LoopsToRefine(List.LoopsToRefine~=0);
    List.RefinedLoopsIt{iteration,:}=List.LoopsToRefine;
end


end