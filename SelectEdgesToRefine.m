function [List,ExitScore] = SelectEdgesToRefine(List,SelectionCriterion,SelectionTol,thresh)
% Based on the ranking of the selection criteria on the essential
% boundaries, it selects the boundaries to be refined and identifies the
% spurious ones

ExitScore = 1;

%%
% The first-ranked element in the list is selected for refinement
List.EdgesToRefine = List.Edge(1,1);

%%
% Sweeps all other edges and checks if the values of their selection
% criteria stand within the defined selection tolerance from the
% first-ranked edge
for i=2:length(List.Edge)
    if List.Edge(i,SelectionCriterion) >= ...
            List.Edge(1,SelectionCriterion)*SelectionTol
        List.EdgesToRefine = cat(2,List.EdgesToRefine,List.Edge(i,1));
    end
end

% Clear rows of Neumann Edges in List.Edge
List.Edge(~any(List.Edge,2),:)=[];

%%
% Analyses the values of the SelectionCriterion and if they are
% smaller than a threshold for all boundaries, it means that the exact
% solution was perfectly recovered at the previous step. If this is the
% case, it will plot the solution and exit
if max(List.Edge(:,SelectionCriterion)) < thresh
    warning('local:NumericalChk',...
        'The values of the selection criterion are smaller than the threshold %g for all boundaries. Suspect that the exact solution was already achieved.\n',thresh);
    ExitScore = 0;
end

%%
% Analyses the ErrorEdgeNorm (column 2 in the List.Edge)  and if it's
% smaller than a threshold, registers it as a spurious mode
for index=1:size(List.Edge,1)
    if List.Edge(index,2)<thresh
        List.SpuriousEdges = cat(1,List.SpuriousEdges,List.Edge(index,1));
    end
end


end