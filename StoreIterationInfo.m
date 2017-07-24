function [Loops,Edges,List] = StoreIterationInfo(Loops,Edges,List,...
    iteration,Energy0)
% Stores a lot of information about the current iteration for
% post-processing and / or for use in the next iteration

%%
% If no loop refinement is required, it registers 'None'
if isempty(List.RefinedLoopsIt{iteration,:})
    List.RefinedLoopsIt{iteration,:}='None';
end

%%
% Stores stuff for the next iteration

List.EdgesOrderIt(:,iteration)=Edges.order(:); % stores the edge's orders
List.LoopsOrderIt(:,iteration)=Loops.order(:); % Make a record per iteration
List.GDL_It(iteration)=sum(Edges.dim(:))+sum(Loops.dim(:));

%%
% Computes the energy variation relative to the previous iteration
if iteration == 1
    List.EnergyVariationIt(iteration)=abs(List.EnergyIt(iteration)-Energy0)/abs(Energy0);
elseif iteration ~= 1
    List.EnergyVariationIt(iteration)=abs(List.EnergyIt(iteration)-List.EnergyIt(iteration-1))/abs(List.EnergyIt(iteration-1));
end

% Computes the error variation relative to the first iteration
if iteration == 1
    List.ErrorEdgeVariationIt(iteration)=1;
elseif iteration ~= 1
    List.ErrorEdgeVariationIt(iteration)=abs(List.ErrorEdgeNormIt(iteration))/abs(List.ErrorEdgeNormIt(1));
end

end