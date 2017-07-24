function [Edges, Loops, Dim] = AssignParts(Edges, Loops)

% The porpuse of this function is to assign to each element and Dirichlet
%side an entry point and a dimension in the global solving system
%------------------------------------------------------------------------

entry = 1;

% initialize the insertion points for loops
Loops.insert = zeros(length(Loops.area),1);


for i = 1:length(Loops.area)
    Loops.insert(i) = entry;
    Loops.dim(i) = 2*Loops.order(i)+1; % calculates the total dimension of each element's basis
    entry = entry + Loops.dim(i);
end

% initialize the insertion points for the edges
Edges.insert = zeros(length(Edges.type),1);


for i= 1:length(Edges.insert)
    if strcmpi(Edges.type(i),'D')
        Edges.insert(i)= entry;
        Edges.dim(i) = Edges.order(i)+1;
        entry = entry + Edges.dim(i);
    end
end
% note that Neumann edges correspond to zero dim, zero insertion
% points

Dim = entry-1;

end
