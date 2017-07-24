function [Loops,Edges] = CheckMinDegrees(Loops,Edges)
% This function checks is where the minimun degrees are checked to see if they meet the requirement
% This is done by comparing the chosen degrees, at an elementary level, of the Dirichlet edges and of the
% domain basis.

Loops.dim = zeros(length(Loops.area),1);
for i=1:length(Loops.area)
    Loops.dim(i) = 2*Loops.order(i)+1; % calculates the total dimension of each element's basis
end
Edges.dim = zeros(length(Edges.type),1);
for i=1:length(Edges.type)
    if strcmpi(Edges.type(i),'D')
        Edges.dim(i) = Edges.order(i)+1;
    end
end

% Computes the indeterminacy number Beta for each element
Beta = zeros(length(Loops.area),1);
for i=1:length(Loops.area)
    Beta(i) = Loops.dim(i) - sum(Edges.dim(Loops.edges(i,:)));
end

for i=1:length(Beta)
    while Beta(i)<=0
        Loops.order(i)=Loops.order(i)+1;
        Loops.dim(i)= 2*Loops.order(i)+1;
        Beta(i) = Loops.dim(i) - sum(Edges.dim(Loops.edges(i,:)));
    end
end
end
