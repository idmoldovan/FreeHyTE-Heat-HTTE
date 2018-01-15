function [Loops,Edges] = CheckMinOrders(Loops,Edges)
% CHECKMINORDERS checks if the initial orders of the domain and Dirichlet 
% boundaries meet the positive kinematic indeterminacy requirement. The 
% kinematic indeterminacy condition is covered in Section 3.4.2 of 
% reference [2], section 4.5.2 of reference [3], and section 3.8 of 
% reference [4]. If the criterion is not satisfied the refinement order is 
% increased in the domain of the element.
%
% CheckMinDegrees is called by MAIN*** 
% Input:
%   Loops, Edges data structures
% Output/Returns to MAIN***
%   Loops, Edges data structures with updated refinement orders (if required)
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

% initialization
Loops.dim = zeros(length(Loops.area),1);
Edges.dim = zeros(length(Edges.type),1);

% for every element, calculates the total dimension of the basis
for i=1:length(Loops.area)
    Loops.dim(i) = 2*Loops.order(i)+1; 
end

% for every essential boundary, calculates the total dimension of the basis
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

% Increase the order of the element basis if Beta <= 0
for i=1:length(Beta)
    while Beta(i)<=0
        Loops.order(i)=Loops.order(i)+1;
        % recompute the dimension of the element's basis
        Loops.dim(i)= 2*Loops.order(i)+1;
        % recomputes the kinematic indeterminacy number, Beta
        Beta(i) = Loops.dim(i) - sum(Edges.dim(Loops.edges(i,:)));
    end
end
end
