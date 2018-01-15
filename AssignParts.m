function [Edges, Loops, Dim] = AssignParts(Edges, Loops)
% ASSIGNPARTS maps the solving system and assigns each (conductivity or
% boundary) block an entry point and a dimension.
%
%
% ASSIGNPARTS is called by MAINREG and MAINTRI. 
% Input: 
%  receives structures Edges, Loops and BConds as input arguments, 
% Return: 
%  Edges and Loops structures, updated with two fields containing the 
%  insertion points of the respective blocks in the solving system, 
%  and their dimensions. It also returns Dim, which is the dimension 
%  of the solving system.
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
% The layout of a solving system with a single finite element a single
% essential boundary is presented below.
%
%
%   _________________
%   |       |       |
%   |   K   |   B   |
%   |_______|_______|
%   |    t  |       |
%   |   B   |   0   |
%   |_______|_______|
%
%
% The conductivity matrix K belonging to element i is inserted on the main 
% diagonal of the solving system at position Loops.insert(i). 
% Its dimension is equal to the order of the domain basis of the respective 
% element, Loops.order(i) ( 2*order of the approximation basis +1). 
%
% Each boundary 'j' has B block. Interior boundaries always have approximations.
% Exterior Dirichlet boundaries have approximations where temperatures
% are enforced. For more insight on how to model different types of 
% boundary conditions, please consult Section 4.6 of reference [3].
% The Bj blocks are inserted on the lines that correspond to the entries of
% the conductivity block of the element(s) the boundary belongs to. 
% The entry line and column of the B block corresponding to boundary j and
% neighbouring element i are Loops.insert(i) and Edges.insert(j).
% The Bj blocks are not square. They have as many lines as the shape
% functions of the neighbouring elements and Edges.dim(i) = Edges.order(i)+1 
% lines.
%
% The layout and storage of the solving system is discussed at length in
% references [2] (Section 6.2) and [4] (Section 4.4).
% 

%% Initialization
% Initializes the current line indicator, entry.
entry = 1;

% Initializes the insertion points for the K blocks (loops)
Loops.insert = zeros(length(Loops.area),1);
% Initializes the insertion points for the boundary blocks
Edges.insert = zeros(length(Edges.type),1);

%% Mapping of the conductivity entry blocks of the system
% The insertion point of the K block and its dismension are
% computed for each element and stored in the Loops structure.
for i = 1:length(Loops.area)
    Loops.insert(i) = entry;
    % calculates the total dimension of each element's basis
    Loops.dim(i) = 2*Loops.order(i)+1; 
    entry = entry + Loops.dim(i);
end

%% Mapping of the boundary entry blocks of the system
% The insertion point of the boundary block and its dimension are
% computed for each edge and stored in the Edges structure.
for i= 1:length(Edges.insert)
    if strcmpi(Edges.type(i),'D')
        Edges.insert(i)= entry;
        Edges.dim(i) = Edges.order(i)+1;
        entry = entry + Edges.dim(i);
    end
end
% note that Neumann edges correspond to zero dim, zero insertion
% points

% Computing the total dimension of the solving system
Dim = entry-1;

end
