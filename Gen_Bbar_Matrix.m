function Abar = Gen_Bbar_Matrix(LHS0,Edges,Loops,Dim,abscissa,weight,index)
% GEN_BBAR_MATRIX computes the BBar matrix/vector for the boundary to be 
% refined and stores it in the augmented matrix ABar. The BBar matrix
% is actually a vector if a single shape function is added to the basis.
% It is defined by expression (5.3) of reference [4].
% GEN_BBAR_MATRIX is called by EDGEREFINEMENT. 
%
% Input:
% the Edges, Loops structures, the matrix of coefficients before the 
% refinement (LSH0), the dimension of the system (Dim), the id of the 
% refined boundary, and the Gauss-Legendre integration parameters 
% abscissa and weights. 
% Output/Returns to EDGEREFINEMENT:
% the Bbar matrix for the refined edge 
% Calls:  BBAR_MATRIX_I to compute BBar matrix
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

% Create the Abar matrix and intialize with zero. The refinement of the 
% edge adds one degree of freedom to the boundary, so the dimension of the 
% Abar matrix is 'Dim+1'. 
Abar=zeros(Dim+1,Dim+1);    

% Insert the "old" A matrix into the new Abar as only the new terms 
% associated with the new dof need to be computed.
Abar(1:Dim,1:Dim)=LHS0;			              

% Identify the refined edge
ii=index;

% Update edges order and dimension of the B matrix
Edges.order(ii)=Edges.order(ii)+1; 
Edges.dim(ii)=Edges.dim(ii)+1;

if strcmpi(Edges.type(ii),'D')
    % LocEdge is a structure where the features of the current
    % edge which are directly useful for the calculation of the
    % Bbar are stored.
    LocEdge = struct('nini',Edges.nini(ii),'nfin',Edges.nfin(ii),...
        'parametric',Edges.parametric(ii,:),'lleft',...
        Edges.lleft(ii),'lright',Edges.lright(ii),...
        'order',Edges.order(ii),'insert',Edges.insert(ii),...
        'dim',Edges.dim(ii));
    
    if LocEdge.lleft  % The edge should always have a left element
        id = LocEdge.lleft;
        sign = 1.;
        
        % LocLoop is a structure where the features of the current element
        % which are directly useful for the calculation of the Bbar are stored.
        LocLoop = struct('id',id,'edges',Loops.edges(id,:),...
            'center',Loops.center(id,:),'order',...
            Loops.order(id),'insert',Loops.insert(id),...
            'dim',Loops.dim(id));
        
        % Computing the B matrix of the edge ii, left loop. Function
        % BBAR_MATRIX_I is a local function, defined below
        Bbardi = sign*Bbar_Matrix_i(LocEdge,LocLoop,abscissa,weight);
        
        % Inserting the new B matrix in the global Abar matrix
        Abar(LocLoop.insert:LocLoop.insert+LocLoop.dim-1,Dim+1) = Bbardi;
        % Inserting the conjugate transposed in the global LHS matrix
        Abar(Dim+1,LocLoop.insert:LocLoop.insert+LocLoop.dim-1) = Bbardi';
        
    end % End left edge
    
    if LocEdge.lright  % Right elements may or may not exist
        id = LocEdge.lright;
        sign = -1.;
        
        % LocLoop is a structure where the features of the current element
        % which are directly useful for the calculation of the Bbar are stored.         
        LocLoop = struct('id',id,'edges',Loops.edges(id,:),...
        'center',Loops.center(id,:),'order',...
            Loops.order(id,1),'insert',Loops.insert(id,1),...
            'dim',Loops.dim(id,1));
        
        % Computing the Bbar matrix of the edge ii, right loop. Function
        % BBAR_MATRIX_I is a local function defined below
        Bbardi = sign*Bbar_Matrix_i(LocEdge,LocLoop,abscissa,weight);
        
        % Inserting the matrix in the global LHS matrix
        Abar(LocLoop.insert:LocLoop.insert+LocLoop.dim-1,Dim+1) = Bbardi;
        % Inserting the conjugate transposed in the global LHS matrix
        Abar(Dim+1,LocLoop.insert:LocLoop.insert+LocLoop.dim-1) = Bbardi';
        
    end % End right loop
end

% Reset edge's order and dimension of the B matrix for current edge
Edges.order(ii)=Edges.order(ii)-1;
Edges.dim(ii)=Edges.dim(ii)-1;

end


%%
function Bbardi = Bbar_Matrix_i(LocEdge,LocLoop,abscissa,weight)
% BBAR_MATRIX_I  local function computes the Bbardi matrix of the edge ii
% The sides are mapped to a [-1,1] interval to perform the integration

%% Initialization
% n + LocLoop.order + 1 is the current line; 
% there is a single column if the basis' order is incremented in 1
n = -LocLoop.order:LocLoop.order;
m = LocEdge.order;                    
a = abscissa;

%% Generating the geometric data
% The following code transforms the abscissa coordinates, expressed in
% the [-1,1] referential, to the polar coordinates required to compute
% the values of the basis functions. The components of the outward
% normal to the boundary in the radial and tangential directions are
% also calculated. 

% Computing the length of the current edge
L = sqrt((LocEdge.parametric(3))^2 + (LocEdge.parametric(4))^2);

% Constructing the 3D matrices containing the n x m x abscissa
% integration grid
[N,M,A]=ndgrid(n,m,a);

%% Getting the r, th, nr, nth for all Gauss points
% Transforming the edge abscissa into local coordinates. The local
% referential is centered in the barycenter of the element, its axes
% aligned with the Cartesian axes of the global referential.
loc_x = LocEdge.parametric(1) - LocLoop.center(1) + 0.5 *...
    (A + 1) * LocEdge.parametric(3); 
loc_y = LocEdge.parametric(2) - LocLoop.center(2) + 0.5 *...
    (A + 1) * LocEdge.parametric(4);

% Transforming the local Cartesian coordinates into polar.
R = sqrt(loc_x.^2 + loc_y.^2);  
Th = atan2(loc_y, loc_x);

%% Computing the basis functions for all integration points
% Temperature basis in polar coordinates
Ustar = conj(R.^abs(N) .* exp(1i*Th.*N));

%% Computing the boundary integral
% The integral is the internal product between the temperature basis U*
% and the flux basis Z (Chebyshev polynomials)
Bdi3D = Ustar .*...
    cos((M).*acos(A)); % Chebyshev functions

% Performing the side integration 
W3D(1,1,:) = weight;
Bbardi = (L/2)*sum(bsxfun(@times,W3D,Bdi3D),3); % computes the integral

end
