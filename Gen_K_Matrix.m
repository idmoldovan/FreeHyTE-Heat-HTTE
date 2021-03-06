function LHS = Gen_K_Matrix(Edges, Loops, LHS, abscissa, weight)
% GEN_K_MATRIX sweeps through the elements and calls the functions 
% that generate the  K block matrix in the LHS
%
% GEN_K_MATRIX is called by MAIN***. 
% Input data: 
% the Edges and Loops structures, the LHS matrix (that is, the matrix of 
% coefficients of the solving system), and the Gauss-Legendre integration 
% parameters abscissa and weights. 
% Output/returns to MAIN***:
% the LHS matrix with the K blocks (conductivity matrix) of all elements 
% inserted at the  correct positions (as determined in ASSIGNPARTS).
%
% BIBLIOGRAPHY
% 1. FreeHyTE Page - https://sites.google.com/site/ionutdmoldovan/freehyte
% 2. Moldovan ID, Cismasiu I - FreeHyTE: theoretical bases and developer�s 
% manual, https://drive.google.com/file/d/0BxuR3pKS2hNHTzB2N2Q4cXZKcGc/view
% 3. FreeHyTE Heat HTTE User's Manual - 
%    https://drive.google.com/drive/folders/0BxuR3pKS2hNHaFhiSjZHOE9TMzg
% 4. Geraldes, MRCB,  Elementos finitos h�bridos-Trefftz adaptativos para
%   problemas de condu��o de calor. MSc Thesis, Universidade Nova de Lisboa,
%   2016 (in Portuguese).
%
%
% GEN_K_MATRIX computes the internal product between the following bases
% expressed in a polar (r,th) referential:
% * the temperature basis U generated by harmonic polynomials, of order n 
%      U = R.^abs(n) .* exp(1i*Th.*n);
% * the boundary normal flux basis N * S, where S is the flux basis is, 
%         S = | Sr | = -k * r^(abs(m)-1)* exp(i*Th.*m) | abs(m)|
%             | St |                                   | i*m   |
% and N is the directory cosine vector,
%         N = | nr  nt |
% In the above expressions, n and m are the line and column of the current 
% term in matrix K, k is the themal conductivity coefficient (see 
% INPUTPROC***), and nr and nt are the radial and tangential components of
% the outward unit normal to the current boundary.
%
% As typical of the Trefftz method, the temperature and flux bases solve
% exactly the governing (Poisson) equation in the domain of each element. 
% As a direct consequence, the stiffness matrix can be calculated using 
% boundary integrals only.
%
% The derivation of the bases from the solution of the Poisson
% equation is presented in reference [4] (Section 3.3).

%% Sweeping the elements
for ii=1:length(Loops.area)
    % LocLoop is a structure where the features of the current
    % element which are directly useful for the calculation of the
    % conductivity block are stored.    
    LocLoop = struct('id',ii,'edges',Loops.edges(ii,:), 'center',...
        Loops.center(ii,:),'order',Loops.order(ii,:),...
        'insert',Loops.insert(ii,:),'dim',Loops.dim(ii,:),...
        'k',Loops.material(ii,1));
    
    % Computing the Ki matrix of element ii. Function K_MATRIX_I is a
    % local function (see below).
    Ki = K_Matrix_i(LocLoop, Edges, abscissa, weight);
    
    % Inserting the matrix in the global LHS matrix. The insertion is made
    % at line & column Loops.insert(ii).
    LHS(LocLoop.insert:LocLoop.insert+LocLoop.dim-1,...
        LocLoop.insert:LocLoop.insert+LocLoop.dim-1) = Ki;
    
end

end


function Ki = K_Matrix_i(LocLoop, Edges, abscissa, weight)
% K_MATRIX_I is a local function that computes the K conductivity 
% block of the LocLoop element. 
% The sides are mapped to a [-1,1] interval to perform the integrations.

%% Initialization 
% Initialization of the Ki block
Ki = zeros(LocLoop.dim);
k = LocLoop.k;

% n + LocLoop.order + 1 is the current line; 
% m + LocLoop.order + 1 is the current column.
n = -LocLoop.order:LocLoop.order;
m = -LocLoop.order:LocLoop.order;

% Iterating on the edges for contour integration
for jj = 1:length(LocLoop.edges) 
    
  % number of the jj-th edge of the loop
  id = LocLoop.edges(jj);  
     
    % LocEdge is a local structure where the features of the current
    % edge which are directly useful for the calculation of the
    % conductivity block are stored.    
    LocEdge =  struct('id',id,'nini',Edges.nini(id), 'nfin',Edges.nfin(id),...
        'parametric',Edges.parametric(id,:),'lleft',Edges.lleft(id),...
        'lright',Edges.lright(id));
    
   
    %% Generating the geometric data
    % The following code transforms the abscissa coordinates, expressed in
    % the [-1,1] referential, to the polar coordinates required to compute
    % the values of the basis functions. The components of the outward
    % normal to the boundary in the radial and tangential directions are
    % also calculated. 
    
    % Computing the length of the current edge
    L = sqrt(LocEdge.parametric(3)^2 + LocEdge.parametric(4)^2); 
 
   
    % Constructing the 3D matrices containing the n x m x abscissa
    % integration grid
    [N,M,A] = ndgrid(n,m,abscissa);
    
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
    
    % Computing the components of the outward normal in the Cartesian
    % directions.
    nx = LocEdge.parametric(4) / L;   
    ny = -1* LocEdge.parametric(3) / L;
    if LocEdge.lright==LocLoop.id  % if the element is on the right,
        nx = -nx;                  % change the sign of the normal
        ny = -ny;
    end
    
    % Computing the components of the outward normal in the polar
    % directions.
    NR = nx * cos(Th) + ny * sin(Th);
    NTh = -1*nx * sin(Th) + ny * cos(Th);
    
    %% Computing the basis functions for all integration points
    % Temperature basis in polar coordinates
    Ustar = conj(R.^abs(N) .* exp(1i*Th.*N));
    
    % Polar components of the flux basis S
    SR = -k * abs(M) .* R.^(abs(M)-1) .* exp(1i*Th.*M);
    STh = -k * 1i*M.*R.^(abs(M)-1) .* exp(1i*Th.*M);
    
    NS = NR.*SR + NTh.*STh;
    
    %% Computing the integral on the side
    % The integral is the internal product between the temperature basis U*
    % and the normal flux basis N*S
    
    NS = NR.*SR + NTh.*STh;
    Ki3D = Ustar .* NS; 
    
    % Performing the side integration and updating the K matrix
    w3D(1,1,:) = weight;
    Ki = Ki + L/2 * sum(bsxfun(@times,Ki3D,w3D),3); 
    
end

end
