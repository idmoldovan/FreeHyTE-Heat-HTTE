function LHS = Gen_B_Matrix(Edges,Loops,LHS,abscissa,weight)
% GEN_B_MATRIX iterates through the edges and calls the functions that 
% generate the B block of the boundary matrix in the LHS. 
%
% B is called by MAIN***. It receives as input data the Edges and Loops
% structures, the LHS matrix (that is, the matrix of coefficients of the
% solving system), and the Gauss-Legendre integration parameters abscissa
% and weights. It returns to MAIN*** the LHS matrix with the B blocks of
% all elements inserted at the correct positions (as determined in
% ASSIGNPARTS).
%
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
% GEN_B_MATRIX computes the internal product between the following bases, 
% expressed in a polar (r,th) referential:
% * the temprature basis U, generated by harmonic polynomials, of order n  
%      U = r^abs(n) .* exp(i*th*n);
% * the boundary flux basis Z, defined by Chebyshev polynomials,
%      Z  = cos(m*ArcCos(abscissa))
% where n and m are the line and column of the current term in matrix B.
%
% Further details on the structure of the solving system are presented in 
% reference [2] (Section 6.2). The transformation of referentials is
% covered in reference [4] (Appendix C). 
%

%% Start iteration on the Dirichlet edges 
for ii=1:length(Edges.type)
    
    % Boundary blocks are constructed for "Dirichlet" (and interior)
    % boundaries only.
    if strcmpi(Edges.type(ii),'D')
       % LocEdge is a local structure where the features of the current
       % edge which are directly useful for the calculation of the
       % boundary block are stored.
       LocEdge = struct('nini',Edges.nini(ii),'nfin',Edges.nfin(ii),...
            'parametric',Edges.parametric(ii,:),'lleft',...
            Edges.lleft(ii),'lright',Edges.lright(ii),...
            'order',Edges.order(ii),'insert',Edges.insert(ii),...
            'dim',Edges.dim(ii));
        
        % Generating the boundary block corresponding to the current edge
        % and its LEFT finite element
        if LocEdge.lleft 
            id = LocEdge.lleft;
            sign = 1.;
            % LocLoop is a local structure where the features of the left
            % element which are directly useful for the calculation of the
            % boundary block are stored.
            LocLoop = struct('id',id,'edges',Loops.edges(id,:),...
                'center',Loops.center(id,:),'order',...
                Loops.order(id),'insert',Loops.insert(id),...
                'dim',Loops.dim(id));
            
            % Launching the function that computes the boundary
            % block B for the current edge ii and its left element
            Bdi = sign*B_Matrix_i(LocEdge,LocLoop,abscissa,weight);
            
            % Inserting the matrix in the global LHS matrix
            LHS(LocLoop.insert:LocLoop.insert+LocLoop.dim-1,...
                LocEdge.insert:LocEdge.insert+LocEdge.dim-1) = -Bdi;
            % Inserting the conjugate transposed in the global LHS matrix
            LHS(LocEdge.insert:LocEdge.insert+LocEdge.dim-1,...
                LocLoop.insert:LocLoop.insert+LocLoop.dim-1) = -Bdi';
            
        else % there should always be a left element
            error('local:consistencyChk',...
                'No left loop for edge %d. \n', ii);
        end
        
        % Generating the boundary block corresponding to the current edge
        % and its RIGHT finite element. 
        % NOTE: right elements exist for internal boundries only
      
        if LocEdge.lright
            id = LocEdge.lright;
            sign = -1.;
            % LocLoop is a local structure where the features of the right
            % element which are directly useful for the calculation of the
            % boundary block are stored.
             LocLoop = struct('id',id,'edges',Loops.edges(id,:),...
                'center',Loops.center(id,:),'order',...
                Loops.order(id,1),'insert',Loops.insert(id,1),...
                'dim',Loops.dim(id,1));
            
            % Computing the B matrix of the edge ii, right loop
            Bdi = sign*B_Matrix_i(LocEdge,LocLoop,abscissa,weight);
            
            % Inserting the matrix in the global LHS matrix
            LHS(LocLoop.insert:LocLoop.insert+LocLoop.dim-1,...
                LocEdge.insert:LocEdge.insert+LocEdge.dim-1) = -Bdi;
            % Inserting the conjugate transposed in the global LHS matrix
            LHS(LocEdge.insert:LocEdge.insert+LocEdge.dim-1,...
                LocLoop.insert:LocLoop.insert+LocLoop.dim-1) = -Bdi';
            
        end
    end
end
end


%%
function Bdi = B_Matrix_i(LocEdge,LocLoop,abscissa,weight)
% B_MATRIX_I is a local function that computes B boundary block for the LocEdge
% LocEdge edge and LocLoop element. The side is mapped to a [-1,1] interval 
% to perform the integration.

%% Initialization
% n + LocLoop.order + 1 is the current line; 
% m + 1 is the current column
n = -LocLoop.order:LocLoop.order;
m = 0:LocEdge.order;
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
Bdi = (L/2)*sum(bsxfun(@times,W3D,Bdi3D),3);% computes the integral

end

