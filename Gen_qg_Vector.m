function RHS = Gen_qg_Vector(Edges, Loops, BConds, RHS, abscissa, weight)
% GEN_QG_VECTOR iterates through the elements and calls the functions 
% that generate qgamma block of the free vector (RHS) of the solving
% system. qgamma is abbreviated as qg in the code.
%
% GEN_QG_VECTOR is called by MAIN***. 
% Input:
%  data the Edges, Loops structures, the RHS vector 
% (that is, the free vector of the solving system), and the Gauss-Legendre 
% integration parameters abscissa and weights. 
% Output/Returns to MAIN*** :
% the RHS vector with the qg blocks of all elements inserted at the 
% correct positions (as determined in ASSIGNPARTS).
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
% GEN_QG_VECTOR computes the internal product between the temperature basis
% expressed in a polar (r,th) referential as,
%      U = conj(R.^abs(n) .* exp(1i*Th.*n))
% and the applied normal flux field on the Neumann boundaries.
% The boundary normal fluxes are defined in the GUI by its values in an
% arbitrary number of equally-spaced points along the boundary (see
% Section 4.6 of reference [3]). A polynomial interpolation is performed
% between these values to obtain the analytic expression of the applied
% fluxes. The degree of the polynomial is equal to the number of
% flux values, minus one.
%
% Further details on the structure of the solving system are presented in
% reference [4] (Chapters 3 and 4).

%% Sweeping the elements
for ii=1:length(Loops.area)
   % LocLoop is a structure where the features of the current
   % element which are directly useful for the calculation of the
   % qg block are stored.        
    LocLoop = struct('id',ii,'edges',Loops.edges(ii,:), 'center',...
        Loops.center(ii,:),'order',Loops.order(ii),... 
        'insert',Loops.insert(ii),'dim',Loops.dim(ii));
    
    % Computing the qgi vector of element ii. Function qg_Vector_i
    % is a local function defined below. 
    qgi = qg_Vector_i(Edges,LocLoop,BConds, abscissa, weight);
    
    % Inserting the vector in the global RHS vector
    RHS(LocLoop.insert:LocLoop.insert+LocLoop.dim-1) =...
        RHS(LocLoop.insert:LocLoop.insert+LocLoop.dim-1)+ qgi;
end 
end 

%%
function qgi = qg_Vector_i(Edges,LocLoop,BConds, abscissa, weight)
%  QG_VECTOR_I local function computes the qg vector of the element LocLoop
% The sides are mapped to a [-1,1] interval to perform the integration

%% Initialization 
qgi = zeros(LocLoop.dim,1);

% n + LocLoop.order + 1 is the current line;
n = -LocLoop.order:LocLoop.order;

% Iterating on the edges for contour integration
for jj = 1:length(LocLoop.edges) 
    
    id = LocLoop.edges(jj);  % number of the jj-th edge of the loop
    
    % qg vector is only calculated on the Neumann boundaries
    if strcmpi(Edges.type(id),'N')    
 
       % LocEdge is a local structure where the features of the current
       % edge which are directly useful for the calculation of the
       % qg block are stored.     
       LocEdge =  struct('id',id,'nini',Edges.nini(id), 'nfin',Edges.nfin(id),...
            'parametric',Edges.parametric(id,:),'lleft',Edges.lleft(id),...
            'lright',Edges.lright(id));
        
        if LocEdge.lright  % exterior Neumann sides cannot have right loops
            error('local:consistencyChk',...
                'Exterior edge %d cannot have a right element. \n',...
                LocEdge.id);
        end

        % Constructing the matrices containing the n x abscissa
        % integration grid
        [N,A] = ndgrid(n,abscissa);
        
        %% Generating the geometric data
        % The following code transforms the abscissa coordinates, expressed in
        % the [-1,1] referential, to the polar coordinates required to compute
        % the values of the basis functions. 
        
        % Computing the length of the current edge
        L = sqrt(LocEdge.parametric(3)^2 + LocEdge.parametric(4)^2);
        
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
        
        %% Computing the integrands at the integration points
        % Computing the values of the flux at the abscissas:
        % obtaining the equally spaced points on [-1,1] interval where the
        % fluxes are defined and stored in BConds.Neumann
        a = linspace(-1,1,length(BConds.Neumann{id}));
        
        % Obtaining the polynomial that gets the values in BConds.Neumann
        % at the points a
        pol = polyfit(a,BConds.Neumann{id},length(BConds.Neumann{id})-1);
        
        % computing the values of "pol" at the abscissas
        q = polyval(pol,abscissa);
        q = q.'; % non-conjugate transpose
        
        % Temperature basis in polar coordinates
        Ustar = conj(R.^(abs(N)) .* exp(1i*Th.*N));
        
        %% Computing the integral on the side
        % Computing the internal product of the temperature basis with the
        % enforced fluxes
        qgi2D = bsxfun(@times, Ustar, q);
        
        % Computes the integral
        qgi = qgi + L/2 * sum(bsxfun(@times,qgi2D,weight.'),2);
        
    end 
    
end

end 
