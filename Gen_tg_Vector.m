function RHS = Gen_tg_Vector(Edges, BConds, RHS, abscissa, weight)
% GEN_TG_VECTOR Iiterates through the edges and calls the functions that 
% generate tgamma block of the free vector (RHS) of the solving system.
% (The tgamma vector is abbreviated in the code as tg). The tg "blocks"
% are only computed on the exterior Dirichlet boundaries. For all
% other elements, the tg blocks are filled with zeros.
%
% GEN_TG_VECTOR is called by MAIN***.
% Input:
%  data the Edges, Loops and BConds structures, the RHS vector
% (that is, the free vector of the solving system), and the Gauss-Legendre
% integration parameters abscissa and weights.
% Output/Returns to MAIN*** :
% the RHS vector with the tg blocks of all elements inserted at the
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
% GEN_TG_VECTOR computes the internal product between the boundary flux
% basis Z on the Dirichlet edge and the applied temperatures on the same edge. 
% * the boundary flux basis Z is defined by Chebyshev polynomials,
%        Z  = cos(n*ArcCos(abscissa))
% where n is the line of the current entry in vector tg. 
% * the boundary temperature function is defined in the GUI by its values 
% in an arbitrary number of equally-spaced points along the boundary (see
% Section 4.6 of reference [3]). A polynomial interpolation is performed
% between these values to obtain the analytic expression of the applied
% temperature. The degree of the polynomial is equal to the number of
% temperature values, minus one.
%
% Further details on the structure of the solving system are presented in
% reference [4] (Chapters 3 and 4).

%% Sweeping the edges and selecting the exterior Dirichlet boundaries
for ii=1:length(Edges.type)
    
    % Exterior Dirichlet boundaries have no right element
    if (strcmpi(Edges.type(ii),'D') && Edges.lright(ii) == 0 )
        
        % LocEdge is a local structure where the features of the current
        % edge which are directly useful for the calculation of the
        % tg block are stored.
        LocEdge  = struct('id',ii,'nini',Edges.nini(ii),...
            'nfin',Edges.nfin(ii),'parametric',Edges.parametric(ii,:),...
            'order',Edges.order(ii),'insert',Edges.insert(ii),...
            'dim',Edges.dim(ii));
        
        % Computing the tg vector of edge ii
        tgi = tg_Vector_i(LocEdge,BConds, abscissa, weight);
        
        % Inserting the vector in the global RHS vector
        RHS(LocEdge.insert:LocEdge.insert+LocEdge.dim-1) = -tgi;
    end 
    
end 
end 

%%
function tgi = tg_Vector_i(LocEdge, BConds, abscissa, weight)
% TG_VECTOR_I local function computes the tg vector of the LocEdge exterior
% Dirichlet boundary. The edge is mapped to a [-1,1] interval to perform 
% the integration.

% n+1 denotes the current line of the block
n = 0:LocEdge.order;

%% Generating the geometric data
% Computing the length of the current edge
L = sqrt(LocEdge.parametric(3)^2 + LocEdge.parametric(4)^2);

%% Computing the integrands at the integration points
% Computing the values of the normal flux basis
Zstar = conj(cos(bsxfun(@times,n,acos(abscissa))));
Zstar = Zstar.';

% obtaining the equally spaced points on [-1,1] interval where the
% temperatures are defined and stored in BConds.Dirichlet
a = linspace(-1,1,length(BConds.Dirichlet{LocEdge.id}));

% obtaining the polynomial that interpolates the values in BConds.Dirichlet
if (isnan(BConds.Dirichlet{LocEdge.id}))
    error('local:consistencyChk',...
        'No Dirichlet boundary conditions are defined on edge %d. \n',...
        LocEdge.id);
else
    pol = polyfit(a,BConds.Dirichlet{LocEdge.id},...
        length(BConds.Dirichlet{LocEdge.id})-1);
end

% computing the values of the interpolation polynomials at the abscissas
t = polyval(pol,abscissa);
t = t.'; % non-conjugate transpose

%% Computing the integral on the side
% The integral is the internal product between the flux basis and the 
% applied temperature 
tgi2D = bsxfun(@times, Zstar, t);

% computes the integral
tgi = L/2 * sum(bsxfun(@times,tgi2D,weight.'),2); 

end % End function
