function tgbar = Gen_tgbar_Vector(Edges, BConds, abscissa, weight, index)
% GEN_TGBAR_VECTOR generates t_Gamma term present in the RHS of the 
% system (5.10) of reference[4]. The term in the present formulation is 
% a scalar, but the  routine is constructed to accommodate its possible 
% future extension to a vectorial quantity. Its definition is given by
% expression (5.4) of reference [4].
%
% GEN_TGBAR_VECTOR is called by EDGEREFINEMENT. 
%
% Input:
% the Edges, BConds structures, the id of the refined boundary (index), 
% and the Gauss-Legendre integration parameters abscissa and weights. 
% Output/Returns to EDGEREFINEMENT:
%   the tgbar value (vector) for the refined edge.
% Calls: TG_VECTOR_I to compute tgbar for the refined edge
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

% Intialize tpbar with zero. In case of a vector change the first 
% term to make it a vector
tgbar=zeros(1,1); 

% Identify the refined edge
ii=index;

% Exterior Dirichlet boundaries have no right element
if (strcmpi(Edges.type(ii),'D') && Edges.lright(ii) == 0 )
    % LocEdge is a structure where the features of the current
    % edge which are directly useful for the calculation of the
    % tgbar are stored.    
    LocEdge  = struct('id',ii,'nini',Edges.nini(ii),...
        'nfin',Edges.nfin(ii),'parametric',Edges.parametric(ii,:),...
        'order',Edges.order(ii),'insert',Edges.insert(ii),...
        'dim',Edges.dim(ii));
    
    % Computing the tgvbar vector of edge ii. Function TG_VECTOR_I 
    % is a local funtion defined below
    tgbar = tg_Vector_i(LocEdge,BConds, abscissa, weight);
 
end

end

%%
function tgbar = tg_Vector_i(LocEdge, BConds, abscissa, weight)
%  TG_VECTOR_I local function computes the tg vector of the edge ii
% The sides are mapped to a [-1,1] interval to perform the integration

% compute the tg vector for edge LocEdge
% only one degree is incremented on the boundary basis
n = LocEdge.order+1; 

% Compute the length of the edge 
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
q = polyval(pol,abscissa);
q = q.'; % non-conjugate transpose

%% Computing the integral on the side
% The integral is the internal product between the flux basis and the 
% applied temperature 
tgbar2D = bsxfun(@times, Zstar, q);

% computes the integral
tgbar = L/2 * sum(bsxfun(@times,tgbar2D,weight.'),2); 

end

