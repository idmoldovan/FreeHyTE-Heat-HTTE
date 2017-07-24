function tgbar = Gen_tgbar_Vector(Edges, BConds, abscissa, weight, index)
% Generates the t_{\Gamma} term present in the RHS of the new system. It belongs
% to the edge that was incremented with a single dof, it's a scalar,
% but the  routine is constructed to accommodate its possible 
% future extension to a vectorial quantity.

% ***********************************************************************
% Uses 1D data structures for storing n and abscissas
% corresponding to all points/orders that must be computed. 
% For integration it constructs a 2D matrix, with each column corresponding
% to the integrand computed at a Gauss point. The integration is 
% performed as a weighted summation on the columns of this 2D matrix.

tgbar=zeros(1,1); % change the first term to make it a vector

ii=index;

if (strcmpi(Edges.type(ii),'D') && Edges.lright(ii) == 0 )
    LocEdge  = struct('id',ii,'nini',Edges.nini(ii),...
        'nfin',Edges.nfin(ii),'parametric',Edges.parametric(ii,:),...
        'order',Edges.order(ii),'insert',Edges.insert(ii),...
        'dim',Edges.dim(ii));
    
    % Computing the tgvbar vector of edge ii
    tgbar = tg_Vector_i(LocEdge,BConds, abscissa, weight);
 
end

end


function tgbar = tg_Vector_i(LocEdge, BConds, abscissa, weight)

% compute the tg vector for edge LocEdge
n = LocEdge.order+1; % only one degree is incremented on the boundary basis

L = sqrt(LocEdge.parametric(3)^2 + LocEdge.parametric(4)^2); % length


% Integrating on the side 

% Z* -> the order is 'n'
Zstar = conj(cos(bsxfun(@times,n,acos(abscissa))));
Zstar = Zstar.';

% Computing the values of the enforced temperatures at the abscissas:

% obtaining the equally spaced points on [-1,1] interval where the 
% "temperatures" are defined and stored in BConds.Dirichlet

a = linspace(-1,1,length(BConds.Dirichlet{LocEdge.id}));

% Obtaining the polynomial that has the values given in BConds.Dirichlet

if (isnan(BConds.Dirichlet{LocEdge.id}))
    error('local:consistencyChk',...
        'No Dirichlet boundary conditions are defined on edge %d. \n',...
        LocEdge.id);
else
    pol = polyfit(a,BConds.Dirichlet{LocEdge.id},...
        length(BConds.Dirichlet{LocEdge.id})-1);
end

% computing the values of "pol" at the abscissas

q = polyval(pol,abscissa);
q = q.'; % non-conjugate transpose

% this creates a 2D Qi2D matrix, one Gauss point per column
tgbar2D = bsxfun(@times, Zstar, q);

tgbar = L/2 * sum(bsxfun(@times,tgbar2D,weight.'),2); % computes the integral

end

