function RHS = Gen_tg_Vector(Edges, BConds, RHS, abscissa, weight)
% sweeps through the edges and calls the functions that generate
% tgamma vector on the exterior Dirichlet side. tgamma is abbreviated in
% the code as tg.
% tgamma = Integrate[Z * T_gamma]
% ***********************************************************************
% Uses 1D data structures for storing n and abscissas
% corresponding to all points/orders that must be computed. 
% For integration it constructs a 2D matrix, with each column corresponding
% to the integrand computed at a Gauss point. The integration is 
% performed as a weighted summation on the columns of this 2D matrix.

for ii=1:length(Edges.type)
    
    if (strcmpi(Edges.type(ii),'D') && Edges.lright(ii) == 0 )
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

function tgi = tg_Vector_i(LocEdge, BConds, abscissa, weight)

n = 0:LocEdge.order;
L = sqrt(LocEdge.parametric(3)^2 + LocEdge.parametric(4)^2); % length

% Z* -> the order is 'n'
Zstar = conj(cos(bsxfun(@times,n,acos(abscissa))));
Zstar = Zstar.';

% Computing the values of the enforced "temperatures" at the abscissas

% obtaining the equally spaced points on [-1,1] interval where the 
% "temperatures" are defined and stored in BConds.Dirichlet

a = linspace(-1,1,length(BConds.Dirichlet{LocEdge.id}));

% Obtaining the polynomial that has the values given in BConds.Dirichlet
% at the points a

if (isnan(BConds.Dirichlet{LocEdge.id}))
    error('local:consistencyChk',...
        'No Dirichlet boundary conditions are defined on edge %d. \n',...
        LocEdge.id);
else
    pol = polyfit(a,BConds.Dirichlet{LocEdge.id},...
        length(BConds.Dirichlet{LocEdge.id})-1);
end

% computing the values of "pol" at the abscissas
t = polyval(pol,abscissa);
t = t.'; % non-conjugate transpose

% this creates a 2D tgi2D matrix, one Gauss point per column
tgi2D = bsxfun(@times, Zstar, t);

tgi = L/2 * sum(bsxfun(@times,tgi2D,weight.'),2); % computes the integral

end