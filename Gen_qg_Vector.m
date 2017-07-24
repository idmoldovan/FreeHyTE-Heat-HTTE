function RHS = Gen_qg_Vector(Edges, Loops, BConds, RHS, abscissa, weight)
% sweeps through the elements and calls the functions that generate the
% qgamma vector in the RHS. qgamma is abbreviated as qg in the code.
% qgamma = Integrate[Û * q_gamma]
% ************************************************************************
% Uses 1D and 2D data structures for storing n and the abscissas
% corresponding to all points/orders that must be computed. For integration
% it constructs a 3D matrix, with each page corresponding to the integrands 
% computed at a Gauss point. The integration is performed as a weighted 
% summation on the pages of this 3D matrix. 


%%

for ii=1:length(Loops.area)
    
    LocLoop = struct('id',ii,'edges',Loops.edges(ii,:), 'center',...
        Loops.center(ii,:),'order',Loops.order(ii),... 
        'insert',Loops.insert(ii),'dim',Loops.dim(ii));
    
    % Computing the qg vector of element ii
    qgi = qg_Vector_i(Edges,LocLoop,BConds, abscissa, weight);
    
    % Inserting the vector in the global RHS vector
    RHS(LocLoop.insert:LocLoop.insert+LocLoop.dim-1) =...
        RHS(LocLoop.insert:LocLoop.insert+LocLoop.dim-1)+ qgi;
    
end

end

function qgi = qg_Vector_i(Edges,LocLoop,BConds, abscissa, weight)

% computes the qg vector of element ii

qgi = zeros(LocLoop.dim,1);

n = -LocLoop.order:LocLoop.order;

for jj = 1:length(LocLoop.edges)  % contour integration
    
    id = LocLoop.edges(jj);  % number of the jj-th edge of the loop
    
    if strcmpi(Edges.type(id),'N')    
        
        LocEdge =  struct('id',id,'nini',Edges.nini(id), 'nfin',Edges.nfin(id),...
            'parametric',Edges.parametric(id,:),'lleft',Edges.lleft(id),...
            'lright',Edges.lright(id));
        
        if LocEdge.lright  % exterior Neumann sides cannot have right loops
            error('local:consistencyChk',...
                'Exterior edge %d cannot have a right element. \n',...
                LocEdge.id);
        end
        
        
    % generate the 3D matrices
        
        [N,A] = ndgrid(n,abscissa);
        
       
        L = sqrt(LocEdge.parametric(3)^2 + LocEdge.parametric(4)^2); % length
        
        % *****************************************************************
        % Getting the local coordinates
        
        loc_x = LocEdge.parametric(1) - LocLoop.center(1) + 0.5 *...
            (A + 1) * LocEdge.parametric(3);  % x & y in local coord
        loc_y = LocEdge.parametric(2) - LocLoop.center(2) + 0.5 *...
            (A + 1) * LocEdge.parametric(4);
        
        % *****************************************************************
     
         R = sqrt(loc_x.^2 + loc_y.^2);  % polar coordinates, local
         Th = atan2(loc_y, loc_x);
     
        
        % Computing the values of the flux at the abscissas:
        
        % obtaining the equally spaced points on [-1,1] interval where the
        % fluxes are defined and stored in BConds.Neumann
        
        a = linspace(-1,1,length(BConds.Neumann{id}));
        
        % obtaining the polynomial that gets the values in BConds.Neumann
        % at the points a
        
        pol = polyfit(a,BConds.Neumann{id},length(BConds.Neumann{id})-1);
        
        % computing the values of "pol" at the abscissas
        
        q = polyval(pol,abscissa);
        q = q.'; % non-conjugate transpose
        
         Ustar = conj(R.^(abs(N)) .* exp(1i*Th.*N)); 
     
        qgi2D = bsxfun(@times, Ustar, q); % this creates the 2D qgi2D matrix,
                                            % one Gauss point per column
        
        qgi = qgi + L/2 * sum(bsxfun(@times,qgi2D,weight.'),2); % computes the integral
        
       
        
    end
    
end

end