function tpbar = Gen_tpbar_Vector(Edges, Loops, abscissa, weight, index)
% Generates the t_{\Gamma_p} term present in the RHS of the new system. It belongs
% to the edge that was incremented with a single dof, it's a scalar,
% but the  routine is constructed to accommodate its possible 
% future extension to a vectorial quantity.
% t_{\Gamma_p} = Integrate[Znew * Up]

tpbar=zeros(1,1); % change the first term to make it a vector

ii=index;

if strcmpi(Edges.type(ii),'D')
    LocEdge = struct('nini',Edges.nini(ii),'nfin',Edges.nfin(ii),...
        'parametric',Edges.parametric(ii,:),'lleft',...
        Edges.lleft(ii),'lright',Edges.lright(ii),...
        'order',Edges.order(ii),'insert',Edges.insert(ii),...
        'dim',Edges.dim(ii));
    
    if LocEdge.lleft
        id = LocEdge.lleft;
        sign = 1.;
        LocLoop = struct('id',id,'edges',Loops.edges(id,:),...
            'center',Loops.center(id,:),'order',...
            Loops.order(id),'insert',Loops.insert(id),...
            'dim',Loops.dim(id),'k',Loops.material(id,1),...
            'Q',Loops.material(id,2));
        
        
        % Computing the tg vector of edge ii
        tpbar = sign*tp_Vector_i(LocEdge,LocLoop, abscissa, weight);
        
    else                    %there should always be a left element
        error('local:consistencyChk',...
            'No left loop for edge %d. \n', ii);
    end
    
    if LocEdge.lright
        id = LocEdge.lright;
        sign = -1.;
        LocLoop = struct('id',id,'edges',Loops.edges(id,:),...
            'center',Loops.center(id,:),'order',...
            Loops.order(id,1),'insert',Loops.insert(id,1),...
            'dim',Loops.dim(id,1),'k',Loops.material(id,1),...
            'Q',Loops.material(id,2));
        
        % Computing the tg vector of edge ii
        tpi = sign*tp_Vector_i(LocEdge,LocLoop, abscissa, weight);
        
        % Inserting the vector in the global RHS vector
        tpbar = tpbar + tpi;
        
    end
    
end

end


function tpi = tp_Vector_i(LocEdge,LocLoop, abscissa, weight)

k = LocLoop.k;
Q = LocLoop.Q;

n = LocEdge.order+1; % only one degree is incremented on the boundary basis
L = sqrt(LocEdge.parametric(3)^2 + LocEdge.parametric(4)^2); % length

% Integrating on the side (the side integration is fully vectorialized)

% Z* -> the order is 'n'
Zstar = conj(cos(bsxfun(@times,n,acos(abscissa))));
Zstar = Zstar.';

% x & y in local coord
loc_x = LocEdge.parametric(1) - LocLoop.center(1) + 0.5 *...
    (abscissa + 1) * LocEdge.parametric(3);
loc_y = LocEdge.parametric(2) - LocLoop.center(2) + 0.5 *...
    (abscissa + 1) * LocEdge.parametric(4);

% *****************************************************************

R = sqrt(loc_x.^2 + loc_y.^2);  % polar coordinates, local
                                %Th is not needed
Up = -(Q/(4*k)).* R.^2;
Up= Up.';

% this creates a 2D tgi2D matrix, one Gauss point per column
tpi2D = bsxfun(@times, Zstar, Up);

tpi = L/2 * sum(bsxfun(@times,tpi2D,weight.'),2); % computes the integral

end