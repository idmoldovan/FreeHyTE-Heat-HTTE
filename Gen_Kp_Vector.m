function RHS = Gen_Kp_Vector(Edges, Loops, RHS, abscissa, weight)
% sweeps through the elements and calls the functions that generate the
% Kp vector in the RHS
%Kp = Integrate[Û * n * k * DeltaUp]
% ************************************************************************

for ii=1:length(Loops.area)
    
    LocLoop = struct('id',ii,'edges',Loops.edges(ii,:), 'center',...
        Loops.center(ii,:),'order',Loops.order(ii),...
        'insert',Loops.insert(ii),'dim',Loops.dim(ii),...
        'Q',Loops.material(ii,2));
    
    % Computing the qg vector of element ii
    Kpi = Kpi_Vector_i(Edges,LocLoop, abscissa, weight);
    
    % Inserting the vector in the global RHS vector
    RHS(LocLoop.insert:LocLoop.insert+LocLoop.dim-1) = ...
        RHS(LocLoop.insert:LocLoop.insert+LocLoop.dim-1) - Kpi;
    
end
end


function Kpi = Kpi_Vector_i(Edges,LocLoop, abscissa, weight)

% computes the Kp vector of element ii

Kpi = zeros(LocLoop.dim,1);

Q = LocLoop.Q;

n = -LocLoop.order:LocLoop.order;

for jj = 1:length(LocLoop.edges)  % contour integration
    
    id = LocLoop.edges(jj);  % number of the jj-th edge of the loop
    
    LocEdge =  struct('id',id,'nini',Edges.nini(id), 'nfin',Edges.nfin(id),...
        'parametric',Edges.parametric(id,:),'lleft',Edges.lleft(id),...
        'lright',Edges.lright(id));
        
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
    
    
    nx = LocEdge.parametric(4) / L;   % normal in (local/global) x & y
    ny = -1* LocEdge.parametric(3) / L;
    
    if LocEdge.lright==LocLoop.id  % if the element is on the right,
        nx = -nx;                  % change the sign of the normal
        ny = -ny;
    end
    
    Ustar = conj(R.^(abs(N)) .* exp(1i*Th.*N)); %Calculating base U*
    
    SPR = (Q/2)*R; %SPTh=0
    
    NR = nx * cos(Th) + ny * sin(Th); % normal in local r-th (only R is needed)
    
    NSP = NR.*SPR;
    
    Kpi2D = bsxfun(@times, Ustar, NSP);
    
    Kpi = Kpi + L/2 * sum(bsxfun(@times,Kpi2D,weight.'),2);
    
end
end