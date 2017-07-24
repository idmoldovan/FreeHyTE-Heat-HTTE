function Ep = Energy_p( Loops, Edges, X, abscissa, weight)
% sweeps through the elements and calls the function that generates the 'Ep'
% value which is a section in the energy computation.
% Ep = Integrate [Up*n*S] * X

%%
%Initializing Ep 
Ep=0;

for ii=1:length(Loops.area)
    
    LocLoop = struct('id',ii,'edges',Loops.edges(ii,:), 'center',...
        Loops.center(ii,:),'order',Loops.order(ii),...
        'insert',Loops.insert(ii),'dim',Loops.dim(ii),...
        'k',Loops.material(ii,1),'Q',Loops.material(ii,2));
 
        Epi = Epi_Vector_i(Edges,LocLoop, abscissa, weight);
        Xi(:,1) = X(LocLoop.insert:LocLoop.insert+LocLoop.dim-1); 
        Ep= Ep + Epi*Xi;                                          
        
        clear Xi;
        
    
end
end

%%
function Epi = Epi_Vector_i(Edges,LocLoop, abscissa, weight)

Epi = zeros(1,LocLoop.dim);

k = LocLoop.k;
Q = LocLoop.Q;

n = -LocLoop.order:LocLoop.order;

for jj = 1:length(LocLoop.edges)  % contour integration
    
    id = LocLoop.edges(jj);  % number of the jj-th edge of the loop
    
    LocEdge =  struct('id',id,'nini',Edges.nini(id), 'nfin',Edges.nfin(id),...
        'parametric',Edges.parametric(id,:),'lleft',Edges.lleft(id),...
        'lright',Edges.lright(id));
    
    [A,N] = ndgrid(abscissa,n);
    
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
    
    % normal in local r-th
    NR = nx * cos(Th) + ny * sin(Th);
    NTh = -1*nx * sin(Th) + ny * cos(Th);
    
    
    % S -> the order is 'm'
    SR = -k * abs(N) .* R.^(abs(N)-1) .* exp(1i*Th.*N);
    STh = -k * 1i*N.*R.^(abs(N)-1) .* exp(1i*Th.*N);
    
    NS = NR.*SR + NTh.*STh;  %Normal to each boundry
    
    Up = -(Q/(4*k)).* R.^2;  %Particular solution
    
    Ekpi2D = bsxfun(@times, Up, NS);
    
    Epi = Epi + L/2 * sum(bsxfun(@times,Ekpi2D,weight),1);
    
end
end