function Abar = Gen_Bbar_Matrix(LHS0,Edges,Loops,Dim,abscissa,weight,index)
% Computes the BBar matrix for the boundary to be refined and stores it in
% the augmented matrix ABar
%%
Abar=zeros(Dim+1,Dim+1);    % Creating the Abar matrix, the new temporary system.
% '+1' comes from the assumption that only one
% degree will be added to the boundary.
Abar(1:Dim,1:Dim)=LHS0;			              % Allocating the "old" A matrix inside.

ii=index;

Edges.order(ii)=Edges.order(ii)+1; % increments the boundary degree
Edges.dim(ii)=Edges.dim(ii)+1;

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
            'dim',Loops.dim(id));
        
        
        %%
        % Computing the B matrix of the edge ii, left loop
        Bbardi = sign*Bbar_Matrix_i(LocEdge,LocLoop,abscissa,weight);
        
        % Inserting the new B matrix in the global Abar matrix
        Abar(LocLoop.insert:LocLoop.insert+LocLoop.dim-1,Dim+1) = Bbardi;
        % Inserting the conjugate transposed in the global LHS matrix
        Abar(Dim+1,LocLoop.insert:LocLoop.insert+LocLoop.dim-1) = Bbardi';
        
        
    end
    
    if LocEdge.lright
        id = LocEdge.lright;
        sign = -1.;
        LocLoop = struct('id',id,'edges',Loops.edges(id,:),...
            'center',Loops.center(id,:),'order',...
            Loops.order(id,1),'insert',Loops.insert(id,1),...
            'dim',Loops.dim(id,1));
        
        %%
        % Computing the Bbar matrix of the edge ii, right loop
        Bbardi = sign*Bbar_Matrix_i(LocEdge,LocLoop,abscissa,weight);
        
        % Inserting the matrix in the global LHS matrix
        Abar(LocLoop.insert:LocLoop.insert+LocLoop.dim-1,Dim+1) = Bbardi;
        % Inserting the conjugate transposed in the global LHS matrix
        Abar(Dim+1,LocLoop.insert:LocLoop.insert+LocLoop.dim-1) = Bbardi';
        
        
    end
end

Edges.order(ii)=Edges.order(ii)-1; % un-incrementing the order & dimension
Edges.dim(ii)=Edges.dim(ii)-1;

end
%%
function Bbardi = Bbar_Matrix_i(LocEdge,LocLoop,abscissa,weight)

n = -LocLoop.order:LocLoop.order;
m = LocEdge.order;                    % varies between 0 and M
a = abscissa;

L = sqrt((LocEdge.parametric(3))^2 + (LocEdge.parametric(4))^2); % length of the side

% Creating the 3D grids...
[N,M,A]=ndgrid(n,m,a);

% Getting the r, th, nr, nth for all Gauss points

loc_x = LocEdge.parametric(1) - LocLoop.center(1) + 0.5 *...
    (A + 1) * LocEdge.parametric(3);  % x & y in local ccord
loc_y = LocEdge.parametric(2) - LocLoop.center(2) + 0.5 *...
    (A + 1) * LocEdge.parametric(4);

R = sqrt(loc_x.^2 + loc_y.^2);  % polar coordinates, local
Th = atan2(loc_y, loc_x);


% U* -> the order is 'n'
Ustar = conj(R.^abs(N) .* exp(1i*Th.*N));

% this creates the 3D Bdi matrix
Bdi3D = Ustar .*...
    cos((M).*acos(A)); % Chebyshev functions


W3D(1,1,:) = weight;
Bbardi = (L/2)*sum(bsxfun(@times,W3D,Bdi3D),3);

end