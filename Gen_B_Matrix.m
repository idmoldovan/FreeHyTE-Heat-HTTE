function LHS = Gen_B_Matrix(Edges,Loops,LHS,abscissa,weight)
% sweeps through the edges and calls the functions that generate the
% B matrix in the LHS
% B = Integrate [ Û * Z ]

% The Chebyshev polynomial T(n,x), or Chebsyhev
% polynomial of the first kind, may be defined,
% for 0<=n, and -1<=x<=+1 by:
% -> cos(t)=x; t=acos(x);
% -> T(n,x) = cos(n*t); T(n,x) = cos(n*acos(x));
% For any value of x, T(n,x) may be evaluated by
% three term recurrence:
% T(0,x)=1;
% T(1,x)=x;
% T(n+1,x)=(2*x)*T(n,x)-T(n-1,x);

%%

for ii=1:length(Edges.type)
    
    if strcmpi(Edges.type(ii),'D')
        LocEdge = struct('nini',Edges.nini(ii),'nfin',Edges.nfin(ii),...
            'parametric',Edges.parametric(ii,:),'lleft',...
            Edges.lleft(ii),'lright',Edges.lright(ii),...
            'order',Edges.order(ii),'insert',Edges.insert(ii),...
            'dim',Edges.dim(ii));
        
        if LocEdge.lleft %there should always be a left element
            id = LocEdge.lleft;
            sign = 1.;
            LocLoop = struct('id',id,'edges',Loops.edges(id,:),...
                'center',Loops.center(id,:),'order',...
                Loops.order(id),'insert',Loops.insert(id),...
                'dim',Loops.dim(id));
            
            % Computing the B matrix of the edge ii, left loop
            Bdi = sign*B_Matrix_i(LocEdge,LocLoop,abscissa,weight);
            
            % Inserting the matrix in the global LHS matrix
            LHS(LocLoop.insert:LocLoop.insert+LocLoop.dim-1,...
                LocEdge.insert:LocEdge.insert+LocEdge.dim-1) = -Bdi;
            % Inserting the conjugate transposed in the global LHS matrix
            LHS(LocEdge.insert:LocEdge.insert+LocEdge.dim-1,...
                LocLoop.insert:LocLoop.insert+LocLoop.dim-1) = -Bdi';
            
        else
            error('local:consistencyChk',...
                'No left loop for edge %d. \n', ii);
        end
        
        if LocEdge.lright %there's only right elements when internal boundries exist
            id = LocEdge.lright;
            sign = -1.;
            LocLoop = struct('id',id,'edges',Loops.edges(id,:),...
                'center',Loops.center(id,:),'order',...
                Loops.order(id,1),'insert',Loops.insert(id,1),...
                'dim',Loops.dim(id,1));
            
            % Computing the B matrix of the edge ii, right loop
            Bdi = sign*B_Matrix_i(LocEdge,LocLoop,abscissa,weight);
            
            % Inserting the matrix in the global LHS matrix
            LHS(LocLoop.insert:LocLoop.insert+LocLoop.dim-1,...
                LocEdge.insert:LocEdge.insert+LocEdge.dim-1) = -Bdi;
            % Inserting the conjugate transposed in the global LHS matrix
            LHS(LocEdge.insert:LocEdge.insert+LocEdge.dim-1,...
                LocLoop.insert:LocLoop.insert+LocLoop.dim-1) = -Bdi';
            
        end
    end
end
end
%%
function Bdi = B_Matrix_i(LocEdge,LocLoop,abscissa,weight)

n = -LocLoop.order:LocLoop.order;
m = 0:LocEdge.order;
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
Bdi = (L/2)*sum(bsxfun(@times,W3D,Bdi3D),3);

end

