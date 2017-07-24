function LHS = Gen_K_Matrix(Edges, Loops, LHS, abscissa, weight)
% sweeps through the elements and calls the functions that generate the
% K matrix in the LHS
% K = Integrate [ Û n k DeltaU ]
% ************************************************************************

for ii=1:length(Loops.area)
    
    LocLoop = struct('id',ii,'edges',Loops.edges(ii,:), 'center',...
        Loops.center(ii,:),'order',Loops.order(ii,:),...
        'insert',Loops.insert(ii,:),'dim',Loops.dim(ii,:),...
        'k',Loops.material(ii,1));
    
    % Computing the Ddd matrix of element ii
    Ki = K_Matrix_i(LocLoop, Edges, abscissa, weight);
    
    % Inserting the matrix in the global LHS matrix
    LHS(LocLoop.insert:LocLoop.insert+LocLoop.dim-1,...
        LocLoop.insert:LocLoop.insert+LocLoop.dim-1) = Ki;
    
end

end

function Ki = K_Matrix_i(LocLoop, Edges, abscissa, weight)

% computes the Ddd matrix of element ii

Ki = zeros(LocLoop.dim);

k = LocLoop.k;

n = -LocLoop.order:LocLoop.order;
m = -LocLoop.order:LocLoop.order;

for jj = 1:length(LocLoop.edges)  % contour integration
    
    id = LocLoop.edges(jj);  % number of the jj-th edge of the loop
    LocEdge =  struct('id',id,'nini',Edges.nini(id), 'nfin',Edges.nfin(id),...
        'parametric',Edges.parametric(id,:),'lleft',Edges.lleft(id),...
        'lright',Edges.lright(id));
    
    L = sqrt(LocEdge.parametric(3)^2 + LocEdge.parametric(4)^2); % length
    
    % *****************************************************************
    
    % generate the 3D matrices
    [N,M,A] = ndgrid(n,m,abscissa);
    
    % Getting the r, th, nr, nth for all Gauss points
    
    loc_x = LocEdge.parametric(1) - LocLoop.center(1) + 0.5 *...
        (A + 1) * LocEdge.parametric(3);  % x & y in local ccord
    loc_y = LocEdge.parametric(2) - LocLoop.center(2) + 0.5 *...
        (A + 1) * LocEdge.parametric(4);
    
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
    
    % Integrating on the side (the side integration is fully vectorialized)
    
    % U* -> the order is 'n'
    Ustar = conj(R.^abs(N) .* exp(1i*Th.*N));
    
    % S -> the order is 'm'
    SR = -k * abs(M) .* R.^(abs(M)-1) .* exp(1i*Th.*M);
    STh = -k * 1i*M.*R.^(abs(M)-1) .* exp(1i*Th.*M);
    
    %Calculating DeltaU
    NS = NR.*SR + NTh.*STh;
    
    Ki3D = Ustar .* NS; % this creates the 3D Ki matrix,
                        % one Gauss point per page
    
    w3D(1,1,:) = weight;
    
    Ki = Ki + L/2 * sum(bsxfun(@times,Ki3D,w3D),3); % computes the integral
    
    
end

end