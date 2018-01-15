function Ep = Energy_p(Loops, Edges, X, abscissa, weight)
% ENERGY_P iterates on the elements and calls the function that generates  
% the Ep part of the solution energy.
%
% Input: Loops and Egdes structures, solution X, and the Gauss-Legendre
% integration parameters abscissa and weights. 
% Output:  Ep part of the solution energy
% Calls: EPI_VECTOR_I function
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
% The computation of the solution energy is covered in Section 3.9 of
% reference [4]. Its calculation requires no domain intergation. Using 
% the notations of reference [4], 
% Ep = Contour_Integral[Up*n*S] * X
%

%% Initialization.

Ep=0;

for ii=1:length(Loops.area)
    
    % LocLoop is a structure where the features of the current
    % element which are directly useful for the calculation of the
    % Ep are stored.    
    LocLoop = struct('id',ii,'edges',Loops.edges(ii,:), 'center',...
        Loops.center(ii,:),'order',Loops.order(ii),...
        'insert',Loops.insert(ii),'dim',Loops.dim(ii),...
        'k',Loops.material(ii,1),'Q',Loops.material(ii,2));
 
        %% Computing the Ep energy of loop ii. 
        % Function EPI_VECTOR_I is a local function defined below. 
        Epi = Epi_Vector_i(Edges,LocLoop, abscissa, weight);
        
        % Extracting the solution corresponding to the LocLoop from the
        % global solution vector.
        Xi(:,1) = X(LocLoop.insert:LocLoop.insert+LocLoop.dim-1); 
        
        % Incrementing the global Ep energy with the LocLoop contribution.
        Ep= Ep + Epi*Xi;                                          
        
        clear Xi;
 
end
end

%%
function Epi = Epi_Vector_i(Edges,LocLoop, abscissa, weight)
% EPI_VECTOR_I local function computes the Epi for the loop ii
% The sides are mapped to a [-1,1] interval to perform the integration

%% Initialization
Epi = zeros(1,LocLoop.dim);

k = LocLoop.k;
Q = LocLoop.Q;

% Vector containing the orders of the basis
n = -LocLoop.order:LocLoop.order;

%% Sweeping the edges for contour integration
for jj = 1:length(LocLoop.edges)  
    
    id = LocLoop.edges(jj);  % number of the jj-th edge of the loop
    
    % LocEdge is a structure where the features of the current
    % edge which are directly useful for the calculation of the
    % of the Epi energy
    LocEdge =  struct('id',id,'nini',Edges.nini(id), 'nfin',Edges.nfin(id),...
        'parametric',Edges.parametric(id,:),'lleft',Edges.lleft(id),...
        'lright',Edges.lright(id));
    
    % Constructing the integration grid
    [A,N] = ndgrid(abscissa,n);
    
    %% Generating the geometric data
    % The following code transforms the abscissa coordinates, expressed in
    % the [-1,1] referential, to the polar coordinates required to compute
    % the values of the basis functions.
    
    % Computing the length of the current edge
    L = sqrt(LocEdge.parametric(3)^2 + LocEdge.parametric(4)^2); 
    
    % Transforming the edge abscissa into local coordinates. The local
    % referential is centered in the barycenter of the element, its axes
    % aligned with the Cartesian axes of the global referential.
    loc_x = LocEdge.parametric(1) - LocLoop.center(1) + 0.5 *...
        (A + 1) * LocEdge.parametric(3);  
    loc_y = LocEdge.parametric(2) - LocLoop.center(2) + 0.5 *...
        (A + 1) * LocEdge.parametric(4);
    
    % Transforming the local Cartesian coordinates into polar.
    R = sqrt(loc_x.^2 + loc_y.^2);  
    Th = atan2(loc_y, loc_x);
    
    % Computing the components of the outward normal in the Cartesian
    % directions. 
    nx = LocEdge.parametric(4) / L;   % normal in (local/global) x & y
    ny = -1* LocEdge.parametric(3) / L;
    if LocEdge.lright==LocLoop.id  % if the element is on the right,
        nx = -nx;                  % change the sign of the normal
        ny = -ny;
    end
    
    % Computing the components of the outward normal in the polar
    % directions.
    NR = nx * cos(Th) + ny * sin(Th);
    NTh = -1*nx * sin(Th) + ny * cos(Th);
    
    % Polar components of the flux basis S
    SR = -k * abs(N) .* R.^(abs(N)-1) .* exp(1i*Th.*N);
    STh = -k * 1i*N.*R.^(abs(N)-1) .* exp(1i*Th.*N);
    % Boundary normal component of the flux basis
    NS = NR.*SR + NTh.*STh;  
    
    % Computing the particlular solution basis
    Up = -(Q/(4*k)).* R.^2;  
    
   % computes the integral  
    Ekpi2D = bsxfun(@times, Up, NS); 
    Epi = Epi + L/2 * sum(bsxfun(@times,Ekpi2D,weight),1);
    
end
end
