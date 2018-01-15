function RHS = Gen_tp_Vector(Edges, Loops, RHS, abscissa, weight)
% GEN_TP_VECTOR iterates through the edges and calls the functions that 
% generate the tparticular vector on the exterior Dirichlet sides. 
% tparticular is abbreviated in the code as tp.
% 
% GEN_TP_VECTOR is called by MAIN***. 
% Input:
% the Edges and Loops data structures, the RHS vector (that is, the free 
% vector of the solving system), and the Gauss-Legendre integration 
% parameters abscissa and weights. 
% Output/Returns to MAIN*** :
% the RHS vector with the tp blocks of all Dirichelet edges inserted at the 
% correct positions (as determined in ASSIGNPARTS).
%
% tp = Integrate[Z * Up]
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
% GEN_TP_VECTOR computes the internal product between the boundary flux
% basis Z on the Dirichlet and interior edges and internally generated 
% temperature. 
% * the boundary flux basis Z is defined by Chebyshev polynomials,
%        Z  = cos(n*ArcCos(abscissa))
% * the internally generated temperature is given by,
%        Up = -(Q/(4*k)).* R.^2
% In the above expressions, n is the line of the current term in vector tp, 
% k is the thermal conductivity coefficient (see INPUTPROC***), Q is the
% internally generated heat, and R is the radial component of the polar
% referential.
%
% Further details on the structure of the solving system are presented in
% reference [4] (Chapters 3 and 4).

%% Sweeping the edges
for ii=1:length(Edges.type)
    
    % tp is computed on the Dirichlet (and interior) edges
    if strcmpi(Edges.type(ii),'D')
        % LocEdge is a structure where the features of the current
        % edge which are directly useful for the calculation of the
        % tp block are stored.
        LocEdge = struct('nini',Edges.nini(ii),'nfin',Edges.nfin(ii),...
            'parametric',Edges.parametric(ii,:),'lleft',...
            Edges.lleft(ii),'lright',Edges.lright(ii),...
            'order',Edges.order(ii),'insert',Edges.insert(ii),...
            'dim',Edges.dim(ii));
        
        if LocEdge.lleft
            id = LocEdge.lleft;
            sign = 1.;
            % LocLoop is a structure where the features of the current
            % element which are directly useful for the calculation of the
            % tp block are stored.
            LocLoop = struct('id',id,'edges',Loops.edges(id,:),...
                'center',Loops.center(id,:),'order',...
                Loops.order(id),'insert',Loops.insert(id),...
                'dim',Loops.dim(id),'k',Loops.material(id,1),...
                'Q',Loops.material(id,2));
            
            
            % Computing the tg vector of edge ii. Function TP_VECTOR_I
            % is a local function defined below.
            tpi = sign*tp_Vector_i(LocEdge,LocLoop, abscissa, weight);
            
            % Inserting the vector in the global RHS vector
            RHS(LocEdge.insert:LocEdge.insert+LocEdge.dim-1) =...
                RHS(LocEdge.insert:LocEdge.insert+LocEdge.dim-1)+tpi;
            
        else                    %there should always be a left element
            error('local:consistencyChk',...
                'No left loop for edge %d. \n', ii);
        end
        
        if LocEdge.lright
            id = LocEdge.lright;
            sign = -1.;
            % LocLoop is a structure where the features of the current
            % element which are directly useful for the calculation of the
            % tp block are stored.
            LocLoop = struct('id',id,'edges',Loops.edges(id,:),...
                'center',Loops.center(id,:),'order',...
                Loops.order(id,1),'insert',Loops.insert(id,1),...
                'dim',Loops.dim(id,1),'k',Loops.material(id,1),...
                'Q',Loops.material(id,2));
            
            % Computing the tg vector of edge ii. Function TP_VECTOR_I
            % is a local function defined below.
            tpi = sign*tp_Vector_i(LocEdge,LocLoop, abscissa, weight);
            
            % Inserting the vector in the global RHS vector
            RHS(LocEdge.insert:LocEdge.insert+LocEdge.dim-1) =...
                RHS(LocEdge.insert:LocEdge.insert+LocEdge.dim-1)+tpi;
        end
    end
    
end
end

%%
function tpi = tp_Vector_i(LocEdge,LocLoop, abscissa, weight)
%  TP_VECTOR_I local function computes the tp vector of the edge ii
% The sides are mapped to a [-1,1] interval to perform the integration

%% Initialization
k = LocLoop.k;
Q = LocLoop.Q;

% n+1 is the current line
n = 0:LocEdge.order;

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
    (abscissa + 1) * LocEdge.parametric(3);
loc_y = LocEdge.parametric(2) - LocLoop.center(2) + 0.5 *...
    (abscissa + 1) * LocEdge.parametric(4);

% Transforming the local Cartesian coordinates into polar.
R = sqrt(loc_x.^2 + loc_y.^2);

%% Computing the integrands at the integration points
% Computing the values of the normal flux basis
Zstar = conj(cos(bsxfun(@times,n,acos(abscissa))));
Zstar = Zstar.';

% Generating the particular solution Up
Up = -(Q/(4*k)).* R.^2;
Up= Up.';

%% Computing the integral on the side
% Computing the internal product
tpi2D = bsxfun(@times, Zstar, Up);

% Computes the integral
tpi = L/2 * sum(bsxfun(@times,tpi2D,weight.'),2); 

end
