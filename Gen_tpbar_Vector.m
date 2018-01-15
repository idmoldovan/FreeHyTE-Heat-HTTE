function tpbar = Gen_tpbar_Vector(Edges, Loops, abscissa, weight, index)
% GEN_TPBAR_VECTOR computes the tpbar term present in the RHS of the 
% system (5.10) of reference[4]. The term in the present formulation is 
% a scalar, but the  routine is constructed to accommodate its possible 
% future extension to a vectorial quantity. Its definition is given by
% expression (5.5) of reference [4].
%
% GEN_TPBAR_VECTOR is called by EDGEREFINEMENT. 
%
% Input:
% the Edges, Loops structures, the id of the refined boundary (index), 
% and the Gauss-Legendre integration parameters abscissa and weights. 
% Output/Returns to EDGEREFINEMENT:
%   the tpbar value (vector) for the refined edge.
% Calls: TP_VECTOR_I to compute tgbar for the refined edge
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

% Intialize tpbar with zero. In case of a vector change the first 
% term to make it a vector
tpbar=zeros(1,1); 

% Identify the refined edge
ii=index;

% tpbar is computed on the Dirichlet (and interior) edges
if strcmpi(Edges.type(ii),'D')
    % LocEdge is a structure where the features of the current
    % edge which are directly useful for the calculation of the
    % tpbar are stored.    
    LocEdge = struct('nini',Edges.nini(ii),'nfin',Edges.nfin(ii),...
        'parametric',Edges.parametric(ii,:),'lleft',...
        Edges.lleft(ii),'lright',Edges.lright(ii),...
        'order',Edges.order(ii),'insert',Edges.insert(ii),...
        'dim',Edges.dim(ii));
    
    % All boundaries must have a left element
    if LocEdge.lleft
        id = LocEdge.lleft;
        sign = 1.;
        
        % LocLoop is a structure where the features of the current element 
        % which are directly useful for the calculation of the tpbar are stored.         
        LocLoop = struct('id',id,'edges',Loops.edges(id,:),...
            'center',Loops.center(id,:),'order',...
            Loops.order(id),'insert',Loops.insert(id),...
            'dim',Loops.dim(id),'k',Loops.material(id,1),...
            'Q',Loops.material(id,2));
               
        % Computing the tpbar vector of edge ii. Function TP_VECTOR_I 
        % is a local function defined below. 
        tpbar = sign*tp_Vector_i(LocEdge,LocLoop, abscissa, weight);
        
    else                    %there should always be a left element
        error('local:consistencyChk',...
            'No left loop for edge %d. \n', ii);
    end
    
    % There may or may not be a right element
    if LocEdge.lright
        id = LocEdge.lright;
        sign = -1.;
        
        % LocLoop is a structure where the features of the current element
        % which are directly useful for the calculation of the tpbar are stored.   
        LocLoop = struct('id',id,'edges',Loops.edges(id,:),...
            'center',Loops.center(id,:),'order',...
            Loops.order(id,1),'insert',Loops.insert(id,1),...
            'dim',Loops.dim(id,1),'k',Loops.material(id,1),...
            'Q',Loops.material(id,2));
        
        % Computing the tpbar vector of edge ii. Function TP_VECTOR_I 
        % is a local function defined below. 
        tpi = sign*tp_Vector_i(LocEdge,LocLoop, abscissa, weight);
        
        % Updating tpbar
        tpbar = tpbar + tpi;
        
    end
    
end

end


%%
function tpi = tp_Vector_i(LocEdge,LocLoop, abscissa, weight)
% TP_VECTOR_I local function computes the tp term of edge LocEdge and loop
% LocLoop. The sides are mapped to a [-1,1] interval to perform the 
% integration

k = LocLoop.k;
Q = LocLoop.Q;

% only one degree is incremented on the boundary basis
n = LocEdge.order+1; 

%% Generating the geometric data
% The following code transforms the abscissa coordinates, expressed in
% the [-1,1] referential, to the polar coordinates required to compute
% the values of the basis functions.

% Compute the length of the edge LocEdge
L = sqrt(LocEdge.parametric(3)^2 + LocEdge.parametric(4)^2); 

% Transforming the edge abscissa into local coordinates. The local
% referential is centered in the barycenter of the element, its axes
% aligned with the Cartesian axes of the global referential. 
loc_x = LocEdge.parametric(1) - LocLoop.center(1) + 0.5 *...
    (abscissa + 1) * LocEdge.parametric(3);
loc_y = LocEdge.parametric(2) - LocLoop.center(2) + 0.5 *...
    (abscissa + 1) * LocEdge.parametric(4);

% Transforming the local Cartesian coordinates into polar. Th is not needed
R = sqrt(loc_x.^2 + loc_y.^2);  

%% Computing the integrands at the integration points
% Computing the values of the normal flux basis
Zstar = conj(cos(bsxfun(@times,n,acos(abscissa))));
Zstar = Zstar.';

% Generating the particular solution basis Up
Up = -(Q/(4*k)).* R.^2;
Up= Up.';

%% Computing the integral on the side
% Computing the internal product
tpi2D = bsxfun(@times, Zstar, Up);

% computes the integral
tpi = L/2 * sum(bsxfun(@times,tpi2D,weight.'),2); 

end
