function MainTri
% MainTri main program of the FreeHyTE - Steady State Heat Conduction 
% module with triangular meshes
% MAINTRI is called upon exiting the data input (GUI) phase of the module.
% It is used to launch all functions required for the solution of the 
% heat transfer problem and centralize all data they provide. 
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

%% ******** PRE-PROCESSING ********
% Launch pre-processing routine: InputProcReg
% Processes the input data from GUI, generates the mesh and topological data
% structures
% * NGP is the number of Gauss points for the line integration;
% * Nodes is a (NNODE x 2) matrix, where NNODE is the number of nodes in 
% mesh. It stores the coordinates of each node;
% * Edges, Loops and BConds are data structures storing information
% on the edges, finite elements (loops) and boundary conditions,
% respectively. They are documented in Section 5.3 of reference [2];
[NGP, Nodes, Edges, Loops, BConds] = InputProcTri;

% Check if the initial domain basis aproximation satisfies the 'kinematic'
% ideterminacy condition and correct if necessary. The kinematic
% indeterminacy condition is covered in Section 3.4.2 of reference [2],
% section 4.5.2 of reference [3], and section 3.8 of reference [4].
[Loops,Edges] = CheckMinOrders(Loops,Edges); 

% Load the analysis type option
% RunOption = 1 - single analysis
% RunOption = 2 - adaptive algorithm
load('HeatStructDef','RunOption');

%% ******** INITIAL ANALYSIS ********
% The "INITIAL ANALYSIS" is completed in two cases: 
%   RunOption = 1 for a single analysis case (without adaptive analysis) 
%   RunOption = 2 for adaptive analysis 
%
%
% ASSIGNPARTS maps the finite element solving system and assigns entries
% and dimensions to elements and sides. The information is used by the
% functions that generate the blocks of the solving system to insert them
% at the right positions.
% * Dim is the total dimension of the finite element solving system;
% * the mapping of the solving system is covered in Section 6.2 of
% reference [2].
[Edges,Loops,Dim] = AssignParts(Edges, Loops); 

% Stores the coordinate of the beginning of the Edges' part in the global 
% system
Edges.Edgeinsert=Loops.insert(end)+Loops.dim(end); 

% Initialization of the matrix of coefficients and the free vector
% * LHS is the matrix of coefficients of the solving system;
% * RHS is the free vector of the solving system.
LHS = zeros(Dim);
RHS = zeros(Dim,1);

% Initialization of Gauss-Legendre weights & abscissas (on a -1:1
% interval). gauleg is a external routine, written by Greg von Winckel.
[abscissa, weight] = gauleg(NGP, -1, 1);

%% Generation of the solving system
% Generation & allocation of the blocks in the matrix of coefficients.
% General mapping of the matrix of coefficients of hybrid-Trefftz elements:
%   _____________________________________
%  |       |        |  <------>  |       |
%  |  K    |    B   |            |-Kp-qg |
%  |_______|________|            |_______|
%  |       |        |            |       |
%  |   B   |    0   |            |tg-tp  |
%  |_______|________|            |_______|
%
% The following functions generate the coefficients blocks for each finite
% element and essential boundary and insert them in LHS at the right place,
% according to the mapping information generated in ASSIGNPARTS. 

% The explicit expressions of the stiffness and boundary matrices are given
% in Chapter 3 of reference [4].
LHS = Gen_K_Matrix(Edges, Loops, LHS, abscissa, weight);
LHS = Gen_B_Matrix(Edges, Loops, LHS, abscissa, weight);

% Store the initial dimension of the solving system in GDL0
GDL0 = sum(Edges.dim(:))+sum(Loops.dim(:)); 

% Generation & allocation of the blocks in the free vector.
% The general mapping of the free vector is consistent to that of the
% coefficient matrix (above). 
%
% The following functions generate the free vectors for each finite
% element and essential boundary and insert them at the right place,
% according to the mapping information generated in ASSIGNPARTS. 
%
% The explicit expressions of the free vectors are given in Chapter 3 of 
% reference [4].
RHS = Gen_tg_Vector(Edges, BConds, RHS, abscissa, weight);
RHS = Gen_qg_Vector(Edges, Loops,BConds, RHS, abscissa, weight);
RHS = Gen_Kp_Vector(Edges, Loops, RHS, abscissa, weight);
RHS = Gen_tp_Vector(Edges,Loops, RHS, abscissa, weight);

%% Solution of the solving system
% The solution of the solving system depends on its numerical stability.
% Ill-conditioned systems are solved by a truncated singular value 
% decomposition technique. The truncating singular value is detected by
% analysing the alignment of the singular values of the coefficients
% matrix: a sudden drop in this alignment (decrease of the singular value)
% indicates that the system may be ill-conditioned.
% Systems that reveal no suddem drops in their singular values are scaled
% and solved using the default Matlab solver (mldivide).

% Analysing the singular values alignment to detect outliers. Two
% procedures for the outlier detection can be used:
% * fit - analyses the misalignment of the logarithms of the singular
% values; and,
% * jump - analyses the relative difference (jump) between consecutive
% singular values.
% Details on the procedures are given in Section 6.2.2 of reference [2] and
% section 4.4.4 of reference [4].
% Output data:
% * RemovePow identifies the (logarithm of the) first SV were the system 
% must be trucated; and,
% * entry is used to calculate the number of singular values that were
% eliminated.
[RemovePow,entry] = FindSvdOutliers(LHS,'fit'); %Choose 'fit' or 'jump'

% If SV outliers are found, two situations may occur:
% * if a single analysis is performed (no p-adaptivity), then a warning is
% issued and the solution computed using the Moore-Penrose pseudoinverse of
% the coefficients matrix, truncated at 10^RemovePow; or,
% * if a p-adaptive analysis is performed, then serious instability may be
% caused in subsequent iterations, so the run ends with an error.
if ~isempty(RemovePow)
    
    if RunOption==1 % system is unstable, but a single run is required
        warning('local:NumericalChk',...
            'Warning: %d SVD outliers found. Results may be incorrect. ',...
            Dim-entry+1);
        % Moore-Penrose pseudoinverse-based solution procedure
        X = pinv(LHS,10^RemovePow)*RHS;
    else   % system is unstable, and a p-adaptive run is required
        error('local:NumericalChk',...
            'System has %d SVD outliers. SVD outliers are not accepted in the first run of adaptive analyses. Exiting. ',Dim-entry+1);
    end
    
    % If no SV outliers are found,
    % * the solving system is pre-conditioned using a scaling procedure aimed
    % at reducing its diagonal elements to unity while preserving the symmetry
    % of the original system; and,
    % *the default Matlab solver (mldivide) is used to obtain the solution.
else
    % System scaling procedure:
    % Generating the scaling matrix, Sc. Sc is a diagonal matrix, whose terms
    % are defined as the inverse of the square roots of the diagonal terms of
    % the coefficient matrix.
    Sc = sqrt(diag(LHS)).^-1;
    Sc(isinf(Sc)) = 1;
    Sc = diag(Sc);
    % Scaling LHS and RHS. ScLHS and ScRHS are the scaled versions
    % of the matrix of coefficients and the free vector.
    ScLHS = Sc' * LHS * Sc;
    ScRHS = Sc' * RHS;
    % Solving the scaled system by using default solution procedure
    ScX = ScLHS\ScRHS;
    % Scaling back the solution
    X = Sc * ScX;
end

%% Post-processing for the single run analysis
%
% If a single analysis is required, construct the temperature and flux 
% fields and store them in the Loops structure and a TecPlot file.
% Subsequently, it plots the temperature and heat flux fields.
if(RunOption==1)
    % 'iteration' number in the file name is set to zero to indicate a
    % single run.
    iteration = 0;
    % Construct the temperature and flux fields and plot the solution.
    % COMPUTEFIELDSTRI returns the Loops structure with the values of the
    % fields in the Gauss points of each element.
    [Loops] = ComputeFieldsTri(Nodes, Loops, X, iteration, NGP);
    
    % The solution energy is computed according to Section 3.9 of of 
    % reference [4]. In single run analyses, the solution energy is not
    % directly used.
    Ep = Energy_p(Loops, Edges, X, abscissa, weight);
    Ec = Energy_ctTri(Nodes,Loops,NGP);
    Energy = -(1/2)*(X(1:Edges.Edgeinsert-1))'*...
        LHS(1:Edges.Edgeinsert-1,1:Edges.Edgeinsert-1)*...
        (X(1:Edges.Edgeinsert-1))-Ep+Ec; 
    display(Energy); 
   
    % PLOTFIELDSTRI plots the temperature and flux fields
    PlotFieldsTri(Nodes, Edges, Loops,NGP);
        
    % single run - stops execution and returns
    return; 
end
%% ******** END OF SINGLE RUN ANALYSIS ********





%% ******** ADAPTIVE ALGORITHM INITIALIZATION *******
% The algorithm only gets to this point is an adaptive analysis was
% requested. The initial analysis is completed at this point.

% Saving the coefficients matrix from the initial run to avoid its
% recalculation
LHS0=LHS; 

% Storing the initial refinements. Since no SV outliers were identified in
% the first run, the solving system is assumed stable. Edges.stable and 
% Loops.stable fields of the Edges and Loops structures store the last
% orders of the boundary and domain bases for which no outliers were found.
Edges.stable = Edges.order; 
Loops.stable = Loops.order;

% Computing the initial solution energy, according to Section 3.9 of of 
% reference [4]. The constant part of the solution energy, associated to
% the particular solution, is omitted, since it does not change over the
% iterative refinement process and is thus irrelevant to its convergence. 
Ep = Energy_p(Loops, Edges, X, abscissa, weight);
Energy0 = -(1/2)*(X(1:Edges.Edgeinsert-1))'*...
    LHS(1:Edges.Edgeinsert-1,1:Edges.Edgeinsert-1)*...
    (X(1:Edges.Edgeinsert-1))-Ep;

% EdgesDirichlet is a list with all Dirichlet edges
EdgesDirichlet=strfind(Edges.type','D');

% Initialization
% * ErrorNorm is the value of the error associated to the current stopping
% criterion (see Section 6.3.3 of reference [2]). Depending on the
% stopping criterion, it may represent the average energy variation, 
% computed according over the last AvgNVal iterations (energy convergence
% stopping criterion), or the average boundary residual, also computed over
% the last AvgNVal iterations;
% * iteration is the current iteration;
ErrorNorm = 1;
iteration = 1;

% LIST is a structure where the algorithmic data of the iterative process
% is stored. It consists of two types of variables:
% * variables with names ending in 'It', which are typically lists that
% function as a logbook of the iterative process; and,
% * variables without 'It' which are local to (and overwritten at) each 
% iteration.
% The fields of the LIST structure are fully described in Section 6.2.2 of
% reference [4].
% List.EnergyIt is a vector with the solution energy at each iteration.
List.EnergyIt(iteration)=0;

% Load p-adaptive algortihmic parameters defined by the user in the GUI.
% The definiton of these parameters and their recommended values are given
% in Section 4.5.6 of reference [3].
load('Adaptive','TargetErrorNorm','SelectionTol','MaxOutlierIter',...
    'SelectionCriterion','thresh', 'MinIter','AvgNVal','StoppingCriterion');
% The SelectionCriterion value coming from the GUI must be incremented, as
% the first column in List.Edge is the index.
SelectionCriterion = SelectionCriterion + 1;


%% ******** ADAPTIVE ALGORITHM BODY ********
% This is the main body of the adaptive algorithm, ran iteratively until
% the convergence is reached (or some other stopping criterion is met).
% The functioning of the adaptive p-refinenent process is documented in
% Section 6.3 of reference [2], Sections 4.5.3 to 4.5.6 of reference [3]
% and in Chapters 5, 6 and 8 (and appendices) of reference [4].
%
% The adaptive refinement is performed by increasing the orders of the 
% approximation bases of some boundaries, selected based on either a global
% or a local selection criteria.
% The refinement of the domain bases is a consequence of the boundary
% refinement to ensure the indeterminacy and numerical stability of the
% solving system.

%% Begin iterative process, until the convergence (or some stopping criterion) are reached
while ErrorNorm > TargetErrorNorm % checking for convergence
    
    %% Initialization
    display(iteration); 
    
    % Initializing structures resposible for the edge refinement management
    % * List.Edge is an NBNDx3 matrix, where NBND is the total number of
    % boundaries of the mesh. The columns list the index of the boundary
    % and the values of the two selection criteria. The selection criteria 
    % are the boundary balance and energy residuals, computed as shown in
    % Section 6.3 of reference [2];
    % * List.SpuriousEdges lists the edges whose incrementation would cause
    % spurious (i.e. zero energy) modes;
    % * List.EdgesToRefine lists the edges to refine in the current 
    % iteration.
    List.Edge=[];List.SpuriousEdges=[];List.EdgesToRefine=[];
    
    %% Selection of the boundaries to refine
    % For each essential boundary, the values of the boundary balance and 
    % energy residuals are computed and stored in the List structure 
    % (List.Edges), ranked according to the (user) defined 
    % SelectionCriterion.
    for i=1:length(EdgesDirichlet)
        index=EdgesDirichlet(i);
        [Edges,Loops,List] = EdgeRefinement(Edges,Loops,BConds,...
            LHS0,Energy0,abscissa,weight,X,List,Dim,index,iteration,...
            List.EnergyIt,SelectionCriterion);
    end
    
    % Selects the edges to refine. The boundaries selected for refinement
    % are the boundary with the highest value of the selection criterion
    % (SC), and all boundaries for which the selection criteria are within
    % a specified selection tolerance (SelectionTol) from the highest
    % value. See Section 6.3.2 of reference [2].
    % ExitScore is zero if the selection crietria for all boundaries are
    % inferior to a 'numerical zero' threshold (thresh).
    [List,ExitScore] = SelectEdgesToRefine(List,SelectionCriterion,...
        SelectionTol,thresh);
    
    % If the selection criterion is smaller than the threshold for all 
    % boundaries, it exits the while cycle. The most probable cause is that
    % the exact solution is already contained in the basis, so refinement
    % increments would not bring any improvement.    
    if ~ExitScore   
        break;
    end
        
    % The approximation basis order of all edges marked for refinement are
    % increased by one.
    Edges.order(List.EdgesToRefine)=Edges.order(List.EdgesToRefine)+1;
    Edges.dim(List.EdgesToRefine)=Edges.dim(List.EdgesToRefine)+1;
    
    % The relevant iteration dependent information is stored for posterior 
    % convergence plotting.
    % * List.RefinedEdgesIt lists the edges refined in the current iteration 
    % * List.ErrorEdgeNormIt lists the largest boundary balance residual in
    % the current iteration.
    List.RefinedEdgesIt{iteration,:}=List.EdgesToRefine; 
    List.ErrorEdgeNormIt(iteration)=List.Edge(1,2); 

    
    %% Selection of the domains to refine
    % There are three possible causes for domain refinement: 
    % * the first is the need to refine a boundary whose refinement would 
    % yield spurious modes. This should not normally happen, but it might,
    % if the threshold (thresh) for calling a spurious mode is set too
    % high. If this happens, the elements adjacent to the boundaries to 
    % refine are also refined;
    % * the second is if the indeterminacy number of certain elements 
    % becomes equal or less than zero because of the refinement of adjacent
    % boundaries;
    % * the third is if after the refinement of the edges, the order of an 
    % edge becomes larger then the order of an element adjacent to it. 
    % Then the domain order must also be increased.
    
    % Identify boundaries with spurious modes selected for refinement
    List.SpurEdgeToRefine = intersect(List.EdgesToRefine,List.SpuriousEdges);

    % Initialization of the element (loop) refinement procedure.
    % * List.LoopsToRefine is a vector that contains the elements to refine
    % * RefinedLoopsIt store the elements refined at each iteration.
    List.LoopsToRefine = [];
    List.RefinedLoopsIt{iteration,:} = [];
      
    % Identify the elements (loops) to be refined
    List = SelectLoopsToRefine(Loops,Edges,List,iteration);
    
    % The orders of the approximation bases order of all loops marked for 
    % refinement are increased by one
    Loops.order(List.LoopsToRefine)=Loops.order(List.LoopsToRefine)+1;    

        
    %% Rebuilding and solving the finite element system
    % Rebuild the system according to the newly selected bases orders.
    % If multiple boundaries are selected for refinement, it is possible
    % that the system becomes ill conditioned even if none of the added
    % modes would be spurious if selected on its own. 
    % The numerical instability criterion is checked by searching for the 
    % SV outliers. If outliers are encountered, the order of the elements 
    % adjacent to the refined boundaries are automatically increased to 
    % improve the stability of the system.
    % If SV outliers are present after MaxOutlierIter several attempts, 
    % the algorithm quits. However, the program returns the lsat stable
    % solution.
    
    
    % Initializing parameters
    % * Outliers is a flag signaling the presence of SV outliers;
    % * RevertAndExit is a flag signaling that MaxOutlierIter iterations
    % with SV outliers were performed. The last stable solution is returned
    % in this case;
    % * OutIt is the number of iterations with SV outliers.
    Outliers= 1;
    RevertAndExit = 0;
    OutIt=1;    

    % Begin rebuilding, while checking for outliers in the system.
    while  Outliers && ~RevertAndExit % runs while outliers are detected
        
        % Mapping the solving system. Please see the comments to the 
        % ASSIGNPARTS function.
        [Edges,Loops,Dim] = AssignParts(Edges, Loops);

        % The system construction procedure is identical to that of the
        % first iteration.
        Edges.Edgeinsert=Loops.insert(end)+Loops.dim(end); 
        LHS = zeros(Dim);
        RHS = zeros(Dim,1);
        LHS = Gen_K_Matrix(Edges,Loops, LHS, abscissa, weight);
        LHS = Gen_B_Matrix(Edges, Loops, LHS, abscissa, weight);
        RHS = Gen_qg_Vector(Edges, Loops, BConds, RHS, abscissa, weight);
        RHS = Gen_tg_Vector(Edges, BConds, RHS, abscissa, weight);
        RHS = Gen_Kp_Vector(Edges, Loops, RHS, abscissa, weight);
        RHS = Gen_tp_Vector(Edges,Loops, RHS, abscissa, weight);
                
        % Check if the maximum orders of the domain and boundary bases 
        % were reached.
        [ErrorNorm,Stop]=CheckMaxOrders(Loops,Edges,ErrorNorm);       
        
        % Checking for outliers in the system
        [RemovePow,~] = FindSvdOutliers(LHS,'fit'); %Choose 'fit' or 'jump';
        
        if isempty(RemovePow)
            % if there's no 'RemovePow' it means no outlier was detected, 
            % and if the maximum degree was not reached,
            % then the execution proceeds to the next iteration.                   
            Outliers=0;
            OutIt=0;              

        else % there are outliers, the loops adjacent to the refined 
             % boundaries will aslo be refined
             
            % Reset the list with the loops to refine
            List.LoopsToRefine = [];
            
            % Identify the elements (loops) which are adjacent to the
            % boundaries selected for refinement.
            for i=1:length(List.EdgesToRefine)
                index = List.EdgesToRefine(i);
                List.LoopsToRefine = cat(2,List.LoopsToRefine,...  
                    Edges.lleft(index),Edges.lright(index));       
            end
            
            % Eliminating the repeated loops and the zero entries (the
            % latter correspond to right elements of exterior boundaries.
            List.LoopsToRefine = unique(List.LoopsToRefine);
            List.LoopsToRefine = List.LoopsToRefine(List.LoopsToRefine~=0); 
            
            % Add newly refined elements to the list of those refined in
            % the current iteration
            List.RefinedLoopsIt{iteration,:} = ...
                cat(2, List.RefinedLoopsIt{iteration,:},...
                List.LoopsToRefine);
            
            % Increment the orders of the elements to be refined
            Loops.order(List.LoopsToRefine)=Loops.order(List.LoopsToRefine)+1;
                        
            % Setting flags to be able to diagnose forced exits
            % * If the max no of iterations with outliers was reached, sets
            % RevertAndExit to 2;
            % * If the maximum order of a basis was reached, (Stop == 1),
            % sets RevertAndExit to 1;
            if OutIt==MaxOutlierIter
                RevertAndExit=2;
                Edges.order=Edges.stable;  
                Loops.order=Loops.stable;
            elseif Stop==1 
                RevertAndExit=1;
                Edges.order=Edges.stable;  
                Loops.order=Loops.stable;
            end
            
            OutIt=OutIt+1;
            
        end
        
    end
    
    
    %% Dealing with revert & exit scenarios
    % At this point, we have the following scenarios:
    % - Stop = 0 & RevertAndExit = 0: solve the system & proceed
    % - Stop = 1 & RevertAndExit = 0: max refinement reached and no outliers - solve the system & quit
    % - Stop = 1 & RevertAndExit = 1: max refinement reached and outliers - revert & quit
    % - Stop = 1 & RevertAndExit = 2: max refinement reached and max outlier iterations reached - revert & quit
       
    %% Dealing with revert & exit scenarios
    if RevertAndExit
        % system needs to be remapped for the last stable configuration
        [Edges,Loops] = AssignParts(Edges, Loops); 
        
        % If the orders had to be reverted to the last stable configuration
        % it registers 'Reverted' in the iteration logbooks.
        List.RefinedLoopsIt{iteration,:} = 'Reverted';
        List.RefinedEdgesIt{iteration,:} = 'Reverted';
        
        % Storing iteration dependent stuff for posterior convergence plotting
        % removes the last edge error norm from the list (reverted iteration)
        List.ErrorEdgeNormIt=  List.ErrorEdgeNormIt(1:iteration-1); 
        if RevertAndExit == 1
            warning('The maximum order of refinement was reached and SVD outliers were reported. Reverting to the last stable configuration...');
        else
            warning('The maximum number of consecutive systems with SVD outliers was reached. Reverting to the last stable configuration...');
        end
        
        % the other operations in the while ErrorNorm > TargetErrorNorm 
        % loop are not performed for RevertAndExit~=0
        break; 
    end
    
    
    %% Stable configuration protocol
    % The refinement process of the current iteration is sucessfuly 
    % completed and the relevenat information is stored. 
    Edges.stable = Edges.order;
    Loops.stable = Loops.order;
    LHS0=LHS;
    
    % Preconditioning and solving the system. The system procedure is
    % identical to that of the first iteration.
    Sc = sqrt(diag(LHS)).^-1;
    Sc(isinf(Sc)) = 1;
    Sc = diag(Sc);
    % Scaling LHS and RHS. ScLHS and ScRHS are the scaled versions
    % of the matrix of coefficients and the free vector.
    ScLHS = Sc' * LHS * Sc;
    ScRHS = Sc' * RHS;
    % Solving the scaled system by using default solution procedure
    ScX = ScLHS\ScRHS;
    % Scaling back the solution
    X = Sc * ScX;
    
    % Computes the new solution energy and stores it for the next iteration
    Ep = Energy_p(Loops, Edges, X, abscissa, weight);
    List.EnergyIt(iteration)= -(1/2)*(X(1:Edges.Edgeinsert-1))'*...
        LHS(1:Edges.Edgeinsert-1,1:Edges.Edgeinsert-1)*...
        (X(1:Edges.Edgeinsert-1))-Ep;
    
    % Calls StoreIterationInfo to insert the information about the current
    % iteration into the List structure
    List = StoreIterationInfo(Loops,Edges,List,iteration,Energy0);
    
    % Computes the values of the stopping criteria.
    % Computes the average error over the last AvgNVal+1 iterations
    % to check if it is below the convergence criterion,TargetErrorNorm.
    % At least MinIter are always performed, to avoid spurious early
    % convergence.
    if StoppingCriterion == 1 && iteration>=MinIter
        ErrorNorm=sum(List.EnergyVariationIt((iteration)-AvgNVal:...
            (iteration)))/(AvgNVal+1);
    elseif StoppingCriterion == 2 && iteration>=MinIter
        ErrorNorm=sum(List.ErrorEdgeVariationIt((iteration)-AvgNVal:...
            (iteration)))/(AvgNVal+1);
    end
    
    % If the maximum order was reached, it exists 
    if Stop == 1
        load('Adaptive','MaxOrder');
        for j=1:length(Loops.area)
            if (Loops.order(j) >= MaxOrder)
                warning('The maximum order of refinement was reached for element %d. Try increasing the energy error set for convergece.',j);
            end
        end
        ErrorNorm = 0.0;  
    end
    
    % Updating the iteration number
    iteration=iteration+1;
end

%% ******** POSTPROCESSING ********
% Constructs the final temperature and flux fields, the final order map and
% the convergence graphs.

% After the end of the iterative process, updates the iteration number
iteration=iteration-1;

% Computing the temperature and flux fields and storing them in the
% Loops structure and a TecPlot file.
[Loops] = ComputeFieldsTri(Nodes, Loops, X, iteration, NGP);

% Plotting the temperature and flux fields
PlotFieldsTri(Nodes, Edges, Loops,NGP);

% If no iterations have been performed, due to the initial
% solution having been exact already, it skips the rest. Otherwise, it
% constructs the convergence graphs.
if iteration ~= 0
    
    % Constructing the convergence graphs 
    GDLTotal=cat(2,GDL0,List.GDL_It);
    ConvGraphs(GDLTotal,List.EnergyIt,List.ErrorEdgeNormIt,...
        List.EnergyVariationIt, List.ErrorEdgeVariationIt,Energy0); 
    
end

end

