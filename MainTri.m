function MainTri
%Main processor program, hybrid elements, regular mesh
tic;

%
%% ******** PRE-PROCESSING ********
% Processes the input data, generates the mesh and topological data
% structures
[NGP, Nodes, Edges, Loops, BConds] = InputProcTri;
[Loops,Edges] = CheckMinDegrees(Loops,Edges); % checks and corrects the initial indeterminacy

load('HeatStructDef','RunOption');

%% ******** INITIAL ANALYSIS ********
% This is the first analysis of an adaptive process or the single analysis
% if no adaptive analysis is required

%%
% Assign parts to elements and sizes in the global system
[Edges,Loops,Dim] = AssignParts(Edges, Loops); % Obtaining the insertion points and dimensions of sub-matrices
Edges.Edgeinsert=Loops.insert(end)+Loops.dim(end); % stores the coordinate of the beginning of the edges' part in the global system
LHS = zeros(Dim);
RHS = zeros(Dim,1);

%%
% Initialization of Gauss weights & abscissas (on a -1:1 interval)
[abscissa, weight] = gauleg(NGP, -1, 1);

%%
% Generating & allocating the LHS matrices
LHS = Gen_K_Matrix(Edges, Loops, LHS, abscissa, weight);
LHS = Gen_B_Matrix(Edges, Loops, LHS, abscissa, weight);
% SparsityIndex = (1-(nnz(LHS)/(length(LHS))^2))*100;
GDL0 = sum(Edges.dim(:))+sum(Loops.dim(:)); % Initial nº of DOF

%%
% Generating & allocating the RHS vectors
RHS = Gen_tg_Vector(Edges, BConds, RHS, abscissa, weight);
RHS = Gen_qg_Vector(Edges, Loops,BConds, RHS, abscissa, weight);
RHS = Gen_Kp_Vector(Edges, Loops, RHS, abscissa, weight);
RHS = Gen_tp_Vector(Edges,Loops, BConds, RHS, abscissa, weight);

%%
% Solving the FE system

% First step is checking the singular values of LHS, looking for sudden
% decreases (outliers) that may indicate the system is ill-conditioned or 
% singular. 

%If a single run is required, the system can be singular and the solution
%is computed using the 'pseudo-inverse' method.
%If an adaptative analyses is required, the system can't be singular.
%There's a third possibility, the system isn't sigular, but SVD outliers
%were detected. Then, the 'pseudo-inverse' method is used to
%computate the solution.

[RemovePow,entry] = FindSvdOutliers(LHS,'fit'); %Choose 'fit' or 'jump'

if ~isempty(RemovePow) % if SVD outliers were found
    % checking if the first run system is truly singular
    % cleaning the LHS
    LHCln = LHS; LHCln(abs(LHCln)<eps) = 0;
    
    if ~rcond(LHCln) && RunOption==1 % system is singular, but a single run is required
        cont = input('System is singular. Do you wish to continue (y/n)?','s');
        if ~strcmpi(cont,'y')
            return;
        else
            X = pinv(LHS,10^RemovePow)*RHS;
        end
        
    elseif RunOption==2 % if the first system is singular,
        % an adaptive analysis cannot be performed
        error('local:NumericalChk',...
            'System has %d SVD outliers. SVD outliers are not accepted in the first run of adaptive analyses. Exiting. ',Dim-entry+1);
        
    else
        % system is not singular, but SVD outliers were found
        warning('local:NumericalChk',...
            'Warning: %d SVD outliers found. Results may be incorrect. ',...
            Dim-entry+1);
        X = pinv(LHS,10^RemovePow)*RHS;
    end
    
else
    % System scaling procedure
    % computing the scaling factors
    Sc = sqrt(diag(LHS)).^-1;
    Sc(isinf(Sc)) = 1;
    Sc = diag(Sc);
    % scaling LHS and RHS
    ScLHS = Sc' * LHS * Sc;
    ScRHS = Sc' * RHS;
    % solving the scaled system
    ScX = ScLHS\ScRHS;
    % scaling back
    X = Sc * ScX;
end

%% ******** SINGLE ANALYSIS OPTION ********
% If a single analysis is required, it just computes the final results,
% prints the fields and quits.

if(RunOption==1) % if a single run is required
       
    %%
    % Computing the temperature and flux fields and storing them in the
    % Loops structure and a TecPlot file. To illustrate that a single run
    % was made, the 'iteration' number in the file name is zero.
    iteration = 0;
    [Loops] = ComputeFieldsTri(Nodes, Loops, X, iteration, NGP);
    
    %%
    % Computing the thermal energy
    Ep = Energy_p(Loops, Edges, X, abscissa, weight);
    Ec = Energy_ctTri(Nodes,Loops,NGP);
    
    Energy = -(1/2)*(X(1:Edges.Edgeinsert-1))'*...
        LHS(1:Edges.Edgeinsert-1,1:Edges.Edgeinsert-1)*...
        (X(1:Edges.Edgeinsert-1))-Ep+Ec; %For elements with a rectangular shape
    display(Energy); 
    
    toc;
    
    %%
    % Plotting the temperature and flux fields
    PlotFieldsTri(Nodes, Edges, Loops,NGP);
        
    return; % single run - stops execution and returns
end

%% ******** ADAPTIVE ALGORITHM INITIALIZATION ********
% Storing the initial data, coming from the first execution, loading
% adaptive algorithm control data and initializing iteration parameters.

%%
% Storing the initial system and the respective (stable) refinements
LHS0=LHS; % Saving the initial left hand side system for later reconstruction
Edges.stable = Edges.order; % these are  the last refinements where the system had no zero pivots
Loops.stable = Loops.order;

%%
% Computing the initial Energy (constant part omitted)
Ep = Energy_p(Loops, Edges, X, abscissa, weight);
Energy0 = -(1/2)*(X(1:Edges.Edgeinsert-1))'*...
    LHS(1:Edges.Edgeinsert-1,1:Edges.Edgeinsert-1)*...
    (X(1:Edges.Edgeinsert-1))-Ep;

%%
% Listing all Dirichlet edges 
EdgesDirichlet=strfind(Edges.type','D');

%%
% Initialize & load parameters
ErrorNorm = 1;
iteration = 1;
List.EnergyIt(iteration)=0;
load('Adaptive','TargetErrorNorm','SelectionTol','MaxOutlierIter',...
    'SelectionCriterion','thresh', 'MinIter','AvgNVal','StoppingCriterion');
% The SelectionCriterion value coming from the GUI must be incremented, as
% the first column in List.Edge is the index
SelectionCriterion = SelectionCriterion + 1;




%% ******** ADAPTIVE ALGORITHM BODY ********
% This is the main body of the adaptive algorithm, ran iteratively until
% the convergence criterion is reached (or some other exit criterion is
% met).

while ErrorNorm > TargetErrorNorm % checking for convergence
    
    display(iteration)
    
    %%
    % ------ SELECTION OF THE EDGE TO BE REFINED ---------
    
    %%
    % Initializing edge refinement
    List.Edge=[];List.SpuriousEdges=[];List.EdgesToRefine=[];
    
    %%
    % For each essential edge, launching the function that computes the
    % values of all selection criteria. These values are stored in the
    % List structure (List.EdgesToRefine), ranked according to the desired
    % SelectionCriterion.
    for i=1:length(EdgesDirichlet)
        index=EdgesDirichlet(i);
        [Edges,Loops,List] = EdgeRefinement(Edges,Loops,BConds,...
            LHS0,Energy0,abscissa,weight,X,List,Dim,index,iteration,...
            List.EnergyIt,SelectionCriterion);
    end
    
    %%
    % Selects the edges that will be refined
    [List,ExitScore] = SelectEdgesToRefine(List,SelectionCriterion,SelectionTol,thresh);
    
    if ~ExitScore   % the selection criterion is smaller than the threshold for all boundaries
        break;
    end
    
    
    %%
    % Refining ALL edges within the tolerance
    Edges.order(List.EdgesToRefine)=Edges.order(List.EdgesToRefine)+1;
    Edges.dim(List.EdgesToRefine)=Edges.dim(List.EdgesToRefine)+1;
    
    %%
    % Storing iteration dependent stuff for posterior convergence plotting
    List.RefinedEdgesIt{iteration,:}=List.EdgesToRefine; % stores the edges refined in the current iteration
    List.ErrorEdgeNormIt(iteration)=List.Edge(1,2); % stores the edge error norm corresponding to the first-ranked edge
    %%
    % -------- END OF EDGE REFINEMENT ---------
    
    
    
    %%
    % ----------- BEGIN DOMAIN REFINEMENT -----------
    
    %%
    % There are three possible causes for domain refinement: the first is the
    % need to refine a boundary whose refinement would yield spurious
    % modes. Then, the adjacent elements also need to be refined.
    % Otherwise, we need to refine a domain whose indeterminacy number is
    % equal or less than zero after the refinement of its boundaries.
    %The third cause is after refining the edges, if its final order is larger then
    %the order of the domain adjacent to it, then the domain order must be
    %increased.
    
    %%
    % Identify boundaries with spurious modes selected for refinement
    List.SpurEdgeToRefine = intersect(List.EdgesToRefine,List.SpuriousEdges);
    %%
    % Initialization of loops refinement structure
    List.LoopsToRefine = [];
    List.RefinedLoopsIt{iteration,:} = [];
      
    %%
    % Identifies the elements to be refined, namely loops with boundaries
    % selected for refinement that would cause spurious modes, loops with
    % non-positive kinematic indeterminacy numbers and loops with boundaries
    % of larger refinement than the loops themselves
    List = SelectLoopsToRefine(Loops,Edges,List,iteration);
    
    %%
    % increment the orders of the elements to be refined
    Loops.order(List.LoopsToRefine)=Loops.order(List.LoopsToRefine)+1;
    
    % ----------------- END OF DOMAIN REFINEMENT ---------------- %
    
    
    %%
    %-------------------- SOLVING SYSTEM OPERATIONS ------------------ %
    % Rebuilds the system according to the newly selected refinements.
    % If multiple boundaries are selected for refinement, it is possible
    % that the system becomes ill conditioned even if none of the added
    % modes would be spurious if selected on its own. In order to check for
    % this situation, the conditioning of the p-refined system cannot have
    % SVD outliers, that are harmful to the solution, detected in the sistem.
    %If any outlier is detected, the order of all elements adjacent to
    %the refined boundries are increased until the outliers stop being detected.
    %However, if 'MaxOutlierIter' is reached, the execution ends and returns
    %to the last stable solution.
    
    %Initializing parameters
    Outliers= 1;
    RevertAndExit = 0;
    OutIt=1;
    
    %%
    % Begin rebuilding, while checking for outliers in the sistem
    %to identify ill-conditioning.
    
    while  Outliers && ~RevertAndExit % runs while outliers are detected
        
        %%
        % Constructing the new solving system and updating Edgeinsert
        [Edges,Loops,Dim] = AssignParts(Edges, Loops);
        Edges.Edgeinsert=Loops.insert(end)+Loops.dim(end); % stores the coordinate of the beginning of the edges' part in the global system
        LHS = zeros(Dim);
        RHS = zeros(Dim,1);
        LHS = Gen_K_Matrix(Edges,Loops, LHS, abscissa, weight);
        LHS = Gen_B_Matrix(Edges, Loops, LHS, abscissa, weight);
        RHS = Gen_qg_Vector(Edges, Loops, BConds, RHS, abscissa, weight);
        RHS = Gen_tg_Vector(Edges, BConds, RHS, abscissa, weight);
        RHS = Gen_Kp_Vector(Edges, Loops, RHS, abscissa, weight);
        RHS = Gen_tp_Vector(Edges,Loops, BConds, RHS, abscissa, weight);
                
        % Check if the maximum degrees were reached
        [ErrorNorm,Stop]=CheckMaxDegrees(Loops,Edges,ErrorNorm);       
        
        %%
        % Checking for outliers in the system
        [RemovePow,~] = FindSvdOutliers(LHS,'fit'); %Choose 'fit' or 'jump';
        
        if isempty(RemovePow)
            
            Outliers=0;             %if there's no 'RemovePow' it means no outlier was detected, and if the maximum degree was not reached,
            OutIt=0;                %then the execution proceeds to the next iteration.
            
        else % there are outliers
            List.LoopsToRefine = [];
            for i=1:length(List.EdgesToRefine)
                index = List.EdgesToRefine(i);
                List.LoopsToRefine = cat(2,List.LoopsToRefine,...  % if an outlier is detected, all elements adjacent to the refined boundries
                    Edges.lleft(index),Edges.lright(index));       % are refined
            end
            List.LoopsToRefine = unique(List.LoopsToRefine);
            List.LoopsToRefine = List.LoopsToRefine(List.LoopsToRefine~=0); % eliminate the zero entries (no right element)
            
            %%
            % Add newly refined elements to the list of those refined
            % because of the other possible reasons in the same
            % iteration
            List.RefinedLoopsIt{iteration,:} = ...
                cat(2, List.RefinedLoopsIt{iteration,:},...
                List.LoopsToRefine);
            
            %%
            % Increment the orders of the elements to be refined
            Loops.order(List.LoopsToRefine)=Loops.order(List.LoopsToRefine)+1;
                        
            if OutIt==MaxOutlierIter
                RevertAndExit=2;
                Edges.order=Edges.stable;  %If 'MaxOutlierIter' is reached it returns to the last stable solution
                Loops.order=Loops.stable;
            elseif Stop==1 % if there are SVD outliers and the max degree was reached, just revert & exit
                RevertAndExit=1;
                Edges.order=Edges.stable;  %If 'MaxOutlierIter' and 'MaxDegrees' are reached it returns to the last stable solution
                Loops.order=Loops.stable;
            end
            
            OutIt=OutIt+1;
            
        end
        
    end
    % At this point, we have the following scenarios:
    % - Stop = 0 & RevertAndExit = 0: all fine, solve the system & proceed
    % - Stop = 1 & RevertAndExit = 0: max refinement reached and no outliers - solve the system & quit
    % - Stop = 1 & RevertAndExit = 1: max refinement reached and outliers - revert & quit
    % - Stop = 1 & RevertAndExit = 2: max refinement reached and max outlier iterations reached - revert & quit
    
    
    %%
    % ------ DEALING WITH REVERT & EXIT SCENARIOS ------
    % If the orders had to be reverted to the last stable configuration, it
    % registers 'Reverted'
    if RevertAndExit
        [Edges,Loops] = AssignParts(Edges, Loops); % parts need to be reassigned in the last stable configuration
        List.RefinedLoopsIt{iteration,:} = 'Reverted';
        List.RefinedEdgesIt{iteration,:} = 'Reverted';
        %%
        % Storing iteration dependent stuff for posterior convergence plotting
        List.ErrorEdgeNormIt=  List.ErrorEdgeNormIt(1:iteration-1); % removes the last edge error norm from the list (reverted iteration)
        if RevertAndExit == 1
            warning('The maximum order of refinement was reached and SVD outliers were reported. Reverting to the last stable configuration...');
        else
            warning('The maximum number of consecutive systems with SVD outliers was reached. Reverting to the last stable configuration...');
        end
        break; % all of the others operations in the while ErrorNorm > TargetErrorNorm loop are not performed for RevertAndExit~=0
    end
    
    
    %%
    % ------ STABLE CONFIGURATION PROTOCOL ------
    %%
    % Stores stuff
    Edges.stable = Edges.order;
    Loops.stable = Loops.order;
    LHS0=LHS;
    
    %%
    % Solves for X
    % System scaling procedure
    % computing the scaling factors
    Sc = sqrt(diag(LHS)).^-1;
    Sc(isinf(Sc)) = 1;
    Sc = diag(Sc);
    % scaling LHS and RHS
    ScLHS = Sc' * LHS * Sc;
    ScRHS = Sc' * RHS;
    % solving the scaled system
    ScX = ScLHS\ScRHS;
    % scaling back
    X = Sc * ScX;
    
    %%
    % Computes the new energy and stores it for the next iteration
    Ep = Energy_p(Loops, Edges, X, abscissa, weight);
    List.EnergyIt(iteration)= -(1/2)*(X(1:Edges.Edgeinsert-1))'*...
        LHS(1:Edges.Edgeinsert-1,1:Edges.Edgeinsert-1)*(X(1:Edges.Edgeinsert-1))-Ep;
    
    %%
    % Calls StoreIterationInfo to do the listings
    [Loops,Edges,List]=StoreIterationInfo(Loops,Edges,List,...
        iteration,Energy0);
    
    %%
    % Computes the values of the stopping criteria

    % Computes the average the error over the last AvgNVal+1 iterations
    % to check if it is below the convergence criterion,TargetErrorNorm
    if StoppingCriterion == 1 && iteration>=MinIter
        ErrorNorm=sum(List.EnergyVariationIt((iteration)-AvgNVal:...
            (iteration)))/(AvgNVal+1);
    elseif StoppingCriterion == 2 && iteration>=MinIter
        ErrorNorm=sum(List.ErrorEdgeVariationIt((iteration)-AvgNVal:...
            (iteration)))/(AvgNVal+1);
    end
    
    % If the maximum degree was reached, it exists anyway
    if Stop == 1
        load('Adaptive','MaxOrder');
        for j=1:length(Loops.area)
            if (Loops.order(j) >= MaxOrder)
                warning('The maximum order of refinement was reached for element %d. Try increasing the energy error set for convergece.',j);
            end
        end
        ErrorNorm = 0.0;  % to make sure it quits
    end
    
    %%
    % Updating the current iteration
    iteration=iteration+1;
end

%%
% After the end of the iterative process, updates the iteration number
iteration=iteration-1;


%% ******** POSTPROCESSING ********
% Constructing temperature and flux fields, final order map and
% convergence graphs

%%
% ------ Constructing the final order map and the temperature and flux fields ------

%%
% Computing the temperature and flux fields and storing them in the
% Loops structure and a TecPlot file.
[Loops] = ComputeFieldsTri(Nodes, Loops, X, iteration, NGP);

toc;

%%
% Plotting the temperature and flux fields
PlotFieldsTri(Nodes, Edges, Loops,NGP);

%%
% If no iterations have been performed, for instance, due to the initial
% solution having been exact already, it skips the rest. Otherwise, it
% constructs the iteration information table and the convergence graphs.
if iteration ~= 0
    % Constructing the Table Iterations Info
%     TABLE(List.GDL_It,List.BetaIt,List.EnergyIt,List.LoopsOrderIt,...
%         List.EdgesOrderIt,List.RefinedEdgesIt,List.RefinedLoopsIt);
    
    % Constructing the data graphs showing convergence
    GDLTotal=cat(2,GDL0,List.GDL_It);
    GRAPHS(GDLTotal,List.EnergyIt,List.ErrorEdgeNormIt,...
        List.EnergyVariationIt, List.ErrorEdgeVariationIt,Energy0); % Use GDL instead of GDLTotal
    
end

end

