function [List,ExitScore] = SelectEdgesToRefine(List,SelectionCriterion,...
    SelectionTol,thresh)
% SELECTEDGESTOREFINE selects the boundaries to be refined based on the
% ranking of the selection criterion and identifies the boundaries for
% which the addition of a new approximation function would not improve the
% solution (spurious modes). 
% 
% EdgeRefinement is called by MAIN*** 
% Input data:
%   the List structure, the SelectionCriterion, the selection tolerance
%   SelectionTol and the "numerical zero" threshold thresh, selected by the
%   user in the GUI.
% Output/Returns to MAIN***:
%  the updated List structure, updated with the essential boundaries to be
%  refined according to the selection criteria and with the boundaries with
%  spurious modes, and ExitScore flag, which is equal to zero if all
%  boundaries present spurious modes and one otherwise.
% Calls:
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
% All essential boundaries that present values of the selection criteria
% larger than SelectionTol times the maximum value of the selection
% criterion are listed for refinement in List.EdgesToRefine. 
% A boundary with spurious mode is a boundary where the addition of a new
% shape function would not contribute to the improvement of the solution.
% This happens when the free vector of system (5.10) in reference [4] is
% null (that is, inferior to the numerical zero threshold, thresh). These
% boundaries are listed in List.SpuriousEdges.
% If the exact solution of the problem is contained in the basis, the
% addition of new functions to any boundary would result in a spurious
% mode. If this happens, the ExitScore is zero and the iterative process is
% aborted in MAIN***.
%
% Further details on the refinement selection process can be found in
% reference [2] (Section 6.3.2), reference [3] (Section 4.5.4), and
% reference [4] (Sections 5.2 to 5.4, and 6.3).

%% Initialization

ExitScore = 1;

% The first-ranked element in the list is selected for refinement
List.EdgesToRefine = List.Edge(1,1);

%% Building the EdgesToRefine list
% Sweeps all other edges and checks if the values of their selection
% criteria stand within the defined selection tolerance from the
% first-ranked edge. If they do, it adds them to the EdgesToRefine list.
for i=2:length(List.Edge)
    if List.Edge(i,SelectionCriterion) >= ...
            List.Edge(1,SelectionCriterion)*SelectionTol
        List.EdgesToRefine = cat(2,List.EdgesToRefine,List.Edge(i,1));
    end
end

% Clear rows of Neumann Edges in List.Edge
List.Edge(~any(List.Edge,2),:)=[];

%% Check if all modes are spurious
% Analyses the values of the SelectionCriterion and if they are
% smaller than the threshold for all boundaries, it means that the exact
% solution was  recovered at the previous step. If this is the
% case, it returns ExitScore = 0 to MAIN***, to plot the solution and exit.
if max(List.Edge(:,SelectionCriterion)) < thresh
    warning('local:NumericalChk',...
        'The values of the selection criterion are smaller than the threshold %g for all boundaries. Suspect that the exact solution was already achieved.\n',thresh);
    ExitScore = 0;
end

%% Building the SpuriousEdges list
% Analyses the temperature balance improvements caused by the addition of
% new modes to the boundaries (column 2 in the List.Edge) and if they are
% smaller than a threshold, lists the boundary as spurious.
for index=1:size(List.Edge,1)
    if List.Edge(index,2)<thresh
        List.SpuriousEdges = cat(1,List.SpuriousEdges,List.Edge(index,1));
    end
end


end
