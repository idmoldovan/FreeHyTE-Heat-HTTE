function [Edges,Loops,List] = EdgeRefinement(Edges,...
    Loops,BConds,LHS0,Energy0,abscissa,weight,X,List,Dim,index,...
    iteration,EnergyIteration,SelectionCriterion)

% Abar is the old LHS matrix (LHS0) with the new B* corresponding to the
% newly refined boundaries (the LHS of equation 4.9)
Abar = Gen_Bbar_Matrix(LHS0,Edges,Loops,Dim,abscissa,weight,index);

% Computes the t_{\Gamma} of equation 4.9. If the edge is incremented with
% a single dof, it's a scalar, but the routine is constructed to
% accommodate its possible extension to a vectorial quantity.
tgbar = Gen_tgbar_Vector(Edges,BConds,abscissa,weight,index);
tpbar = Gen_tpbar_Vector(Edges,Loops,abscissa,weight,index);
tgbar = tgbar - tpbar;

% Compute xdot according to expression (4.13)
if (rcond(LHS0)<eps)
    xdot = -pinv(LHS0)*Abar(1:Dim,Dim+1);
else
    xdot = -LHS0\(Abar(1:Dim,Dim+1));
end


% Compute Error in the Edge according to equation (4.17)
ErrorEdge = (tgbar-Abar(Dim+1,1:Edges.Edgeinsert-1)*...
         X(1:Edges.Edgeinsert-1));     
     
% Compute Ybar according to expression 4.15 (is assumed scalar as one
% single dof is added to the side)
if abs(ErrorEdge) % if ErrorEdge ~= 0, (Abar(Dim+1,1:Dim)*xdot) CANNOT BE ZERO
    Ybar = ErrorEdge/(Abar(Dim+1,1:Dim)*xdot);
    % Error edge norm is normalized to the length of the edge
    ErrorEdgeNorm = abs(ErrorEdge)/(sqrt(Edges.parametric(index,3)^2 + Edges.parametric(index,4)^2));
else
    % if ErrorEdge = 0, the added dof brings no change to the solution
    List.Edge(index,1)=index;
    List.Edge(index,2)=0;
    List.Edge(index,3)=0;
    return;
end

% Compute Delta X "solution variation" according to equation (4.12)
DeltaX = xdot*Ybar;

% Compute the final energy, with the new variation
Ep = Energy_p(Loops, Edges, X+DeltaX , abscissa, weight);
Energy = -(1/2)*(X(1:Edges.Edgeinsert-1)+DeltaX(1:Edges.Edgeinsert-1))'*...
         Abar(1:Edges.Edgeinsert-1,1:Edges.Edgeinsert-1)*...
         (X(1:Edges.Edgeinsert-1)+DeltaX(1:Edges.Edgeinsert-1))-Ep;
% Compute the energy difference, 'DeltaU'
if iteration == 1
    DeltaU = abs(Energy0-Energy);
elseif iteration > 1
    DeltaU = abs(EnergyIteration(iteration-1)-Energy);
end
% List and tag all energy variation for each edge refinement, and tag loops
% for possible spurious modes and/or kinematic indetermenacy
List.Edge(index,1)=index;
List.Edge(index,2)=ErrorEdgeNorm;
List.Edge(index,3)=DeltaU;

% Sorting the list according to the selection criterion 
[List.Edge(:,SelectionCriterion),indice]=sort(List.Edge(:,SelectionCriterion),'descend'); 

% rewritting the reminder of the list according to the order determined above
for ii=1:3
    if ii ~= SelectionCriterion
        List.Edge(:,ii)=List.Edge(indice(:),ii);
    end
end
