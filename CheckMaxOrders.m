function [ErrorEnergyNorm,Stop]=CheckMaxOrders(Loops,Edges,ErrorEnergyNorm)
% CHECKMAXORDERS checks if the maximum order of refinement, defined by the 
% user in the GUI, was reached for any element or essential boundary. 
%
% CheckMaxDegrees is called by MAIN*** 
% Input:
%  Loops, Edges data structures and ErrorEnergyNorm. The maximum acceptable
%  order, MaxOrder, is also loaded from the GUI.
% Output/Returns to MAIN***
%  ErrorEnergyNorm (set to zero value if maximum order of the refinement
%  was reached, to abort the iterative process)
%  Stop code: 
%    is set to zero - if none of the elements/edges reached the maximum
%    order  
%    is set to "1" - if the maximum order was reached
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
% The excessive increase of the order of p-refinement may results in a loss
% of accuracy, because of computational instabilities. 
% To avoid the numerical problems the user is required to define in the GUI 
% the maximum allowable order for the bases. If this order is reached, the
% execution ends with "Maximum order criterion" and the most refined
% solution is output. 

%% Initialization
% Loading the maximum allowable order
load('Adaptive','MaxOrder');
Stop = 0;

%% Checking the maximum order criteria
% For the elements
for j=1:length(Loops.area)
    if (Loops.order(j) >= MaxOrder) 
          ErrorEnergyNorm=0.0;
          Stop=1;
    end
end

% For the boundaries
for i=1:length(Edges.type)
      if (Edges.order(i) >= MaxOrder)
          ErrorEnergyNorm=0.0;
          Stop=1;
      end
end

end

