function ConvGraphs(GDLTotal,EnergyIt,ErrorEdgeNormIt,...
    EnergyVariationNormIt,ErrorEdgeVariationNormIt,Energy0)
% CONVGRAPHS plots the variation of the values of the two selection
% criteria and the two convergence criteria over the course of the
% iterative process.
%
% CONVGRAPHS is called by MAIN***
%
% Input:
%  list with the total number of DOFs at each iteration (GDLTotal), the
%  four plotting quantities described below, and the solution energy in the
%  first iteration, Energy0
% Output:
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
% The four plots are briefly described next:
%  * Energy vs DOF: Plots the variation of the solution energy (EnergyIt) 
%  with the total number of degrees of freedom of the model. Each dot
%  corresponds to an iteration. The abscissa is logarithmic. The
%  stabilization of this quantity is an important indicator of the
%  convergence;
%  * Edge error density vs DOF: Plots the maximum boundary balance residual
%  (ErrorEdgeNormIt) at each iteration against the total number of degrees
%  of freedom of the model. Both axes are logarithmic. This quantity is
%  used as a selection criterion. The blue line in the plot is a linear fit
%  of the boundary residual decay and the function y(x) is its mathematical
%  expression;
%  * Normalized energy variation norm vs DOF: Plots the energy residual
%  (EnergyVariationNormIt) between two successive iterations against the
%  total number of degrees of freedom of the model. Both axes are
%  logarithmic. This quantity is used both as a refinement and as a
%  stopping criterion;
%  * Normalized edge error density vs DOF: Plots the maximum boundary
%  balance residual, normalized to its value in the first iteration
%  (ErrorEdgeVariationNormIt) against the total number of degrees of
%  freedom of the model. Both axes are logarithmic. This quantity is
%  used as a stopping criterion. 
%
% Further information on the refinement and stopping criteria can be found
% in reference [2] (Section 6.3.2), reference [3] (Section 4.5.4), and
% reference [4] (Sections 5.2 to 5.4, and 6.3).
%

%% Initialization
% catenating the information from the initial run
EnergyIt=cat(2,Energy0,EnergyIt);
ErrorEdgeNormIt=cat(2,0,ErrorEdgeNormIt);
EnergyVariationNormIt=cat(2,0,EnergyVariationNormIt);
ErrorEdgeVariationNormIt=cat(2,0,ErrorEdgeVariationNormIt);

%% First Plot: Energy vs DOF
Fig = figure('Position',[100 50 1200 650]);
set(Fig,'name','Graphs','numbertitle','off','color','w')
Fig1 = subplot(2,2,1);
plot(log10(GDLTotal),real(EnergyIt),'--rs','LineWidth',2,...
    'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',5);
title('Energy vs DOF','FontWeight','bold');
xlabel('(Degrees of freedom) log N');
ylabel('(Energy) U');
grid on;

%% Second Plot: Edge error density vs DOF
% Getting the linear fit of the boundary residual decay
p2=polyfit(log10(GDLTotal(length(ErrorEdgeNormIt)-nnz(ErrorEdgeNormIt)+1:end)),...
    log10(ErrorEdgeNormIt(length(ErrorEdgeNormIt)-nnz(ErrorEdgeNormIt)+1:end)),1);
p22=polyval(p2,log10(GDLTotal(length(ErrorEdgeNormIt)-nnz(ErrorEdgeNormIt)+1:end)));

% getting the anlytic expression of the linear fit
C2={'y','=',num2str(p2(1)),'x','+',num2str(p2(2))};
str2=sprintf('%s',C2{:});

% Plotting the data
Fig2=subplot(2,2,2);
plot(log10(GDLTotal),log10(ErrorEdgeNormIt),'--rs','LineWidth',2,'MarkerEdgeColor','k',...
    'MarkerFaceColor','g','MarkerSize',5);
hold on;
plot(log10(GDLTotal(length(ErrorEdgeNormIt)-nnz(ErrorEdgeNormIt)+1:end)),p22,'LineWidth',3);
text(0.6,0.95,str2,'FontWeight','bold','Units','normalized','BackgroundColor','w');
% text(0.6,0.85,str22,'FontWeight','bold','Units','normalized','BackgroundColor','w');
title('Edge error density vs DOF','FontWeight','bold');
xlabel('(Degrees of Freedom) log N');
ylabel('(Edge error density) log \epsilon_{\Gamma}');
grid on;

%%  Third Plot: Normalized energy variation norm vs DOF
Fig3=subplot(2,2,3);
plot(log10(GDLTotal),log10(EnergyVariationNormIt),'--rs','LineWidth',2,'MarkerEdgeColor','k',...
    'MarkerFaceColor','g','MarkerSize',5);
title('Normalized energy variation vs DOF','FontWeight','bold');
xlabel('(Degrees of freedom) log N');
ylabel('(Normalized energy variation) log \epsilon_{U}');
grid on;

%% Fourth Plot: Normalized edge error density vs DOF
Fig4=subplot(2,2,4);
plot(log10(GDLTotal),log10(ErrorEdgeVariationNormIt),'--rs','LineWidth',2,'MarkerEdgeColor','k',...
    'MarkerFaceColor','g','MarkerSize',5);
title('Normalized edge error density vs DOF','FontWeight','bold');
xlabel('(Degrees of freedom) log N');
ylabel('(Normalized edge error density) log \epsilon_{\Gamma} \cdot \epsilon_{\Gamma_0}^{-1}');
grid on;

end