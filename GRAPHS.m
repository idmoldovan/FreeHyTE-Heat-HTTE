function GRAPHS(GDLTotal,EnergyIt,ErrorEdgeNormIt,...
    EnergyVariationNormIt,ErrorEdgeVariationNormIt,Energy0)

EnergyIt=cat(2,Energy0,EnergyIt);
ErrorEdgeNormIt=cat(2,0,ErrorEdgeNormIt);
EnergyVariationNormIt=cat(2,0,EnergyVariationNormIt);
ErrorEdgeVariationNormIt=cat(2,0,ErrorEdgeVariationNormIt);

%% First Plot
Fig = figure('Position',[100 50 1200 650]);
set(Fig,'name','Graphs','numbertitle','off','color','w')
Fig1 = subplot(2,2,1);
plot(log10(GDLTotal),real(EnergyIt),'--rs','LineWidth',2,'MarkerEdgeColor','k',...
    'MarkerFaceColor','g','MarkerSize',5);
title('Energy vs DOF','FontWeight','bold');
xlabel('(Degrees of freedom) log N');
ylabel('(Energy) U');
grid on;

%% Second Plot
p2=polyfit(log10(GDLTotal(length(ErrorEdgeNormIt)-nnz(ErrorEdgeNormIt)+1:end)),...
    log10(ErrorEdgeNormIt(length(ErrorEdgeNormIt)-nnz(ErrorEdgeNormIt)+1:end)),1);
p22=polyval(p2,log10(GDLTotal(length(ErrorEdgeNormIt)-nnz(ErrorEdgeNormIt)+1:end)));
C2={'y','=',num2str(p2(1)),'x','+',num2str(p2(2))};
str2=sprintf('%s',C2{:});
AverageY=sum(log10(ErrorEdgeNormIt(length(ErrorEdgeNormIt)-...
    nnz(ErrorEdgeNormIt)+1:end)))/...
    length(ErrorEdgeNormIt(length(ErrorEdgeNormIt)-nnz(ErrorEdgeNormIt)+1:end));
SQtotal=sum((log10(ErrorEdgeNormIt(length(ErrorEdgeNormIt)-nnz(ErrorEdgeNormIt)+1:end))-AverageY).^2);
ValReta=p2(1)*log10(GDLTotal(length(ErrorEdgeNormIt)-nnz(ErrorEdgeNormIt)+1:end))+p2(2);
SQres=sum((log10(ErrorEdgeNormIt(length(ErrorEdgeNormIt)-nnz(ErrorEdgeNormIt)+1:end))-ValReta).^2);
Rsquared=1-(SQres/SQtotal);
% str22=sprintf('R^2 = %s',num2str(Rsquared));
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

%%  Third Plot
Fig3=subplot(2,2,3);
plot(log10(GDLTotal),log10(EnergyVariationNormIt),'--rs','LineWidth',2,'MarkerEdgeColor','k',...
    'MarkerFaceColor','g','MarkerSize',5);
title('Normalized energy variation vs DOF','FontWeight','bold');
xlabel('(Degrees of freedom) log N');
ylabel('(Normalized energy variation) log \epsilon_{U}');
grid on;

%% Fourth Plot
Fig4=subplot(2,2,4);
plot(log10(GDLTotal),log10(ErrorEdgeVariationNormIt),'--rs','LineWidth',2,'MarkerEdgeColor','k',...
    'MarkerFaceColor','g','MarkerSize',5);
title('Normalized edge error density vs DOF','FontWeight','bold');
xlabel('(Degrees of freedom) log N');
ylabel('(Normalized edge error density) log \epsilon_{\Gamma} \cdot \epsilon_{\Gamma_0}^{-1}');
grid on;

end