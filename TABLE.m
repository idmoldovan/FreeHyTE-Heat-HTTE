function TABLE(GDL,Beta,...
    EnergyIt,LoopsOrderIt,EdgesOrderIt,...
    RefinedEdgesIt,RefinedLoopsIt) 


Data=cell(length(GDL),8);
for i= 1:length(GDL)
    Data{i,1}=i;
    Data{i,2}=num2str(GDL(i));
    Data{i,3}=num2str(Beta(i,:));
    Data{i,4}=num2str(RefinedLoopsIt{i,:});
    Data{i,5}=num2str(LoopsOrderIt(:,i)');
    Data{i,6}=num2str(RefinedEdgesIt{i,:});
    Data{i,7}=num2str(EdgesOrderIt(:,i)');
    Data{i,8}=num2str(EnergyIt(i));
end
cnames={'Iteration','DOF','Beta','Refined Loops','Loops Order',...
    'Refined Edge','Edges Order','Energy'};
maxLen = zeros(1,8);
for i=1:8
    for j=1:length(GDL)
        len = length(Data{j,i});
        if (len>maxLen(i))
            maxLen(i) = len;
        end
    end
end
cellMaxLen = num2cell(ceil(maxLen*4.2));
for i=1:8
    if cellMaxLen{i}<length(cnames{i})
        cellMaxLen{i}=length(cnames{i})+10;
    end
end

f = figure('Position',[50 100 ...
    70+50+cellMaxLen{3}+20+100+cellMaxLen{5}+20+100+cellMaxLen{7}+70+100 450],...
    'color','w');
set(f,'name','Table of iterations info','numbertitle','off');
table = uitable('Parent',f,'units','pixels');
set(table,'Data',Data,'ColumnName',cnames,...
    'RowName',[],'Position',[10 25 70+50+cellMaxLen{3}+20+100+cellMaxLen{5}+20+100+cellMaxLen{7}+70+60 400],...
    'ColumnWidth',{70 50 cellMaxLen{3}+20 100 cellMaxLen{5}+20 100 cellMaxLen{7} 70}); 
%set(f,'menubar','none');
end