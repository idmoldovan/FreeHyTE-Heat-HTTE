function PlotFieldsTri(Nodes, Edges, Loops,NoDiv)
% Flux postproessing - Computes the final fluxes, in x and y directions

% Plotting the Finite Element Mesh (4x)
% Initialization of the required matrices

% get the coordinates of the mesh points
xmesh = [Nodes(Edges.nini(:),1) Nodes(Edges.nfin(:),1)];
ymesh = [Nodes(Edges.nini(:),2) Nodes(Edges.nfin(:),2)];

% draw the mesh
Fig = figure;
set(Fig,'name','Refinement Orders, Temperature and Flux Fields',...
    'numbertitle','off','color','w');

%% ------------------ Refinement orders ---------------------------

delta = sqrt((max(max(xmesh))-min(min(xmesh)))^2 + ...
    (max(max(ymesh))-min(min(ymesh)))^2);

Orders = subplot(2,2,1);
hold on; title('Orders of the Loops and Edges Bases','FontWeight','bold');
axis([min(min(xmesh)) max(max(xmesh)) min(min(ymesh)) max(max(ymesh))]);
for ii=1:length(Edges.type)
    line(xmesh(ii,:),ymesh(ii,:),'color','k','LineWidth',1);
end
daspect([1 1 1]);


% Inserting the Edge's Order
for i = 1:length(Edges.type)
    EX = Nodes([Edges.nini(i) Edges.nfin(i)],1);
    EY = Nodes([Edges.nini(i) Edges.nfin(i)],2);
    pos = [sum(EX)/2+0.0*delta,sum(EY)/2+0.0*delta] ;
    text(pos(1),pos(2),int2str(Edges.order(i)),'fontsize',8, ...
        'fontweight','bold','color','g');
end
% Inserting the Loop's Order
for i = 1:length(Loops.area)
    pos = polygonCentroid(Nodes(Loops.nodes(i,:),:));
    text(pos(1),pos(2),int2str(Loops.order(i,1)),'fontsize',8, 'fontweight','bold','color','b');
end
%
axis off;

% ------------------ Temperature ---------------------------
Temp = subplot(2,2,2);
hold on; title('Temperature Field','FontWeight','bold');
axis([min(min(xmesh)) max(max(xmesh)) min(min(ymesh)) max(max(ymesh))]);
for ii=1:length(Edges.type)
    line(xmesh(ii,:),ymesh(ii,:),'color','b');
end
daspect([1 1 1]);

% ------------------ Flux x direction ---------------------------
Qx = subplot(2,2,3);
hold on; title('Flux field, x direction','FontWeight','bold');
axis([min(min(xmesh)) max(max(xmesh)) min(min(ymesh)) max(max(ymesh))]);
for ii=1:length(Edges.type)
    line(xmesh(ii,:),ymesh(ii,:),'color','b');
end
daspect([1 1 1]);
% ------------------ Flux y direction ---------------------------
Qy = subplot(2,2,4);
hold on; title('Flux field, y direction','FontWeight','bold');
axis([min(min(xmesh)) max(max(xmesh)) min(min(ymesh)) max(max(ymesh))]);
for ii=1:length(Edges.type)
    line(xmesh(ii,:),ymesh(ii,:),'color','b');
end
daspect([1 1 1]);

%%
% Start plotting

for ii=1:length(Loops.area)
    
    LocLoop = struct('id',ii,'nodes',Loops.nodes(ii,:),'center',...
        Loops.center(ii,:),'T',Loops.T{ii}, 'Qx',Loops.Qx{ii},...
        'Qy',Loops.Qy{ii});
    
    % Getting coordinates of the nodes of the element (global)
    LocNodes = Nodes(LocLoop.nodes(:),:);
    
    [GlobalX,GlobalY,~,~]=triquad(NoDiv,LocNodes);
    
    subplot(Temp);
    contourf(GlobalX, GlobalY, real(LocLoop.T), 50, 'edgecolor','none');
    colormap jet;
    
    subplot(Qx);
    contourf(GlobalX,GlobalY,real(LocLoop.Qx), 50, 'edgecolor','none');
    
    subplot(Qy);
    contourf(GlobalX,GlobalY,real(LocLoop.Qy), 50, 'edgecolor','none');
    
end

subplot(Temp); colorbar; axis off;    % put on the legends
subplot(Qx); colorbar; axis off;
subplot(Qy); colorbar; axis off;

end