function PlotFieldsReg(Nodes, Edges, Loops,abscissa)
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

delta = sqrt(max(max(xmesh))^2+max(max(ymesh))^2);

Orders = subplot(2,2,1);
hold on; title('Orders of the Loops and Edges Bases','FontWeight','bold');
axis([-0.05 1.05*max(max(xmesh)) -0.05 1.05*max(max(ymesh))]);
for ii=1:length(Edges.type)
    line(xmesh(ii,:),ymesh(ii,:),'color','k','LineWidth',1);
end
daspect([1 1 1]);


% Inserting the Edge's Order
    for i = 1:length(Edges.type)
        EX = Nodes([Edges.nini(i) Edges.nfin(i)],1);
        EY = Nodes([Edges.nini(i) Edges.nfin(i)],2);
        pos = [sum(EX)/2+0.01*delta,sum(EY)/2+0.01*delta] ;
        text(pos(1),pos(2),int2str(Edges.order(i)),'fontsize',8, ...
            'fontweight','bold','color','g');
    end
% Inserting the Loop's Order
    for i = 1:length(Loops.area)
        EX = Nodes(Loops.nodes(i,:),1) ;EY = Nodes(Loops.nodes(i,:),2) ;
        pos = [sum(EX)/4,sum(EY)/4] ;
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

%% Generating the mesh of output points. A very small shift is induced
% the limits of the elements in order to avoid the duplication of
% boundary points on adjacent elements, which may compromise the
% visibility of the continuity gaps. The shift is not symmetric to avoid
% potential problems with the flux field when x and y are absolute zero.

%pts = linspace(-(1-1.e-3),(1-1.1e-3),NoDiv+1);
pts = (abscissa)';


%%
for ii=1:length(Loops.area)
    
    LocLoop = struct('id',ii,'nodes',Loops.nodes(ii,:),'center',...
        Loops.center(ii,:),'T',Loops.T{ii}, 'Qx',Loops.Qx{ii},...
        'Qy',Loops.Qy{ii});
    
    % Computing the length of the sides of the element in x and y
    % direction. sze simply collects the distances between the two
    % most far apart points of the element in x and y direcntions.
    % ASSUMES THAT THE ELEMENT IS RECTANGULAR!
    sze = max(Nodes(LocLoop.nodes(:),:))-min(Nodes(LocLoop.nodes(:),:));
    
    x = 1/2*sze(1)*pts;
    y = 1/2*sze(2)*pts;
    [xd,yd] = ndgrid(x,y);
    
    
    % Start plotting fields on element
    
    GlobalX = xd + LocLoop.center(1);
    GlobalY = yd + LocLoop.center(2);
    
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