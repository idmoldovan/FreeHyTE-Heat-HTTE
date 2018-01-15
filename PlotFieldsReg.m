function PlotFieldsReg(Nodes, Edges, Loops,abscissa)
% PLOTFIELDSREG plots the final temperature and flux fields.
%
% PlotFieldsReg is called by MAINREG
%
% Input:
%  structures Edges and Loops where the final temperature and flux fields
%  stored, the node position list Node and the abscissa of the
%  Gauss-Legendre integration points. 
% Output:
%  4 figures - three color maps and one plot with final refinement orders 
%
% This is a plotting-only function. The values of the final temperature and
% heat flux values were computed and stored in the Loops structure in
% COMPUTEFIELDSREG.
%

% Getting the (global, Cartesian) coordinates of the mesh points
xmesh = [Nodes(Edges.nini(:),1) Nodes(Edges.nfin(:),1)];
ymesh = [Nodes(Edges.nini(:),2) Nodes(Edges.nfin(:),2)];

%% Drawing the mesh, to prepare the plotting of the fields
% The figure Fig contains 4 plots, namely one map containing the refinement
% orders, one colour map of the temperature field and two colour maps of
% the flux fields in x and y directions.
Fig = figure;
set(Fig,'name','Refinement Orders, Temperature and Flux Fields',...
    'numbertitle','off','color','w');

%% Plotting the refinement order map
% getting the global dimension of the domain, for scaling purposes
delta = sqrt(max(max(xmesh))^2+max(max(ymesh))^2);

% Plotting the approximation order map
Orders = subplot(2,2,1);
hold on; title('Orders of the Loops and Edges Bases','FontWeight','bold');
axis([-0.05 1.05*max(max(xmesh)) -0.05 1.05*max(max(ymesh))]);
% Drawing the mesh, one edge at a time
for ii=1:length(Edges.type)
    line(xmesh(ii,:),ymesh(ii,:),'color','k','LineWidth',1);
end
daspect([1 1 1]);

% Writing the edges Order
for i = 1:length(Edges.type)
    % getting the limits of the edge
    EX = Nodes([Edges.nini(i) Edges.nfin(i)],1);
    EY = Nodes([Edges.nini(i) Edges.nfin(i)],2);
    % computing the insertion point for the text
    pos = [sum(EX)/2+0.01*delta,sum(EY)/2+0.01*delta] ;
    % writing the text
    text(pos(1),pos(2),int2str(Edges.order(i)),'fontsize',8, ...
        'fontweight','bold','color','g');
end
% Writing the Loops Order
for i = 1:length(Loops.area)
    % getting the limits of the element
    EX = Nodes(Loops.nodes(i,:),1);EY = Nodes(Loops.nodes(i,:),2) ;
    % computing the insertion point for the text
    pos = [sum(EX)/4,sum(EY)/4] ;
    % writing the text
    text(pos(1),pos(2),int2str(Loops.order(i,1)),'fontsize',8,...
        'fontweight','bold','color','b');
end
%
axis off;


%% Temperature plot
% Preparing the mesh for the temperature plot
Temp = subplot(2,2,2);
hold on; title('Temperature Field','FontWeight','bold');
axis([min(min(xmesh)) max(max(xmesh)) min(min(ymesh)) max(max(ymesh))]);
% Drawing the mesh, one edge at a time
for ii=1:length(Edges.type)
    line(xmesh(ii,:),ymesh(ii,:),'color','b');
end
daspect([1 1 1]);

%% Heat flux plot, in X
% Preparing the mesh for the Qx plot
Qx = subplot(2,2,3);
hold on; title('Flux field, x direction','FontWeight','bold');
axis([min(min(xmesh)) max(max(xmesh)) min(min(ymesh)) max(max(ymesh))]);
% Drawing the mesh, one edge at a time
for ii=1:length(Edges.type)
    line(xmesh(ii,:),ymesh(ii,:),'color','b');
end
daspect([1 1 1]);

%% Heat flux plot, in Y
% Preparing the mesh for the Qy plot
Qy = subplot(2,2,4);
hold on; title('Flux field, y direction','FontWeight','bold');
axis([min(min(xmesh)) max(max(xmesh)) min(min(ymesh)) max(max(ymesh))]);
% Drawing the mesh, one edge at a time
for ii=1:length(Edges.type)
    line(xmesh(ii,:),ymesh(ii,:),'color','b');
end
daspect([1 1 1]);

%% Initialization of the output points
% The plotting points are, by default, the Gauss points
pts = (abscissa)';


%% Plotting of the temperature and heat flux fields
% Sweeping the elements to draw the colour maps
for ii=1:length(Loops.area)
    
    % LocLoop is a local structure where the features of the current
    % element which are directly useful for the plotting of the output
    % fields are stored.
    LocLoop = struct('id',ii,'nodes',Loops.nodes(ii,:),'center',...
        Loops.center(ii,:),'T',Loops.T{ii}, 'Qx',Loops.Qx{ii},...
        'Qy',Loops.Qy{ii});
    
    %% Generating the geometric data
    % Computing the length of the sides of the element in x and y
    % direction. sze simply collects the distances between the two
    % most far apart points of the element in x and y direcntions.
    sze = max(Nodes(LocLoop.nodes(:),:))-min(Nodes(LocLoop.nodes(:),:));
    
    % Generating the output grid in local coordinates.
    x = 1/2*sze(1)*pts;
    y = 1/2*sze(2)*pts;
    [xd,yd] = ndgrid(x,y);
    
    % Start plotting fields on element
    GlobalX = xd + LocLoop.center(1);
    GlobalY = yd + LocLoop.center(2);
    
    %% Plotting the temperature and flux fields as colour maps
    
    subplot(Temp);
    contourf(GlobalX, GlobalY, real(LocLoop.T), 50, 'edgecolor','none');
    colormap jet;
    
    subplot(Qx);
    contourf(GlobalX,GlobalY,real(LocLoop.Qx), 50, 'edgecolor','none');
    
    subplot(Qy);
    contourf(GlobalX,GlobalY,real(LocLoop.Qy), 50, 'edgecolor','none');
    
end

% Plotting the legends
subplot(Temp); caxis(caxis); colorbar; axis off;    
subplot(Qx); caxis(caxis); colorbar; axis off;
subplot(Qy); caxis(caxis); colorbar; axis off;

end
