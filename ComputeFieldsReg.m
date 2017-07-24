function [Loops] = ComputeFieldsReg(Nodes, Loops, X, iteration, abscissa)
% Postprocessing - Computes the final temperatures and flluxes, writes them
% in a TecPlot format file and delivers the data for posterior plotting in Matlab

%% Pre-flight

% Creating the sub-folder to store the analysis results
load('HeatStructDef','DirName','FileName','RunOption',...
    'EdgesOrder','LoopsOrder'); %loads the folder and the filename defined in the GUI

if  ~isempty(DirName) % if the user left DirName blank, does not generate written output
    if RunOption == 1 % Single run
        % Generating the file name
        FileNameExt = sprintf('_ND%dNB%d.dat',LoopsOrder,EdgesOrder);
        FileName = strcat(FileName(1:end-4),FileNameExt); % end-4 removes the mat extension
        % link the path to the filename
        UFilename = fullfile(DirName,FileName);
        FileU = fopen(UFilename,'w');
        fprintf(FileU,'TITLE="%s"\n',FileName);
        fprintf(FileU,'VARIABLES="X", "Y", "T", "Qx", "Qy"\n');
        fprintf(FileU,'ZONE T="ND = %d; NB = %d"\n',LoopsOrder,EdgesOrder);
    else % Adpative run
        % Generating the file name
        FileNameExt = sprintf('_It%d.dat',iteration);
        FileName = strcat(FileName(1:end-4),FileNameExt); % end-4 removes the mat extension
        % link the path to the filename
        UFilename = fullfile(DirName,FileName);
        FileU = fopen(UFilename,'w');
        fprintf(FileU,'TITLE="%s"\n',FileName);
        fprintf(FileU,'VARIABLES="X", "Y", "T", "Qx", "Qy"\n');
        fprintf(FileU,'ZONE T="Iteration = %d"\n',iteration);
    end
end

%% Generating the mesh of output points. A very small shift is induced
% to the limits of the elements in order to avoid the duplication of
% boundary points on adjacent elements, which may comprimise the
% visibility of the continuity gaps. The shift is not symmetric to avoid
% potential problems with the flux field when x and y are absolute zero.
%pts = linspace(-(1-1.e-3),(1-1.1e-3),NoDiv+1);
pts = abscissa';


%% ************ BEGIN SWEEPING THE ELEMENTS *************
for ii = 1:length(Loops.area)
    %
    LocLoop = struct('id',ii,'nodes',Loops.nodes(ii,:),'center',...
        Loops.center(ii,:),'order',Loops.order(ii),...
        'insert',Loops.insert(ii),'dim',Loops.dim(ii));
    
    k = Loops.material(ii,1);
    Q = Loops.material(ii,2);
    
    % Computing the length of the sides of the element in x and y
    % direction. sze simply collects the distances between the two
    % most far apart points of the element in x and y direcntions.
    % ASSUMES THAT THE ELEMENT IS RECTANGULAR!
    sze = max(Nodes(LocLoop.nodes(:),:))-min(Nodes(LocLoop.nodes(:),:));
          
    % Generating the output grid (local coordinates)
    m=-LocLoop.order:LocLoop.order;
    x=1/2*sze(1)*pts;  
    y=1/2*sze(2)*pts;
    [Xd,Yd,M] = ndgrid(x,y,m);
    
   %polar coordinates, local
   
    R = sqrt(Xd.^2 + Yd.^2);  
    Th = atan2(Yd, Xd);
        
    % loads the solution vector for the current loop
    Xi(1,1,:) = X(LocLoop.insert:LocLoop.insert+LocLoop.dim-1);
    
    % Temperature shape functions and the temperature in the grid points
    Ustar = R.^abs(M) .* exp(1i*Th.*M);
    Up = -(Q/(4*k)).* R(:,:,1).^2;
    T = sum(bsxfun(@times,Ustar,Xi),3)+Up;
    
    % Derivative of the temperature shape functions in r and the Sr
    Sr = (-1)*k* abs(M) .* R.^(abs(M)-1) .* exp(1i*Th.*M); 
    SPr = (Q/2)*R(:,:,1);
    
    % Derivative of the temperature shape functions in th and the Sth
    Sth = (-1)*k*1i*M.*R.^(abs(M)-1) .* exp(1i*Th.*M);
    
    qr= sum(bsxfun(@times,Sr,Xi),3)+SPr;
    qth= sum(bsxfun(@times,Sth,Xi),3)+0;
    
    % Creating qx and qy
    Qx= cos(Th(:,:,1)).*qr - sin(Th(:,:,1)).*qth;
    Qy= sin(Th(:,:,1)).*qr + cos(Th(:,:,1)).*qth;
    
    
    clear Xi m;
   
    
    % Writing the TecPlot data file
    if ~isempty(DirName) % if the user left DirName blank, does not generate written output    
        
        % Computing the global coordinates
        [xd,yd] = ndgrid(x,y);
        GlobalX = xd + LocLoop.center(1);
        GlobalY = yd + LocLoop.center(2);
        
        % Writting the fields in a TecPlot compatible file
        for kk = 1:length(abscissa)           
            for ll = 1:length(abscissa)        
                fprintf(FileU,'%0.6e %0.6e %0.6e %0.6e %0.6e \n',...
                    GlobalX(kk,ll), GlobalY(kk,ll), real(T(kk,ll)), ...
                    real(Qx(kk,ll)), real(Qy(kk,ll)));
            end
        end
    end
    
    Loops.T(ii) = {T};
    Loops.Qx(ii) = {Qx};
    Loops.Qy(ii) = {Qy};
    
end

%%
if ~isempty(DirName)
    fclose(FileU);
end

end