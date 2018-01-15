function [Loops] = ComputeFieldsTri(Nodes, Loops, X, iteration, NoDiv)
% COMPUTEFIELDSTRI computes the final temperature and flux fields, and
% stores them in the Loops structure and in an output file, if requested
% by the user.
%
% COMPUTEFIELDSTRI is called by MAINTRI.
% Input:
% It receives structures Edges and  Loops, the solution vector X
% (calculated in MAINTRI), the node list Nodes and the iteration number.
% Besides the input arguments, COMPUTEFIELDSTRI loads from HeatStructDef.mat
% information regarding the file where the results should be stored.
% Ouput:
% The Loops structures updated with the temperatures and flux values in
% the (NGP+1)^2 plotting points.
% NGP is the number of Gauss-Legendre integration points
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
% The estimates of the temperature and flux fields are obtained by
% substituting the solution X of the solving system in the domain
% approximations of the respective fields. This is done element by element.
% The estimates are used for:
% * storing the temperature and flux estimates, calculated in the
% (NGP+1)^2 plotting points, in Loop structure returned to MAINTRI,
%  for some posterior plotting;
% * storing the temperature and flux estimates, calculated in the
% (NGP+1)^2 plotting points in a result file which can be accessed by
% third party post-processing software (it is pre-formatted for TecPlot).
% To avoid overloading the disk with (fairly large) undesired result files,
% the result file is only created if the input was saved in the first GUI
% (see Section 5.2.2 of reference [3]).
%
% The computation of the field estimates is further covered in Section 7.1
% of reference [2] and Section 4.5 of reference [4]. The output options are
% explained in detail in Section 7.2 of reference [2].

%% Managing the result file
% Loading the save file data from HeatStructDef.
% * DirName is the folder to save the result file;
% * FileName is the name of the file where the analysis data was saved;
% * RunOption in case of adaptive algortihm the file name is created for
% each iteration
% * EdgesOrder and LoopsOrder are the orders of the approximation bases of
% the edges and elements, as defined in HeatStructDef GUI.
load('HeatStructDef','DirName','FileName','RunOption','EdgesOrder','LoopsOrder');

% if the user left DirName blank, does not generate written output
if  ~isempty(DirName)
    if RunOption == 1 % Single run
        % Generating the file name
        % adding information on the refinement orders to the FileName
        FileNameExt = sprintf('_ND%dNB%d.dat',LoopsOrder,EdgesOrder);
        % removing .mat extension (end-4)
        FileName = strcat(FileName(1:end-4),FileNameExt);
        % link the path (DirName) to the file name (FileName). This generates
        % the UFilename string, which allows the creation of the result file
        UFilename = fullfile(DirName,FileName);
        % open result file for writing
        FileU = fopen(UFilename,'w');
        
        % The following lines are the header needed to prepare the result file
        % to be readable with TecPlot. If you wish to use another software for
        % post-processing, you probably need to change the header.
        % However, the format of the data should be fairly general. The results
        % are written as a matrix of NEL*(NGP+1)^2 lines (where NEL is the
        % total number of finite elements and (NGP+1)^2 the total number of
        % result points per element), and 5 columns. For each result point, the
        % columns list the x and y coordinates, in the global referential,
        % followed by the values of the temperature and flux fields.
        fprintf(FileU,'TITLE="%s"\n',FileName);
        fprintf(FileU,'VARIABLES="X", "Y", "T", "Qx", "Qy"\n');
        fprintf(FileU,'ZONE T="ND = %d; NB = %d"\n',LoopsOrder,EdgesOrder);
        
        
    else % Adaptive run
        % Generating the file name adding information on the current iteration
        FileNameExt = sprintf('_It%d.dat',iteration);
        % removing .mat extension (end-4)
        FileName = strcat(FileName(1:end-4),FileNameExt);
        % link the path to the filename
        UFilename = fullfile(DirName,FileName);
        % The following lines are the header needed to prepare the result file
        % to be readable with TecPlot. If you wish to use another software for
        % post-processing, you probably need to change the header.
        % However, the format of the data should be fairly general. The results
        % are written as a matrix of NEL*(NGP+1)^2 lines (where NEL is the
        % total number of finite elements and (NGP+1)^2 the total number of
        % result points per element), and 5 columns. For each result point, the
        % columns list the x and y coordinates, in the global referential,
        % followed by the values of the temperature and flux fields.
        FileU = fopen(UFilename,'w');
        fprintf(FileU,'TITLE="%s"\n',FileName);
        fprintf(FileU,'VARIABLES="X", "Y", "T", "Qx", "Qy"\n');
        fprintf(FileU,'ZONE T="Iteration = %d"\n',iteration);
    end
end

%% ************ BEGIN SWEEPING THE ELEMENTS *************
% Sweeping the elements to compute the solution fields
for ii = 1:length(Loops.area)
    
    %% Initialization
    % LocLoop is a local structure where the features of the current
    % element which are directly useful for the calculation of the output
    % fields are stored.
    LocLoop = struct('id',ii,'nodes',Loops.nodes(ii,:),'center',...
        Loops.center(ii,:),'order',Loops.order(ii),...
        'insert',Loops.insert(ii),'dim',Loops.dim(ii));
    
    k = Loops.material(ii,1);
    Q = Loops.material(ii,2);
    
    %% Generating the geometric data
    % Getting coordinates of the nodes of the element (global)
    LocNodes = Nodes(LocLoop.nodes(:),:);
    
    % Generating the output points in the global Cartesian referential.
    % The output points belong to the Gauss-Legendre quadrature of the
    % triangular element.
    [GlobalX,GlobalY,~,~]=triquad(NoDiv,LocNodes);
    % Generating the output grid in local coordinates.
    x = GlobalX - LocLoop.center(1);
    y = GlobalY - LocLoop.center(2);
    
    % Vector containing the orders of the basis
    m(1,1,:) = -LocLoop.order:LocLoop.order;
    
    % Transforming the local Cartesian coordinates into polar.
    r = sqrt(x.^2 + y.^2);
    t = atan2(y,x);
    % Generating the 3D matrices, for consistency with the programming
    % strategy used in the regular meshes.
    R = repmat(r,[1 1 size(m,3)]);
    Th = repmat(t,[1 1 size(m,3)]);
    M = repmat(m,[NoDiv NoDiv 1]);
    
    % loads the solution vector for the current loop
    Xi(1,1,:) = X(LocLoop.insert:LocLoop.insert+LocLoop.dim-1);
    
    %% Computing the basis functions
    % Computing temperature shape functions and the temperature values
    % in the grid points in the polar referential.
    % For a full description of the basis, please refer to Section 3 of
    % reference [4].
    Ustar = R.^abs(M) .* exp(1i*Th.*M);
    Up = -(Q/(4*k)).* R(:,:,1).^2;
    T = sum(bsxfun(@times,Ustar,Xi),3)+Up;
    
    % Computing the flux shape functions, in the R direction...
    Sr = (-1)*k* abs(M) .* R.^(abs(M)-1) .* exp(1i*Th.*M);
    SPr = (Q/2)*R(:,:,1);
    
    % ... and in THETA (no internally generated flux in THETA)
    Sth = (-1)*k*1i*M.*R.^(abs(M)-1) .* exp(1i*Th.*M);
    
    % Computing the flux values qr and qth in the grid points, in the polar
    % referential.
    qr= sum(bsxfun(@times,Sr,Xi),3)+SPr;
    qth= sum(bsxfun(@times,Sth,Xi),3)+0;
    
    % Transforming the flux field qr and qth from the polar to the Cartesian
    % referential. All pages in Th are equal, so the first one is selected 
    % to compute the  normal cosines.
    Qx= cos(Th(:,:,1)).*qr - sin(Th(:,:,1)).*qth;
    Qy= sin(Th(:,:,1)).*qr + cos(Th(:,:,1)).*qth;
    
    % Clearing the Xi variables for reuse in the next element
    clear Xi m;
    
    % Writing the fields in the TecPlot compatible file, if requested by
    % the user. The results are written as a matrix of NEL*(NGP+1)^2 lines
    % (where NEL is the total number of finite elements and (NGP+1)^2 the
    % total number of result points per element), and 5 columns. For each
    % result point, the columns list the x and y coordinates, in the global
    % referential, followed by the values of the temperature and flux
    % fields.
    
    % if the user left DirName blank, does not generate written output
    if ~isempty(DirName)
        % Writting the fields in a TecPlot compatible file
        for kk = 1:NoDiv
            for ll = 1:NoDiv
                fprintf(FileU,'%0.6e %0.6e %0.6e %0.6e %0.6e \n',...
                    GlobalX(kk,ll), GlobalY(kk,ll), real(T(kk,ll)), ...
                    real(Qx(kk,ll)), real(Qy(kk,ll)));
            end
        end
    end

    % Storing the temperature and heat flux values in the Loops data
    % structure.    
    Loops.T(ii) = {T};
    Loops.Qx(ii) = {Qx};
    Loops.Qy(ii) = {Qy};
    
end

%%
if ~isempty(DirName)
    fclose(FileU);
end

end
