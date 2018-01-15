function Ect = Energy_ctTri(Nodes,Loops,NGP)
% ENERGY_CTTRI iterates on the elements and calls the function that 
% generates the constant part of the solution energy.
%
% Input: The Node list, Loops structure and geometrical information loaded 
% from the HeatStructDef GUI. The number of Gauss integration points.
% Output:  Ect, the constant part of the solution energy
%
% Iterates on the elements and calls the function that generates the 
% constant part of the element energy for triangular meshes
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
% The computation of the solution energy is covered in Section 3.9 of
% reference [4]. Its calculation involves domain integration. Using the 
% notations of reference [4], 
% Ect = Integrate [Up*k*Up] = Q^2/(4k) * Integrate[R^2],
% where Q and k are the (constant) internal heat generation and
% conductivity, R is the radial coordinate.

%Initializing
Ect=0;

%% Sweeping through the elements
for ii=1:length(Loops.area)
    
   % LocLoop is a structure where the features of the current
    % element which are directly useful for the calculation of the
    % Ect are stored.    
    LocLoop = struct('id',ii,'nodes',Loops.nodes(ii,:),'center',...
        Loops.center(ii,:));
    
    % Loading the material characteristics from the Loops structure
    k = Loops.material(ii,1);
    Q = Loops.material(ii,2);

    % Getting coordinates of the nodes of the element (global)
    LocNodes = Nodes(LocLoop.nodes(:),:);
    
    % Generating the Gauss points and weights(global and local). TRIQUAD is
    % a library function created by Greg von Winckel
    [XG,YG,Wx,Wy]=triquad(NGP,LocNodes);
    
    % Computing the local coordinates of the Gauss points
    X = XG - LocLoop.center(1);
    Y = YG - LocLoop.center(2);
    
    % Computing the square of the radial coordinate 
    R2 = X.^2 + Y.^2;  

    % Computing the energy using domain integration
    Ecti = Q^2/(4*k) * Wx'*R2*Wy;
    
    % Incrementing the global Ect energy with the current loop's 
    % contribution.
    Ect= Ect + Ecti;                                          
    
end

end
