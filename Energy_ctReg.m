function Ect = Energy_ctReg(Loops)
% ENERGY_CTREG iterates on the elements and calls the function that 
% generates the constant part of the solution energy, with analytic 
% expression for rectangular domains.
%
% Input: The Loops structure and geometrical information loaded from 
% the HeatStructDef GUI
% Output:  Ect, the constant part of the solution energy
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
% reference [4]. Its calculation is analytic. Using the notations of 
% reference [4], 
% Ect = Q^2/(96k)*l*b*(l^2 + b^2)
% where Q and k are the (constant) internal heat generation and
% conductivity, and l and b are the dimensions of the element.

%% Initialization
% Loading the global dimensions of the structure and the number of elements
% in the Cartesian directions from the HeatStructDef GUI.
load('HeatStructDef','L','B','Nx','Ny');
Ect=0;
l = L/Nx; 
b = B/Ny;

%% Sweeping through the elements
for ii=1:length(Loops.area)
    
    k = Loops.material(ii,1);
    Q = Loops.material(ii,2);
    
    % Computing the energy using the analytic expression
    Ecti = Q^2/(96*k) * l * b * (l^2 + b^2);
    
    % Incrementing the global Ect energy with the current loop's 
    % contribution.
    Ect= Ect + Ecti;                                          
    
end
end
