function Ect = Energy_ctReg(Loops)
% sweeps through the elements and calls the function that generates the 
% constant part of the elenergy,
% Ect = Integrate [Up*Up]
% with analytic expression for rectangular domains.
%%
load('HeatStructDef','L','B','Nx','Ny');
%Initializing
Ect=0;
l = L/Nx;
b = B/Ny;

for ii=1:length(Loops.area)
    
    k = Loops.material(ii,1);
    Q = Loops.material(ii,2);
    
    Ecti = Q^2/(96*k) * l * b * (l^2 + b^2);
    Ect= Ect + Ecti;                                          
    
end
end