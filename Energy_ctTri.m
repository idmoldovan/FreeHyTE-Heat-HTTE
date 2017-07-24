function Ect = Energy_ctTri(Nodes,Loops,NGP)
% sweeps through the elements and calls the function that generates the 
% constant part of the elenergy,
% Ect = Integrate [Up*k*Up] = Q^2/(4k) * Integrate[r^2].


%Initializing
Ect=0;

for ii=1:length(Loops.area)
    
    LocLoop = struct('id',ii,'nodes',Loops.nodes(ii,:),'center',...
        Loops.center(ii,:));
    
    k = Loops.material(ii,1);
    Q = Loops.material(ii,2);

    % Getting coordinates of the nodes of the element (global)
    LocNodes = Nodes(LocLoop.nodes(:),:);
    
    % Generating the Gauss points and weights(global and local)
    [XG,YG,Wx,Wy]=triquad(NGP,LocNodes);
    X = XG - LocLoop.center(1);
    Y = YG - LocLoop.center(2);
    
    R2 = X.^2 + Y.^2;  % polar coordinates, local

    Ecti = Q^2/(4*k) * Wx'*R2*Wy;
    
    Ect= Ect + Ecti;                                          
    
end

end