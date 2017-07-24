function [ErrorEnergyNorm,Stop]=CheckMaxDegrees(Loops,Edges,ErrorEnergyNorm)
% This function checks if the maximum order of refinement was reached for
% any element or essential boundary. It stops the code if it was.

load('Adaptive','MaxOrder');
Stop = 0;

for j=1:length(Loops.area)
    if (Loops.order(j) >= MaxOrder) 
          ErrorEnergyNorm=0.0;
          Stop=1;
    end
end

for i=1:length(Edges.type)
      if (Edges.order(i) >= MaxOrder)
          ErrorEnergyNorm=0.0;
          Stop=1;
      end
end

end

