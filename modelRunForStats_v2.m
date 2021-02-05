function [pStats, xOut, R0, Dyn] = modelRunForStats_v2(nnMatrix,fracInf,VGR,VI,basalDeathRate,infDeathRate,seed)
%% Init cells, small frac infected
nCells = size(nnMatrix,1);
infected = zeros(nCells,1);%<fracInf;
infected(seed)=1; %always start with center pixel infected
x= [~infected, infected, zeros(nCells,1), zeros(nCells,1)]';
x = x(:);

% Run model

[xOut, ~, R0, Dyn] = SSA_forSpatialSIDGrids(x,nnMatrix,[],VGR,VI,basalDeathRate,infDeathRate);

pStats = sum(reshape(xOut,4,[]),2)./sum(sum(reshape(xOut,4,[]),2));
