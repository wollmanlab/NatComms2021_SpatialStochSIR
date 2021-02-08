function [xOut, loadWhenInf, Rt] = SSA_forSpatialSIDGrids_SimpleAnnotated(x,neighMat,...
    VGR,VI,basalDeathRate,infDeathRate)
% inputs are:
%   x - nSpecies * nCells boolean vector i.e. 
%       [1 0 0 0 , 0 1 0 0, 1 0 0 0,...]  corresponds to [alive, infected, alive, ...]
%   neighMat - nReactions X nReactions matrix specifying which cells are nearest neighbors
%   VGR - Viral growth rate in units of infection EC50(S4C-D)
%   VI - Viral infectivity (S4C-D)
%   basalDeathRate - net uninfected death rate for a specific TNF concentration
%   infDeathRate - net infected death rate for a specific TNF concentration


% Define propensity functions in terms of parameters, k, and states, x.
% each cell can be in 4 states - SIDF
w_k_x = @(x,k) [k(1)*x(1), k(2)*x(1), k(3)*x(2)];
% Specify stoichiometry for a single cell
S = [ -1 -1  0 ;...
    1  0 -1 ;...
    0  1  0 ;
    0  0  1];
%reactions: infection, false positive death, infected death

global nReactions nSpecies
nReactions = 3;
nSpecies = 4;

%number of cells
nC = numel(x)/nSpecies;

%make full stoichiometry matrix, block diagonal of single cell stoch mat
SCell = repmat({S}, 1, nC);
S = sparse(blkdiag(SCell{:}));

tstop = 48;%finaltime = 48 hours


[xOut, loadWhenInf, Rt] = Run_SSA(w_k_x,S,x,tstop,nC, neighMat,VGR,VI,basalDeathRate, infDeathRate); % call code to run stochastic simulation.

function [x, loadWhenInf, Rt] = Run_SSA(prop_fun,S,x0,tstop,nC,neighMat,VGR,VI,basalDeathRate, infDeathRate)

global nReactions nSpecies

t=0;
x = x0;     %% Specify initial conditions


%Init params
k=zeros(1,nC*nReactions);
k(2:nReactions:end)=basalDeathRate;%basal death rate
k(3:nReactions:end) = infDeathRate; %virus induced death

loadWhenInf = cell(size(neighMat,1),1);
%vector of virus infection. 0 for healthy or dead cells. initial infection=1 for infected cells.
v=x(2:nSpecies:end);


while t<tstop
    %update k
    NNvLoad = neighMat*(v(:).*x(2:nSpecies:end)); %nearest neighbor viral load in units of EC50
    k(1:nReactions:end) = VI./(1+1./NNvLoad);  %infection probability follows hill function (S4.D)
    
    %Propensity functions: all the probabilities of all the reactions
    w = reshape(cell2mat(arrayfun(@(y,z) prop_fun(y{:},z{:}), mat2cell(x(:),nSpecies*ones(1,nC)), mat2cell(k(:),nReactions*ones(1,nC)),'uniformoutput', false))',1,[]);%sorry for the oneliner
    w0 = sum(w);                               % sum of the prop. functions
    dt = 1/w0*log(1/rand);          %when's the next reaction?
    t = t+dt;                       % update time of next reaction, exp dist
    if t<=tstop
        r2w0=rand*w0;               % generate second random number and multiply by prop. sum
        i=1;                                          % initialize reaction counter
        while sum(w(1:i))<r2w0             % what's the next reaction? increment counter until sum(w(1:i)) exceeds r2w0
            i=i+1;
        end
        x = x+S(:,i);                                 % update the configuration (apply reaction)
        
        if mod(i,3)==1 %for Rt, if a cell gets infected, keep it's nn viral load
            loadWhenInf{ceil(i/3)} = neighMat(ceil(i/3),:).*v';
        end
        v = v+VGR*dt*x(2:nSpecies:end);                      % increase viral load for infected live cells
    end
    
end

NNLoad = loadWhenInf;
C1 = cellfun(@(x) x./sum(x), NNLoad, 'UniformOutput', false);
indsOfCellsThatGotInfected = find(~cellfun(@isempty, NNLoad));

cellsIInfected = zeros(size(NNLoad));

for i=indsOfCellsThatGotInfected'
    [~,ind,prob] = find(C1{i});
    cellsIInfected(ind) = cellsIInfected(ind) + prob';
end

a = k(1:nReactions:end)./ infDeathRate; %probability of getting infected devided by probability of infecting cell to die gives Rt estimate
a(a==0)=[];
Rt = mean(a);
if isnan(Rt)
    R0=0;
end

