function [xOut, loadWhenInf, R0, InfectedDyn] = SSA_forSpatialSIDGrids(x,neighMat,Grid,VGR,VI,basalDeathRate,infDeathRate, varargin)
% Define propensity functions in terms of parameters, k, and states, x.
% each cell can be in 4 states - SIDF
w_k_x = @(x,k) [k(1)*x(1), k(2)*x(1), k(3)*x(2)];
% Specify stoichiometry
S = [ -1 -1  0 ;...
       1  0 -1 ;...
       0  1  0 ;
       0  0  1];

%reactions: infection, false positive death, infected death

%number of cells
global nReactions nSpecies
nReactions = 3;
nSpecies = 4;

nC = numel(x)/nSpecies;

%make full stoichiometry matrix, block diagonal of single cell stoch mat
SCell = repmat({S}, 1, nC);
S = sparse(blkdiag(SCell{:}));

tstop = ParseInputs('tstop', 48, varargin);%finaltime


[xOut, loadWhenInf, R0, InfectedDyn] = Run_SSA(w_k_x,S,x,tstop,nC, neighMat,Grid,VGR,VI,basalDeathRate, infDeathRate,varargin{:}); % call code to run stochastic simulation.

function [x, loadWhenInf, R0, Dyn] = Run_SSA(prop_fun,S,x0,tstop,nC,neighMat,Grid,VGR,VI,basalDeathRate, infDeathRate, varargin)


global nReactions nSpecies
Dyn = struct('t',[],'infected',[],'infDead',[]);

t=0;
x = x0;     %% Specify initial conditions

arg.plot=ParseInputs('plot', false,varargin);
arg.save=ParseInputs('save', false,varargin);
arg.savpath=ParseInputs('savpath', '',varargin);
if ~isempty(arg.savpath)
    arg.save = true;
end
arg.TNF=ParseInputs('TNF', '',varargin);

if arg.save && isempty(arg.savpath)
    error('If you want to save you must specify a savpath argument')
end

%Init params
k=zeros(1,nC*nReactions);
k(2:nReactions:end)=basalDeathRate;%basal death rate
k(3:nReactions:end) = infDeathRate; %virus induced death 


loadWhenInf = cell(size(neighMat,1),1);
%vector of virus infection. 0 for healthy or dead cells. initial infection=1 for infected cells.
v=x(2:nSpecies:end);

if arg.plot
    figure('color','w','Position',[100,100, 450, 450])
    ax1 = axes('position',[0.02,0.1,0.5,0.5]);
    set(gcf, 'renderer', 'painters');
    pointSize = 40;
    lineSize = 2;
end

if arg.save
    outputVideo = VideoWriter(sprintf('%s%s',arg.savpath,'.avi'));
    outputVideo.FrameRate = 16;
    open(outputVideo);
end
frameNum=-1;
while t<tstop
    
    if arg.plot
        
        liveColor = [64,224,208]/255; %cyan
        virusColor = [148, 0, 211]/255; %purple
        deathColor = [255,165,0]/255; %orange
        deadinfColor = [0.4 0.4 0.4]; % pink
        if round(4*t)>frameNum
            frameNum=round(4*t);
            h(1) = scatter(Grid(logical(x(1:4:end)),1),Grid(logical(x(1:4:end)),2),pointSize,liveColor,'filled')
            h(1).MarkerEdgeColor = liveColor;
            hold on;
            h(2) = scatter(Grid(logical(x(2:4:end)),1),Grid(logical(x(2:4:end)),2),pointSize,virusColor,'filled')
            h(2).MarkerEdgeColor = virusColor;
            h(3) = scatter(Grid(logical(x(3:4:end)),1),Grid(logical(x(3:4:end)),2),pointSize,deathColor,'filled')
            h(3).MarkerEdgeColor = deathColor;
            h(4) = scatter(Grid(logical(x(4:4:end)),1),Grid(logical(x(4:4:end)),2),pointSize,deathColor,'filled')
            h(4).MarkerEdgeColor = virusColor;
            for ii=1:4
                h(ii).LineWidth = lineSize;
            end
            
            
            axis equal;
            set(ax1,'xcolor','none','ycolor','none');

            tl = title(['T = ' sprintf('%.2f',frameNum/4) 'h']);
            tl.Position(2)=10.5;
            tl.Color = 'k';
            
            if t==0
                if ~isempty(arg.TNF)
                   an = annotation('textbox',[0.5, 0.7, 0.14, 0.05],'String',{['TNF=' num2str(arg.TNF) 'ng/ml']},'LineStyle','none', 'Fontsize', 24,'color','k');
                end
                
                hl = legend('   Healthy','   Infected','   Dead Bystander','   Dead post infection');
                hl.Location = 'northeastoutside';
                hl.Position = [0.63    0.4    0.3    0.2046];
                hl.Box = 'off';

                hl.FontSize = 14;
                for i=1:4;
                    hlobj(i).Color='w';
                    hlobj(i).FontSize=20;
                end
                for i=5:8;
                    hlobj(i).Children.MarkerSize=12;
                end

            end
            
            drawnow
            shg
            if arg.save
                frame = getframe(gcf);
                im = frame2im(frame);
                writeVideo(outputVideo,im);
            end
        end
    else
        if round(4*t)>frameNum
            frameNum=round(4*t);
            Dyn.infected = [Dyn.infected, numel(find(x(2:4:end)))];
            Dyn.infDead = [Dyn.infDead, numel(find(x(4:4:end)))];
            Dyn.t = [Dyn.t, t];
        end
    end
    
    %update k
    NNvLoad = neighMat*(v(:).*x(2:nSpecies:end));
    k(1:nReactions:end) = VI./(1+1./NNvLoad);  %infection \propto virus load of NN
    
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
        x = x+S(:,i);                                 % update the configuration
        if mod(i,3)==1
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

a = k(1:nReactions:end)./ infDeathRate;%k(3:nReactions:end);
a(a==0)=[];
R0 = mean(a);
if isnan(R0)
    R0=0;
end

if arg.save
    close(outputVideo);
end



