%% define colors
liveColor = [64,224,208]/255; %cyan
virusColor = [148, 0, 211]/255; %purple
deathColor = [255,165,0]/255; %orange
deadinfColor = [0.4 0.4 0.4]; % pink

cTNF = ([100*(1/2).^[0:8], 0]);

cmTNF = @(x) makeColorMap([0,100,0]/255,[152,251,152]/255,x);

cmMAC = @gray;
%%
   set(0,'DefaultTextInterpreter', 'tex')
set(0, 'DefaultAxesFontName', 'Arial')
set(0, 'DefaultUIControlFontName', 'Arial')
set(0,'defaulttextfontname','Arial');
set(groot,'defaultFigureColor','w')
set(groot,'defaultAxesColor','w')
set(groot,'defaultAxesTickLength',[0.03 0.01])
set(groot,'defaultLineLineWidth',2)

%% Experimental death rates for infected and uninfected cells taken from experiment TNFTitr_HighMOI_Dec122019_2019Dec12
infDeathRatesExp = [0.1481    0.1339    0.1354    0.1225    0.1224    0.1082    0.1137    0.1024    0.0867    0.0159];
uninfDeathRatesExp = [0.0343    0.0281    0.0240    0.0227    0.0154    0.0047    0.0011    0.0018    0.0015   -0.0023];

%timepoints
t = 0:1/3:48;
%TNF concentrations used
cTNF = ([100*(1/2).^[0:8], 0]);
%% This is how we calculate rates by interpolation of experimental results:
TNF = 100; %input TNF
infDeathRate = interp1(cTNF,filtfilt([1,1,1],3,infDeathRatesExp),TNF);
basalDeathRate = interp1(cTNF,filtfilt([1,1,1],3,infDeathRatesExp),TNF);

%% Make triangular grid
A = triangleGrid([-7 -7 7 7], [0,0], 1);
B = triangleGrid([-9 -9 9 9], [0,0], 1);%for presentation purposes

nCells = size(A,1);
%find all nearest neighnors and create sparse matrix
[IDX, D] = rangesearch(A, A,1.01);
IDX = cellfun(@(x,y) x(y>0), IDX,D,'uniformoutput',false);
pointsToMatch = arrayfun(@(x,y) repmat(x,y,1),1:numel(IDX), cellfun(@numel,IDX)','UniformOutput', false);
nnMatrix = sparse(cat(1,pointsToMatch{:}), cat(2,IDX{:})',1);

pointSize = 40;
lineSize = 2;
%% Init cells, center cell infected
infected = zeros(nCells,1);
infected(sum(A==0, 2)==2)=1;

x= [~infected, infected, zeros(nCells,1), zeros(nCells,1)]';
x = x(:);

figure('color','w','Position',[100, 100,600, 300])
ax1 = axes('position',[0.1,0.1,0.4,0.8]);
h1 = voronoi(B(:,1),B(:,2));
    h1(1).Color = 'none';
    h1(2).Color = [0 0 0 0.2];
    hold on
h(1) = scatter(A(logical(x(1:4:end)),1),A(logical(x(1:4:end)),2),pointSize,liveColor,'filled')
h(1).MarkerEdgeColor = liveColor;
hold on;
h(2) = scatter(A(logical(x(2:4:end)),1),A(logical(x(2:4:end)),2),pointSize,virusColor,'filled')
h(2).MarkerEdgeColor = virusColor;
h(3) = scatter(A(logical(x(3:4:end)),1),A(logical(x(3:4:end)),2),pointSize,deathColor,'filled')
h(3).MarkerEdgeColor = deathColor;
h(4) = scatter(A(logical(x(4:4:end)),1),A(logical(x(4:4:end)),2),pointSize,virusColor,'filled')
h(4).MarkerEdgeColor = deathColor;
axis equal
set(ax1,'xcolor','w','ycolor','w','xlim', [-7, 7],'ylim', [-7, 7])
shg


hl = legend(h, 'Healthy','Infected','False positive','Dead post infection');
hl.Box = 'off';
hl.Position=[0.5752    0.5755    0.3617    0.3117];
%% Run model once for example
close all
VGR = 0.01;%viral growth rate in units of infection ec50, See Table 1
VI = 1.5; %viral infectivity, See table 1

TNF = str2double(  inputdlg('Input TNF concentration in ng/ml'));

infDeathRate = max(interp1(cTNF,filtfilt([1,1,1],3,infDeathRatesExp),TNF),0);
basalDeathRate = max(interp1(cTNF,filtfilt([1,1,1],3,uninfDeathRatesExp),TNF),0);

[xOut, NNLoad, R0, Dyn] = SSA_forSpatialSIDGrids(x,nnMatrix,A,VGR,VI,basalDeathRate,infDeathRate);

sum(reshape(xOut,4,[]),2)./sum(sum(reshape(xOut,4,[]),2));

%
figure('color','w','Position',[100, 100,600, 300])
ax1 = axes('position',[0.1,0.1,0.4,0.8]);
h1 = voronoi(B(:,1),B(:,2));
    h1(1).Color = 'none';
    h1(2).Color = [0 0 0 0.2];
    hold on
h(1) = scatter(A(logical(xOut(1:4:end)),1),A(logical(xOut(1:4:end)),2),pointSize,liveColor,'filled')
h(1).MarkerEdgeColor = liveColor;
hold on;
h(2) = scatter(A(logical(xOut(2:4:end)),1),A(logical(xOut(2:4:end)),2),pointSize,virusColor,'filled')
h(2).MarkerEdgeColor = virusColor;
h(3) = scatter(A(logical(xOut(3:4:end)),1),A(logical(xOut(3:4:end)),2),pointSize,deathColor,'filled')
h(3).MarkerEdgeColor = deathColor;
h(4) = scatter(A(logical(xOut(4:4:end)),1),A(logical(xOut(4:4:end)),2),pointSize,virusColor,'filled')
h(4).MarkerEdgeColor = deathColor;
axis equal
set(ax1,'xcolor','w','ycolor','w','xlim', [-7, 7],'ylim', [-7, 7])
shg


an = annotation('textbox',[0.7, 0.1, 0.14, 0.4],'String',{['TNF=' num2str(TNF) 'ng/ml']},'LineStyle','none');

shg
hl = legend(h, 'Healthy','Infected','False positive','Dead post infection');
hl.Box = 'off';
hl.Position=[0.5752    0.5755    0.3617    0.3117];



%% Make triangular grid
A = triangleGrid([-7 -7 7 7], [0,0], 1);
B = triangleGrid([-9 -9 9 9], [0,0], 1);%for presentation purposes

nCells = size(A,1);
%find all nearest neighnors and create sparse matrix
[IDX, D] = rangesearch(A, A,1.01);
IDX = cellfun(@(x,y) x(y>0), IDX,D,'uniformoutput',false);
pointsToMatch = arrayfun(@(x,y) repmat(x,y,1),1:numel(IDX), cellfun(@numel,IDX)','UniformOutput', false);
nnMatrix = sparse(cat(1,pointsToMatch{:}), cat(2,IDX{:})',1);

pointSize = 40;
lineSize = 2;
%% Run model 500 times for each condition to get statistics
%% This takes ~20min running on 32 cores
VGR = 0.01;%viral growth rate, See Table 1, Figure S4C-D
VI = 1.5; %viral infectivity, See Table 1, Figure S4C-D

nRuns = 500;
TNF = [0 fliplr(100*(1/sqrt(2)).^[0:17])];
Stats = {}
xs = {}
R0s = {}
Dynamics = {}
for j=1:numel(TNF)
    infDeathRate = max(interp1(cTNF,filtfilt([1,1,1],3,infDeathRatesExp),TNF(j)),0);
    basalDeathRate = max(interp1(cTNF,filtfilt([1,1,1],3,uninfDeathRatesExp),TNF(j)),0);
    pStats = cell(nRuns,1);
    pxs = cell(nRuns,1);
    Dyns = cell(nRuns,1);

    parfor i=1:nRuns
        [pStats{i}, pxs{i} Rs{i}, Dyns{i}] = modelRunForStats_v2(nnMatrix,[],VGR,VI,basalDeathRate,infDeathRate, find(sum(A==0, 2)==2));
    end
    Stats{j} = cat(2,pStats{:})';
    xs{j} = cat(2,pxs{:})';
    R0s{j} = cat(2,Rs{:})';
    Dynamics{j} = cat(1,Dyns{:})';
end




%% Plot of average model results, spread and death - Figure 4.B
figure('color','w','Position',[100, 100,900, 600])
rang = [1 3 6];
for i=1:numel(rang)
    xOut = mean(xs{rang(i)});
    ax1 = axes('position',[0.03+0.31*(i-1),0.52,0.28,0.4]);

    h1 = voronoi(B(:,1),B(:,2));
    h1(1).Color = 'none';
    h1(2).Color = [0 0 0 0.2];
    hold on
    vCM = makeColorMap([1,1,1],virusColor,1001);
    h = scatter(A(:,1),A(:,2),150*(xOut(2:4:end)+xOut(4:4:end)).^2+0.000001,vCM(interp1(0:0.001:1,0:1000,(xOut(2:4:end))./max(xOut(2:4:end)+0.00001),'nearest')+1,:) ,'filled');%xOut(2:4:end)
    h.MarkerFaceAlpha = 1
    hold on;
    dCM = makeColorMap([1,1,1], [0 0 0],1001);
    h = scatter(A(:,1),A(:,2),150*xOut(4:4:end).^2+0.000001,dCM(interp1(0:0.001:1,0:1000,xOut(4:4:end)./max(xOut(4:4:end)+0.00001),'nearest')+1,:) ,'filled');%xOut(2:4:end)
    h.MarkerFaceAlpha = 1
%
     dCM = makeColorMap([1,1,1], liveColor,1001);
     h = scatter(A(:,1),A(:,2),150*xOut(1:4:end)+0.0001,dCM(interp1(0:0.001:1,0:1000,xOut(1:4:end)./max(xOut(1:4:end)+0.00001),'nearest')+1,:) ,'filled');%xOut(2:4:end)
     h.MarkerFaceAlpha = 0.5
    
    
    axis equal
    set(ax1,'xcolor','none','ycolor','none','XLim',[-7,7],'YLim',[-7,7]);
end

rang = [9 13 19];
for i=1:numel(rang)
    xOut = mean(xs{rang(i)});
    ax1 = axes('position',[0.03+0.31*(i-1),0.03,0.28,0.4]);

    h1 = voronoi(B(:,1),B(:,2));
    h1(1).Color = 'none';
    h1(2).Color = [0 0 0 0.2];
    hold on
    vCM = makeColorMap([1,1,1],virusColor,1001);
    h = scatter(A(:,1),A(:,2),150*(xOut(2:4:end)+xOut(4:4:end)).^2+0.000001,vCM(interp1(0:0.001:1,0:1000,(xOut(2:4:end))./max(xOut(2:4:end)+0.00001),'nearest')+1,:) ,'filled');%xOut(2:4:end)
    h.MarkerFaceAlpha = 1
    hold on;
    dCM = makeColorMap([1,1,1], [0 0 0],1001);
    h = scatter(A(:,1),A(:,2),150*xOut(4:4:end).^2+0.000001,dCM(interp1(0:0.001:1,0:1000,xOut(4:4:end)./max(xOut(4:4:end)+0.00001),'nearest')+1,:) ,'filled');%xOut(2:4:end)
    h.MarkerFaceAlpha = 1
% 
     dCM = makeColorMap([1,1,1], liveColor,1001);
     h = scatter(A(:,1),A(:,2),150*xOut(1:4:end)+0.0001,dCM(interp1(0:0.001:1,0:1000,xOut(1:4:end)./max(xOut(1:4:end)+0.00001),'nearest')+1,:) ,'filled');%xOut(2:4:end)
     h.MarkerFaceAlpha = 0.5
    
    axis equal
    set(ax1,'xcolor','none','ycolor','none','XLim',[-7,7],'YLim',[-7,7]);
end



%% Plot Model statistics - Figure 4.D
close all
LinearRange = 0.4

figure('color','w','Position',[100,100, 900, 450])
ax = axes('Position', [0.15, 0.14, 0.75/2, 0.75])
a = cellfun(@(x) x(:,1), Stats,'uniformOutput',false)
plot(asinh(TNF/LinearRange), median(cat(2,a{:})),'Color', liveColor);
set(gca,'YLim',[0,100])
ylabel('Percent of cells')
hold on
b = cellfun(@(x) x(:,2), Stats,'uniformOutput',false)
plot(asinh(TNF/LinearRange), median(cat(2,b{:})),'Color',virusColor);
set(gca,'YLim',[0,100])
hold on
c = cellfun(@(x) x(:,4), Stats,'uniformOutput',false)
plot(asinh(TNF/LinearRange), median(cat(2,c{:})),'Color',deadinfColor);
set(gca,'YLim',[0,100])


d = cellfun(@(x) x(:,3), Stats,'uniformOutput',false)
plot(asinh(TNF/LinearRange), median(cat(2,d{:})),'Color',deathColor);
set(gca,'YLim',[0,100])

hold on

h = boxplot(cat(2,a{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' },'Color', liveColor);
shg
xlabel('[TNF] ng/ml')
set(gca,'YLim',[-1,100]/100)
h = boxplot(cat(2,b{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' },'Color', virusColor);
shg
set(gca,'YLim',[-1,100]/100)


h = boxplot(cat(2,c{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' },'Color', deadinfColor);
shg
xlabel('[TNF] ng/ml')
set(gca,'YLim',[-1,100]/100)
h = boxplot(cat(2,d{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels',  {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' },'Color', deathColor);
shg
set(gca,'YLim',[-1,100]/100)

hleg = legend('   TN - Healthy','   FN - Infected','   TP - Dead infected','   FP - Dead uninfected');
hleg.Position = [0.5818    0.6425    0.2726    0.2046];
hleg.Box = 'off';

ax.YTick = 0:0.2:1;
ax.YTickLabel = 100*[0:0.2:1];

a1 = sort(- (1:9)'*10.^(1:2));
a2 =(1:9)'*10.^(1:5);
a = [sort(a1(:)); 0; a2(:)]/1000;
emptyStr = {''; ''; '';''; ''; '';''; ''};
ax.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
ax.XTick = asinh(a/LinearRange);

xl = xlabel('[TNF] ng/ml')
xl.Position(2) = xl.Position(2)


%% Infected over time means - Figure 4.F
tzeva = cmTNF(17);
figure('color','w','Position',[100,100, 450, 450])
ax = axes('Position', [0.2, 0.15, 0.7, 0.7])
LinearRangeY = 5;

for i=1:17
[XX, YY] = BinData_v3(cat(2,Dynamics{i}.t), asinh(cat(2,Dynamics{i}.infected)/LinearRangeY), -0.5:49)
h = plot(XX, YY);
h.Color = tzeva(i,:)
hold on
end
shg
xlabel('time (h)')
ylabel('mean # infected cells')
ax.XLim=[0, 48]

a1 = sort(- (1:9)'*10.^(1:2));
a2 =(1:9)'*10.^(1:5);
a = [sort(a1(:)); 0; a2(:)]/1000;
emptyStr = {''; ''; '';''; ''; '';''; ''};
XTickLabel = [emptyStr; {' '}; emptyStr; {' '}; '0'; {' '}; emptyStr; '  '; emptyStr; '1'; emptyStr; '10'; emptyStr; '100'; emptyStr];
YTicks = asinh(a/LinearRangeY);
set(gca,  'ytick', YTicks, 'yticklabel', XTickLabel)


%% infected cells over time stats and mean - Figure S4.G
figure('color','w','Position',[100,100, 450, 450])
ax = axes('Position', [0.22, 0.15, 0.7, 0.7])
LinearRangeY = 5;

i=1

for j=1:numel(Dynamics{i})
plot(cat(2,Dynamics{i}(j).t), asinh(cat(2,Dynamics{i}(j).infected)/LinearRangeY),'-','Color',[0.5,0.5,0.5], 'LineWidth', 0.25)
hold on
end
[XX, YY] = BinData_v3(cat(2,Dynamics{i}.t), cat(2,Dynamics{i}.infected), -0.5:49)
plot(XX, asinh(YY/LinearRangeY),'r');
shg
ax.XLim=[0, 48]
ax.YLim = asinh([0, 150]/LinearRangeY)
xlabel('time (h)')
ylabel('# infected cells')

a1 = sort(- (1:9)'*10.^(1:2));
a2 =(1:9)'*10.^(1:5);
a = [sort(a1(:)); 0; a2(:)]/1000;
emptyStr = {''; ''; '';''; ''; '';''; ''};
XTickLabel = [emptyStr; {' '}; emptyStr; {''}; '0'; {''}; emptyStr; ' '; emptyStr; '    1'; emptyStr; '   10'; emptyStr; '    100'; emptyStr];
YTicks = asinh(a/LinearRangeY);
set(gca,  'ytick', YTicks, 'yticklabel', XTickLabel)


%% Rt distribution - Figure 4.G
% Interpertation : a live cell is Rt times more likely of get infected by 
% a neighbor than for the infected neighbor to die. 
figure('color','w','Position',[100,100, 450, 450])
ax = axes('Position', [0.23, 0.15, 0.7, 0.7])
LinearRange = 0.4;

LinearRangeY = 1;

h1= violin(cellfun(@(x) asinh(x/1), R0s,'UniformOutput', false),'x',asinh(TNF/LinearRange),'support',[-1,10000],'bw',0.05,'facecolor', [0.6 0 0.6], 'medc',[0 0 0])

ylabel('$R_t = \frac{P_{\textrm{getting infected by your infected neighbor}}}{P_{\textrm{your infected neighbor will die}}}$','Interpreter','latex','Fontsize', 16)
xlabel('[TNF] (ng/ml)')

hold on
plot(asinh([-1, 200]./LinearRange), asinh([1,1]./LinearRangeY),'r')
a1 = sort(- (1:9)'*10.^(1:2));
a2 =(1:9)'*10.^(1:5);
a = [sort(a1(:)); 0; a2(:)]/1000;
emptyStr = {''; ''; '';''; ''; '';''; ''};
XTickLabel = [emptyStr; {' '}; emptyStr; {''}; '0'; {''}; emptyStr; ' '; emptyStr; '1'; emptyStr; '10'; emptyStr; '100'; emptyStr];
XTicks = asinh(a/LinearRange);
YTicks = asinh(a/LinearRangeY);

set(gca, 'xtick', XTicks, 'xticklabel', XTickLabel, 'ytick', YTicks, 'yticklabel', XTickLabel)

ax.XLim = asinh([-0.1, 150]/LinearRange)
ax.YLim = asinh([-0.1, 60]/LinearRangeY)


%thanks for reading


% %% Run model with variable VI and VGR - This takes A LONG time
% VGRarray = 0.01*2.^([-4:4]);%viral growth rate, See Table 1, Figure S4C-D
% VIarray = 1.5*2.^([-4:4]); %viral infectivity, See Table 1, Figure S4C-D
% modelResults=struct;
% for jj=1:numel(VGRarray)
%     VGR=VGRarray(jj);
%     for kk=1:numel(VIarray)
%         VI=VIarray(kk)
%         nRuns = 200;
%         TNF = [0 fliplr(100*(1/sqrt(2)).^[0:17])];
%         Stats = {}
%         xs = {}
%         R0s = {}
%         Dynamics = {}
%         for j=1:numel(TNF)
%             infDeathRate = max(interp1(cTNF,filtfilt([1,1,1],3,infDeathRatesExp),TNF(j)),0);
%             basalDeathRate = max(interp1(cTNF,filtfilt([1,1,1],3,uninfDeathRatesExp),TNF(j)),0);
%             pStats = cell(nRuns,1);
%             pxs = cell(nRuns,1);
%             Dyns = cell(nRuns,1);
%             
%             parfor i=1:nRuns
%                 [pStats{i}, pxs{i} Rs{i}, Dyns{i}] = modelRunForStats_v2(nnMatrix,[],VGR,VI,basalDeathRate,infDeathRate, find(sum(A==0, 2)==2));
%             end
%             Stats{j} = cat(2,pStats{:})';
%             xs{j} = cat(2,pxs{:})';
%             R0s{j} = cat(2,Rs{:})';
%             Dynamics{j} = cat(1,Dyns{:})';
%         end
%         
%         modelResults(jj, kk).VGR=VGR;
%         modelResults(jj, kk).VI=VI;
%         modelResults(jj, kk).Stats=Stats;
%         modelResults(jj, kk).xs=xs;
%         modelResults(jj, kk).R0s=R0s;
%         modelResults(jj, kk).Dynamics=Dynamics;
%     end
% end

%%
%load('/bigstore/GeneralStorage/Alon/Figures/DecisionPaper2019/paramScanModelResults.mat')

%% Sweep Compare figure
% m1s=nan(9,9);
% minds=nan(9,9);
% for i=1:9
%     for j=1:9
%     Stats  = modelResults(i,j).Stats;
%     a = cellfun(@(x) x(:,1), Stats,'uniformOutput',false)
%     [m1,mind]=max(median(cat(2,a{:})));
%     m1s(i,j)=m1;
%     minds(i,j)=mind;
%     end
% end
% 
% %%
% VGRarray = 0.01*2.^([-4:4]);%viral growth rate, See Table 1, Figure S4C-D
% VIarray = 1.5*2.^([-4:4]); %viral infectivity, See Table 1, Figure S4C-D
% 
% 
% figure('color','w','Position',[100,100, 900, 450])
% ax = axes('Position', [0.1, 0.14, 0.6/2, 0.6])
% as=[]
% for i=1:9
%     Stats  = modelResults(5,i).Stats;
%     a =cellfun(@(x) x(:,1), Stats,'uniformOutput',false);
%     as = [as; 100*median(cat(2,a{:}))];
% end
% 
% h = imagesc(as)
% ylabel('Viral infectivity (h^{-1})')
% xlabel('[TNF] (ng/ml)')
% 
% ax.YDir='normal'
% %h.EdgeColor='none'
% colormap('inferno')
% 
% r1 = rectangle('Position',[0.5, 4.5, 19, 1],'EdgeColor',[0 0.7 0],'LineWidth', 1.5)
% set(gca,'XTick', [1:4:numel(TNF) 19],'XTickLabel', [round(100*TNF(1:4:end))./100, 100])
% set(gca,'YTick', 1:2:numel(VIarray),'YTickLabel', round(100*VIarray(1:2:end))./100,'YDir', 'normal')
% 
% 
% 
% 
% ax = axes('Position', [0.54, 0.14, 0.6/2, 0.6])
% as=[]
% for i=1:9
%     Stats  = modelResults(i,5).Stats;
%     a =cellfun(@(x) x(:,1), Stats,'uniformOutput',false);
%     as = [as; 100*median(cat(2,a{:}))];
% end
% 
% h = imagesc(as)
% ylabel('Viral growth rate (a.u./h)')
% xlabel('[TNF] (ng/ml)')
% 
% ax.YDir='normal'
% %h.EdgeColor='none'
% colormap('inferno')
% 
% r1 = rectangle('Position',[0.5, 4.5, 19, 1],'EdgeColor',[0 0.7 0],'LineWidth', 1.5)
% set(gca,'XTick', [1:4:numel(TNF) 19],'XTickLabel', [round(100*TNF(1:4:end))./100, 100])
% set(gca,'YTick', 1:2:numel(VGRarray),'YTickLabel', round(1000*VGRarray(1:2:end))./1000,'YDir', 'normal')
% 
% cb = colorbar('Position',[0.88, 0.14, 0.03, 0.6])
% cb.Ticks=[0:25:100]
% cb.Limits=[0, 100]
% cb.Label.String = '% Healthy'
% cb.Label.FontSize=14
% cb.Label.Position(1) = 0.2
% cb.Label.Color='w'
% 
% %%
% figname = 'Fig_ModelVI_VGRSweepCompare'
% set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off','Renderer', 'painters')
% print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
% print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);



%% Run model with variable VI - Go get some coffee and come back in a few hours/days
VGRarray = 0.01;%viral growth rate, See Table 1, Figure S4C-D
VIarray = 1.5*2.^([-10:10]./2); %viral infectivity, See Table 1, Figure S4C-D
modelResults2=struct;
for jj=1:numel(VGRarray)
    VGR=VGRarray(jj);
    for kk=1:numel(VIarray)
        VI=VIarray(kk)
        nRuns = 200;
        TNF = [0 fliplr(100*(1/sqrt(2)).^[0:17])];
        Stats = {}
        xs = {}
        R0s = {}
        Dynamics = {}
        for j=1:numel(TNF)
            infDeathRate = max(interp1(cTNF,filtfilt([1,1,1],3,infDeathRatesExp),TNF(j)),0);
            basalDeathRate = max(interp1(cTNF,filtfilt([1,1,1],3,uninfDeathRatesExp),TNF(j)),0);
            pStats = cell(nRuns,1);
            pxs = cell(nRuns,1);
            Dyns = cell(nRuns,1);
            
            parfor i=1:nRuns
                [pStats{i}, pxs{i} Rs{i}, Dyns{i}] = modelRunForStats_v2(nnMatrix,[],VGR,VI,basalDeathRate,infDeathRate, find(sum(A==0, 2)==2));
            end
            Stats{j} = cat(2,pStats{:})';
            xs{j} = cat(2,pxs{:})';
            R0s{j} = cat(2,Rs{:})';
            Dynamics{j} = cat(1,Dyns{:})';
        end
        
        modelResults2(jj, kk).VGR=VGR;
        modelResults2(jj, kk).VI=VI;
        modelResults2(jj, kk).Stats=Stats;
        modelResults2(jj, kk).xs=xs;
        modelResults2(jj, kk).R0s=R0s;
        modelResults2(jj, kk).Dynamics=Dynamics;
    end
end
%%
%load('/bigstore/GeneralStorage/Alon/Figures/DecisionPaper2019/paramScanModelResults2.mat')

%% Fig 4I - VI Sweep
VGRarray = 0.01;%viral growth rate, See Table 1, Figure S4C-D
VIarray = 1.5*2.^([-10:10]./2); %viral infectivity, See Table 1, Figure S4C-D

LinearRange = 0.4

figure('color','w','Position',[100,100, 900, 450])
ax = axes('Position', [0.15, 0.14, 0.75/2, 0.75])
as=[]
for i=1:21
    Stats  = modelResults2(1,i).Stats;
    a =cellfun(@(x) x(:,1), Stats,'uniformOutput',false);
    as = [as; 100*median(cat(2,a{:}))];
end
%plot(asinh(TNF/LinearRange), as,'Color', liveColor);

LinearRange=0.4;
LinearRangeY=0.1;

h = imagesc(as)
%[x, y] = meshgrid(asinh(TNF/LinearRange), asinh(VIarray/LinearRangeY));
%z = zeros(size(y))
%h = surface(x,y,z,255*as./100,'EdgeColor','none','FaceColor','texturemap',...
%  'CDataMapping','direct');
axis tight

ylabel('Viral infectivity (h^{-1})')
xlabel('[TNF] (ng/ml)')
ax.YDir='normal'
%h.EdgeColor='none'
colormap('inferno')



set(gca,'XTick', [1:4:numel(TNF) 19],'XTickLabel', [round(100*TNF(1:4:end))./100, 100])
set(gca,'YTick', 1:4:numel(VIarray),'YTickLabel', round(100*VIarray(1:4:end))./100,'YDir', 'normal')



r1 = rectangle('Position',[0.5, 10.5, 19, 1],'EdgeColor',[0 0.7 0],'LineWidth', 1.5)


cb = colorbar('Position',[0.58, 0.14, 0.03, 0.75])

cb.Ticks=[0:25:100]
cb.Limits=[0, 100]
cb.Label.String = '% Healthy'
cb.Label.FontSize=14
cb.Label.Position(1) = 0.2
cb.Label.Color='w'


a1 = annotation('arrow',[0.5450, 0.5250], [0.14+0.75/2, 0.14+0.75/2],'color',[0 0.7 0]);
t1 = text([20.8], [9.5 ], 'measured','FontSize', 11,'EdgeColor', 'none','Rotation', 90, 'color', [0 0.7 0])

%%
figname = 'Fig_ModelVISweep'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off','Renderer', 'painters')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);
