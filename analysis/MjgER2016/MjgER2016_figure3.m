% MjgER2016 Figure3
%
% HP:= Head Pitch
% BP:= Body Pitch
% PFD:= behavior field restricted to DRZ[-0.5,0.5]
%
% Subplots:
%    A. BHV pfd example, large for clear axes 
%        1. PFD behavioral space HPxBP
%        2. PFS theta
%    B. Place field examples
%        1. auto correlogram of unit
%        2. PFD behavioral space HPxBP
%        4. PFS theta
%        5. PFS rear & theta
%        6. PFS hloc & theta
%        7. PFS lloc & theta
%        8. PFS hpause & theta
%        9. PFS lpause & theta
%    C. Eigenvectors of PFD decomposition of HPxBP space
%    D. t-SNE mapping of fscores within HPxBP of first 3 eigenvectors
%    E. ECDF of zscores


MjgER2016_load_data();
%  MjgER2016_load_data:
%
%  Variables:
%      Trials      units      cluSessionMap      pitchReferenceTrial      FigDir
%      sessionListName        sessionList        states                   numStates      
%      interpParPfsp          interpParDfs
%      
%  Functions:
%      reshape_eigen_vector

MjgER2016_figure3_args();
% MjgER2016_figure3_args:
%
% Default arument overrides:
%     fet_HB_pitchB
%     compute_bhv_ratemaps
%     compute_bhv_ratemaps_shuffled
%     compute_bhv_ratemap_erpPCA

% EXAMPLE Trial : jg05-20120312
tind = 20;
Trial = Trials{tind};
stc = Trials{tind}.stc.copy();


% LOAD place restricted behavior fields
bfs   = cf(@(t,u)   compute_bhv_ratemaps(t,u),  Trials, units);
bfsEx = bfs{tind};
bfsTag = {'HPITCHxBPITCH'};


% COMPUTE bfs erpPCA
[eigVecs, eigScrs, eigVars, unitSubset, validDims, zrmMean, zrmStd] = ...
                    compute_bhv_ratemaps_erpPCA(bfs, units);
numComp = size(eigVecs,2);
fpc  = cell([1,numComp]);
for i = 1:numComp,
    fpc{i} = nan(size(validDims));
    fpc{i}(validDims) = eigVecs(:,i);
end
fpcLims = [min(cellfun(@min,fpc)),max(cellfun(@max,fpc))];


% LOAD theta place fields
pft = pfs_2d_theta(Trial);
pft.parameters.states = 'theta-groom-sit';

% LOAD behavior place field
pfs = cat(2,{pft},pfs_2d_states(Trial));
% SORT place field states to match states
pfStates = cf(@(s) ['^',s,'$'],cf(@(p) p.parameters.states,pfs));
for s = 1:numel(pfStates),
    pfStates{s} = strrep(pfStates{s},'&','[&]');
    %pfStates{s} = strrep(pfStates{s},'-','[-]');
end    
for s = 1:numStates,
    psi(s) = find(~cellfun(@isempty,regexp(repmat(states(s),size(pfStates)),pfStates)));
end
pfs = pfs(psi);





bins = ;
% LOAD Behavioral state contours
[stateContourMaps,stateContourHandles] =                           ...
    bhv_contours(sessionListName,                                  ... sessionListName
                 'fet_HB_pitchB',                                  ... featureSet
                 [1,2],                                            ... featureInd
                 {linspace(-2,2,50),linspace(-2,2,50)},            ... featureBin
                 'Ed05-20140529.ont.all',                          ... referenceTrial
                 {'lloc+lpause&theta','hloc+hpause&theta',         ... states
                   'rear&theta'},                                  ...
                 'wcr',                                            ... stateColors
                 '0c641aa30eb6c126a115c42cdfc3fada'                ... tag
);





[fsrcz,FSrC,rmaps,FSCFr,FSrM,FSrS,fsrsMean,fsrsStd,rmapsShuffledMean,rmapsShuffled] = ...
    compute_bhv_ratemaps_erpPCA_scores('tag','68536b9a930748fd25a91dfc8ca84d39');


% GET auxiliary features
clu =  cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:), bfs,units);    
si  =  cf(@(p,u) p.data.si(:,ismember(p.data.clu,u),:),  bfs,units);
tlu =  cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
si  = cat(2, si{:});   
clu = cat(2, clu{:});    
tlu = cat(2, tlu{:});    
clu = [tlu',clu'];


[~,rind] = sortrows(clu);
si  = si(:,rind);

sigUnits = any(abs(fsrcz(:,1:3))>=1.96,2);

cc = eigScrs(:,[2,1,3])+0.75;
cc(~sigUnits,:) = repmat([0.75,0.75,0.75],[sum(~sigUnits),1]);


if ~exist('mapa','var'),
mapa = tsne([FSrC(:,1:3),(si(unitSubset)'-mean(si(unitSubset))')./std(si(unitSubset))'],[],2,3,40);    
% $$$     mapa = tsne([fsrcz(:,1:3),(si(unitSubset)'-mean(si(unitSubset))')./std(si(unitSubset))'],[],2,4,80);
% $$$     mapa = tsne([fsrcz(:,1:3),si(unitSubset)'],[],2,3,50);
end


% LOAD accgs
[accg,tbins] = autoccg(Trial);
cond_round = @(rate) max([round(rate,0),round(rate,1)].*[rate>=10,rate<10]);
nanColor = [0.15,0.15,0.15];

pageWidth  = 21.;
pageHeight = 29.7;
pwidth = 1.25;
pheight = 1.25;
xpad = 0.1;
ypad = 0.1;

xpos = 3.5:(pwidth+xpad):pageWidth;
ypos = fliplr(0.5:(pheight+xpad):pageHeight-3.7);

% SET figure opts
hfig = figure(666003);
clf();
% $$$ [hfig, figOpts, fax, sax] = set_figure_layout(hfig,'A4','portrait','centimeters',...

hfig.Units = 'centimeters';
hfig.Position = [1, 1, pageWidth,pageHeight];
hfig.PaperPositionMode = 'auto';



%% MjgER2016F3A - behavior and place field examples %%
cluMap = [20,119];
yind = 1;
xind = 1;



u = cluMap';
% SET color scale max
maxPfsRate = max(cell2mat(cf(@(p,u) maxRate(p,u,false,'prctile99',0.5),...
                             {pfs{1},bfsEx},repmat({u(2)},[1,2]))));
sp = gobjects([1,0]);
fax = axes('Position',[0,0,1,1],'Visible','off','Units','centimeters');
xlim([0,hfig.Position(3)]);
ylim([0,hfig.Position(4)]);

% ACCG 
% $$$ sp(end+1) = axes('Units','centimeters',...
% $$$                  'Position',[xpos(xind),ypos(yind)+pheight+ypad,pwidth,pheight],...
% $$$                  'FontSize', 8);
% $$$ bar(tbins,accg(:,u(2)));
% $$$ axis('tight');


xind = 1;


sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind),ypos(yind),pwidth*1.25,pheight*1.25],...
                 'FontSize', 8);
%                 'Position',[xpos(xind),ypos(yind),(pwidth+xpad)*2-xpad,(pheight+ypad)*2-ypad],...
bfs{tind}.plot(u(2),'mean',false,[0,maxPfsRate],false,0.5,false,[],@jet,[],nanColor);
text(-2,-0.45,num2str(cond_round(bfs{tind}.maxRate(u(2),false,1))),'FontSize',10,'Color',[1,1,1]);
xlabel({'HP (rad)'});
ylabel({'BP (rad)'});
xind = xind+2;
xlim([-2,nonzeros(xlim.*[0,1])-0.2])
sp(end).Color = nanColor;
title({'Behavior','Rate Map'});    



% PLOT theta example    
sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind),ypos(yind),pwidth*1.25,pheight*1.25],...
                 'FontSize', 8);
%                 'Position',[xpos(xind),ypos(yind),(pwidth+xpad)*2-xpad,(pheight+ypad)*2-ypad],...
plot(pfs{1},u(2),'mean',false,[0,maxPfsRate],true,0.5,false,interpParPfs,@jet,[],nanColor);
text(-490,-380,num2str(cond_round(maxPfsRate)),'FontSize',10,'Color',[1,1,1]);

hold('on');
rmap = plot(pfs{1},u(2),'mean',false,[0,maxPfsRate],true,0.5,false,interpParPfs);
[~,mxp] = pft.maxRate(u(2));
binGrid = cell([1,2]);
[binGrid{:}] = ndgrid(interpParPfs.bins{:});
binGrid = cat(3,binGrid{:});
binDist = sqrt(sum((binGrid-repmat(permute(mxp,[1,3,2]),[cellfun(@numel,interpParPfs.bins),1])).^2,3));
[cont,conHax] = contour(binGrid(:,:,1),binGrid(:,:,2),(rmap./max(rmap(:)))>0.2&binDist<250,[0.5,0.5],'m','LineWidth',1);
[contYMax,contYMaxInd] = max(reshape(binGrid(:,:,2).*double((rmap./max(rmap(:)))>0.2&binDist<250),[],1));
[contYMin,contYMinInd] = min(reshape(binGrid(:,:,2).*double((rmap./max(rmap(:)))>0.2&binDist<250),[],1));
contMaxInd = zeros([1,2]);
[contMaxInd(1),contMaxInd(2)] = ind2sub(cellfun(@numel,interpParPfs.bins),contYMaxInd);
contMinInd = zeros([1,2]);
[contMinInd(1),contMinInd(2)] = ind2sub(cellfun(@numel,interpParPfs.bins),contYMinInd);
contMaxPos = sq(binGrid(contMaxInd(1),contMaxInd(2),:));
contMinPos = sq(binGrid(contMinInd(1),contMinInd(2),:));

line(fax,...
     [sp(end-1).Position(1)+pwidth.*1.25,sp(end).Position(1)+(contMaxInd(1)./numel(interpParPfs.bins{1})).*pwidth.*1.25],...
     [sp(end-1).Position(2)+pheight.*1.25,sp(end).Position(2)+(contMaxInd(2)./numel(interpParPfs.bins{2})).*pheight.*1.25],...
     'Color','m','LineWidth',1);
line(fax,...
     [sp(end-1).Position(1)+pwidth.*1.25,sp(end).Position(1)+(contMinInd(1)./numel(interpParPfs.bins{1})).*pwidth.*1.25],...
     [sp(end-1).Position(2),sp(end).Position(2)+(contMinInd(2)./numel(interpParPfs.bins{2})).*pheight.*1.25],...
     'Color','m','LineWidth',1);
line(fax,...
     [sp(end-1).Position(1),sp(end-1).Position(1)+pwidth.*1.25],...
     [sp(end-1).Position(2),sp(end-1).Position(2)],...
     'Color','m','LineWidth',1);      
line(fax,...
     [sp(end-1).Position(1),sp(end-1).Position(1)+pwidth.*1.25],...
     [sp(end-1).Position(2)+pheight.*1.25,sp(end-1).Position(2)+pheight.*1.25],...
     'Color','m','LineWidth',1);      
line(fax,...
     [sp(end-1).Position(1)+pwidth.*1.25,sp(end-1).Position(1)+pwidth.*1.25],...
     [sp(end-1).Position(2),sp(end-1).Position(2)+pheight.*1.25],...
     'Color','m','LineWidth',1);      
line(fax,...
     [sp(end-1).Position(1),sp(end-1).Position(1)],...
     [sp(end-1).Position(2),sp(end-1).Position(2)+pheight.*1.25],...
     'Color','m','LineWidth',1);      
uistack(fax,'top');

title({'Spatial','Rate Map'});   
chax = colorbar(sp(end));
chax.LineWidth = 1;
chax.Units = 'centimeters';
colormap(sp(end),'jet');
caxis([0,maxPfsRate]);
chax.Position(1) = sp(end).Position(1)+pwidth*1.25+0.1;

sp(end).YTickLabel = {};
sp(end).XTickLabel = {};
line(fax,xpos(xind)+[0,pwidth*1.25/2],ypos(yind).*[1,1]-0.125,'Color','k','LineWidth',2);
text(fax,xpos(xind),ypos(yind)-0.35,'50cm','FontSize',8,'Color',[0,0,0]);

xind = xind+3;





%% MjgER2016F3B - behavior and place field examples %%

% $$$ % Multisession examples
% $$$ % selected units session x unit id
% $$$ clumap = [17,147;...
% $$$            3,171;...
% $$$            3,197;...
% $$$           16, 51;...
% $$$            1, 15;...
% $$$           17,165;...
% $$$           18,129;...
% $$$          ];


cluMap = [20,74;...
          20,73;...
          20,34;...
          20,83;...
          20,79;...
          20,103;...
          20,104;...
          20,109;...
          20,59];

% {accg,BFS{HP,BP},BFS{HP,BS},theta,rear,hloc,lloc,hpause,lpause}




labels = {'HPxBP','Theta','Rear','H Loc','L Loc','H Pause','L Pause'};

% PLOT 
%clf();    
yind = 3;
yinit = yind;
for u = cluMap',

    xind = 1;

% SET color scale max
    maxPfsRate = max(cell2mat(cf(@(p,u) ...
                                 maxRate(p,u,false,'prctile99',0.5),...
                                 [pfs,{bfsEx}],repmat({u(2)},[1,numel(pfs)+1]))));

    pfsMaxRatesMean = cell2mat(cf(@(p,u) ...
                                  max(p.maxRate(u,true,'mean')),...
                                  pfs,repmat({u(2)},[1,numel(pfs)])));
    %bfsMaxRatesMean = cell2mat(cf(@(p,u) max(p.maxRate(u,false,'mean')),bfs,repmat({u(2)},[1,numel(bfs)])));
    bfsMaxRatesMean = mean(bfsEx.maxRate(u,false,1));

    

% DRZFIELDS 
    sp(end+1) = axes('Units','centimeters',...
                     'Position',[xpos(xind),ypos(yind),pwidth,pheight],...
                     'FontSize', 8);

    bfsEx.plot(u(2),'mean',false,[0,maxPfsRate],false,0.5,false,[],@jet,[],nanColor);
    sp(end).YTickLabel = {};
    sp(end).XTickLabel = {};        
    text(-2,-0.475,num2str(cond_round(bfsMaxRatesMean)),'FontSize',8,'Color',[1,1,1]);
    xlim([-2,nonzeros(xlim.*[0,1])-0.2])
    sp(end).Color = nanColor;

    if yind == yinit, title(labels{xind});end        
    xind = xind+1;

    
    for s = 1:numStates,
% PLACEFIELDS MTAApfs
        sp(end+1) = axes('Units','centimeters',...
                         'Position',[xpos(xind),ypos(yind),pwidth,pheight],...
                         'FontSize', 8);

        plot(pfs{s},u(2),'mean',false,[0,maxPfsRate],true,0.5,false,interpParPfs,@jet,[],nanColor);
        sp(end).YTickLabel = {};
        sp(end).XTickLabel = {};
        if yind == yinit, title(labels{xind});end        
        xind = xind+1;
        text(-495,-380,num2str(cond_round(pfsMaxRatesMean(s))),'FontSize',8,'Color',[1,1,1]);
    end
    yind = yind+1;

    
end



%% MjgER2016F3C - behavior field erpPCA eigenvectors
% see: req20180123_vis_HBPITCHxBPITCH_erpPCA.m


% PLOT 
yind = 15
for i = 1:3,
    sp(end+1) = axes('Units','centimeters',...
                     'Position',[xpos(1),ypos(yind)+0.5,pwidth,pheight],...
                     'FontSize', 8);    
    
    imagescnan({bfsEx.adata.bins{:},abs(reshape_eigen_vector(fpc{i},bfsEx))},...
               fpcLims,'linear',false,nanColor,1,1);                % PRINT eigenvectors
    axis('xy');
    axis('tight');
    hold('on');    
    for s = 1:numel(stateContourHandles),                            % OVERLAY state Contours
        copyobj(stateContourHandles{s},sp(end));
    end
    sp(end).YTickLabel = {};
    sp(end).XTickLabel = {};
    sp(end).Color = nanColor;
    ylabel(['F',num2str(i)])
    xlim(bfsEx.adata.bins{1}([1,end]));
    xlim([-2,nonzeros(xlim.*[0,1])-0.2])            
    ylim(bfsEx.adata.bins{2}([1,end]));

    yind = yind + 1;
end



af(@(ax) set(ax,'LineWidth',1), sp);

%% MjgER2016F3D t-SNE mapping of fscores within HPxBP of first 3 eigenvectors
% see: req20180319.m

% CREATED Vars
% 
% NON - ESSENTIAL 
%     clu    - cluster ids
%     tlu    - trial ids
%     rind   - maps cluster order from MTAApfs object to Units order
%     D      - covariance matrix
%     LR     - eigenvectors
%     FSCFr  - pseudo-inverse of LR
%     rsMean   - rmaps shuffled mean
%     rsStd    - rmaps shuffled std
% 
% 
% AUXILLARAY vars
%     pfdShuffled
%     rmapsShuffledMean
%     rmapsShuffled
%     FSrM     - mean of 'shuffled' scores
%     FSrS     - std of 'shuffled' score
%     fsrsMean - 
%     fsrsStd
%     fsrsMean - 
%     fsrsStd
%
%
% OUTPUT Vars
%     rmaps    - matrix[D x S](numeric); rate maps corresponding to the valid eigenvector dims
%     FSrC     - matrix[U x V](Numeric); fscores
%     fsrcz    - matrix[U x V](Numeric); normalized fscores




% $$$ dsm = pdist([FSrC(:,1:3),si(unitSubset)']);
% $$$ mapa = mdscale(dsm,2);
% $$$ figure,plot(mapa(:,1),mapa(:,2),'.');


yind = yind - 1;
xind = 2;
sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind)+0.5,ypos(yind)+0.5,(pwidth+0.1)*3-0.1,(pheight+0.1)*3-0.1],...
                 'FontSize', 8);    

% $$$ figure,scatter(mapa(:,1),mapa(:,2),5,cc,'filled');

% $$$ ss = ones([size(cc,1),1])*10;
% $$$ ss(all(abs(fsrcz(:,1:3))>1.96,2)) = 20;
% $$$ ss(all(abs(fsrcz(:,1:3))<1.96,2)) = 5;
%mapa = tsne([FSrC(:,1:3),(si(unitSubset)'-mean(si(unitSubset))')./std(si(unitSubset))'],[],2,3,40);
%figure();

cla();
mi = [2,1];
%scatter3(FSrC(:,1),FSrC(:,3),FSrC(:,2),15,cc,'o','filled'); 
scatter(mapa(:,mi(1)),mapa(:,mi(2)),15,cc,'o','filled'); 
sp(end).YTickLabel = {};
sp(end).XTickLabel = {};
box('on');
axis('tight')
xlim(xlim.*1.05);
ylim(ylim.*1.05);

cluSessionSubset = cluSessionMap(unitSubset,:);
for u = cluMap'
    uind = find(ismember(cluSessionSubset,u','rows'));
    sigUnits(uind);
end
set(gca(),'Color',nanColor)
    
% $$$ sp(end+1)=subplot(359); hold('on');scatter(mapz(:,1)/10,mapz(:,2)/10,10,cc,'o','filled'); grid('on');
% $$$ title('tsne on zscores')

%plot(mapa(:,1),mapa(:,2),'.');
%mapz = tsne([fsrcz(:,1:3),si(unitSubsets{pfindex})'],[],2,4,25);
%figure,scatter(mapz(:,1),mapz(:,2),10,cc,'filled');


yind = yind-1;
xind = 6;
sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind),ypos(yind),(pwidth+0.1)*3-0.1,(pheight+0.1)*2-0.1],...
                 'FontSize', 8);    

hold('on');
for i = [2,1,3],
    [F,X] = ecdf(fsrcz(:,i));
    cdfplot(X);
end
lax = legend({'rear','low prone','high prone'},'Location','southeastoutside');
lax.Units = 'centimeters';
lax.Position(1) = 0.1;
drawnow();
lax.Position(1) = sum(sp(end).Position([1,3]));
Lines(-1.96,[],'k');
Lines(1.96,[],'k');

Lines(-3.1,[],'k');
Lines(3.1,[],'k');
xlim([-20,20])

