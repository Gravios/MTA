% MjgER2016 Figure3
%
% HP:= Head Pitch
% BP:= Body Pitch
% PFD:= behavior field restricted to DRZ[-0.5,0.5]
%
% Subplots:
%    A. BHV pfd example, large for clear axes 
%        1. auto correlogram of unit
%        2. PFD behavioral space HPxBP
%        3. PFD behavioral space BSxHS
%        4. PFS theta
%    B. Place field examples
%        1. auto correlogram of unit
%        2. PFD behavioral space HPxBP
%        3. PFD behavioral space BSxHS
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
%      Trials
%      units
%      cluSessionMap
%      pitchReferenceTrial
%      FigDir
%      sessionListName
%      sessionList
%      states
%      numStates
%      interpParPfsp
%      interpParDfs
%
%  Functions:
%      reshape_eigen_vector


Trial = Trials{20};
stc = Trials{20}.stc.copy();
% LOAD place fields
pft = pfs_2d_theta(Trial);
pft.parameters.states = 'theta-groom-sit';
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

% LOAD DRZ fields
dfs = req20180123_ver5(Trial);
dfst = {'HPITCHxBPITCH'};%,'HPITCHxBSPEED','BPITCHxBSPEED','BPITCHxHSPEED','HPITCHxRHM'};

[ddz] = compute_ddz(Trial,units{20},pft);

% LOAD accgs
[accg,tbins] = autoccg(Trial);

cond_round = @(rate) max([round(rate,0),round(rate,1)].*[rate>=10,rate<10]);


pageWidth  = 21.;
pageHeight = 29.7;

pwidth = 1.25;
pheight = 1.25;

xpad = 0.1;
ypad = 0.1;

xpos = 3.5:(pwidth+xpad):pageWidth;
ypos = fliplr(0.5:(pheight+xpad):pageHeight-3.7);

% SET figure opts
hfig = figure(666001);
hfig.Units = 'centimeters';
hfig.Position = [1, 1, pageWidth,pageHeight];
hfig.PaperPositionMode = 'auto';

clf();

%% MjgER2016F3A - behavior and place field examples %%


cluMap = [20,119];


yind = 1;
xind = 1;

% SET color scale max
u = cluMap';
maxPfsRate = max(cell2mat(cf(@(p,u) maxRate(p,u,false,'prctile99',0.5),...
                             [pfs(1),dfs],repmat({u(2)},[1,1+numel(dfs)]))));

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


xind = 2;


sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind),ypos(yind),pwidth*1.25,pheight*1.25],...
                 'FontSize', 8);
%                 'Position',[xpos(xind),ypos(yind),(pwidth+xpad)*2-xpad,(pheight+ypad)*2-ypad],...
dfs{1}.plot(u(2),'mean',false,[0,maxPfsRate],false,0.5,false,[],@jet);
text(-2,-0.45,num2str(cond_round(dfs{1}.maxRate(u(2),false,1))),'FontSize',8,'Color',[1,1,1]);
xlabel({'HP (rad)'});
ylabel({'BP (rad)'});
xind = xind+3;
xlim([-2,nonzeros(xlim.*[0,1])-0.2])
sp(end).Color = [0,0,0];
title('Bhv Rate Map');    



% PLOT theta example    
sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind),ypos(yind),pwidth*1.25,pheight*1.25],...
                 'FontSize', 8);
%                 'Position',[xpos(xind),ypos(yind),(pwidth+xpad)*2-xpad,(pheight+ypad)*2-ypad],...
plot(pfs{1},u(2),'mean',false,[0,maxPfsRate],true,0.5,false,interpPar,@jet);
text(-490,-380,num2str(cond_round(maxPfsRate)),'FontSize',10,'Color',[1,1,1]);

hold('on');
rmap = plot(pfs{1},u(2),'mean',false,[0,maxPfsRate],true,0.5,false,interpPar,@jet);
[~,mxp] = pft.maxRate(u(2));
binGrid = cell([1,2]);
[binGrid{:}] = ndgrid(interpPar.bins{:});
binGrid = cat(3,binGrid{:});
binDist = sqrt(sum((binGrid-repmat(permute(mxp,[1,3,2]),[cellfun(@numel,interpPar.bins),1])).^2,3));
[cont,conHax] = contour(binGrid(:,:,1),binGrid(:,:,2),(rmap./max(rmap(:)))>0.2&binDist<250,[0.5,0.5],'m','LineWidth',1);
[contYMax,contYMaxInd] = max(reshape(binGrid(:,:,2).*double((rmap./max(rmap(:)))>0.2&binDist<250),[],1));
[contYMin,contYMinInd] = min(reshape(binGrid(:,:,2).*double((rmap./max(rmap(:)))>0.2&binDist<250),[],1));
contMaxInd = zeros([1,2]);
[contMaxInd(1),contMaxInd(2)] = ind2sub(cellfun(@numel,interpPar.bins),contYMaxInd);
contMinInd = zeros([1,2]);
[contMinInd(1),contMinInd(2)] = ind2sub(cellfun(@numel,interpPar.bins),contYMinInd);
contMaxPos = sq(binGrid(contMaxInd(1),contMaxInd(2),:));
contMinPos = sq(binGrid(contMinInd(1),contMinInd(2),:));

line(fax,...
     [sp(end-1).Position(1)+pwidth.*1.5,sp(end).Position(1)+(contMaxInd(1)./numel(interpPar.bins{1})).*pwidth.*1.5],...
     [sp(end-1).Position(2)+pheight.*1.5,sp(end).Position(2)+(contMaxInd(2)./numel(interpPar.bins{2})).*pheight.*1.5],...
     'Color','m','LineWidth',1);
line(fax,...
     [sp(end-1).Position(1)+pwidth.*1.5,sp(end).Position(1)+(contMinInd(1)./numel(interpPar.bins{1})).*pwidth.*1.5],...
     [sp(end-1).Position(2),sp(end).Position(2)+(contMinInd(2)./numel(interpPar.bins{2})).*pheight.*1.5],...
     'Color','m','LineWidth',1);
line(fax,...
     [sp(end-1).Position(1),sp(end-1).Position(1)+pwidth.*1.5],...
     [sp(end-1).Position(2),sp(end-1).Position(2)],...
     'Color','m','LineWidth',1);      
line(fax,...
     [sp(end-1).Position(1),sp(end-1).Position(1)+pwidth.*1.5],...
     [sp(end-1).Position(2)+pheight.*1.5,sp(end-1).Position(2)+pheight.*1.5],...
     'Color','m','LineWidth',1);      
line(fax,...
     [sp(end-1).Position(1)+pwidth.*1.5,sp(end-1).Position(1)+pwidth.*1.5],...
     [sp(end-1).Position(2),sp(end-1).Position(2)+pheight.*1.5],...
     'Color','m','LineWidth',1);      
line(fax,...
     [sp(end-1).Position(1),sp(end-1).Position(1)],...
     [sp(end-1).Position(2),sp(end-1).Position(2)+pheight.*1.5],...
     'Color','m','LineWidth',1);      
uistack(fax,'top');



title('Spatial Rate Map');   
chax = colorbar(sp(end));
chax.LineWidth = 1;
chax.Units = 'centimeters';
colormap(sp(end),'jet');
caxis([0,maxPfsRate]);
chax.Position(1) = sp(end).Position(1)+pwidth*1.5+0.1;

sp(end).YTickLabel = {};
sp(end).XTickLabel = {};
line(fax,xpos(xind)+[0,pwidth*0.75],ypos(yind).*[1,1]-0.5,'Color','k','LineWidth',1);
text(fax,xpos(xind),ypos(yind)-0.75,'50cm','FontSize',8,'Color',[0,0,0]);

xind = xind+3;

% $$$ sp(end+1) = axes('Units','centimeters',...
% $$$                  'Position',[xpos(xind),ypos(yind),(pwidth+xpad)*2-xpad,(pheight+ypad)*2-ypad],...
% $$$                  'FontSize', 8);
% $$$ dfs{2}.plot(u(2),'mean',false,[0,maxPfsRate],false,0.5,false,[],@jet);
% $$$ text(-1.95,-1.65,num2str(cond_round(dfsMaxRatesMean(2))),'FontSize',12,'Color',[1,1,1]);
% $$$ ylabel({'Body Speed','(log10(cm/s))'});





%% MjgER2016F3B - behavior and place field examples %%


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

% {accg,DFS{HP,BP},DFS{HP,BS},theta,rear,hloc,lloc,hpause,lpause}




labels = {'HPxBP','Theta','Rear','HLoc','LLoc','HPause','LPause'};

% PLOT 
%clf();    
yind = 3;
for u = cluMap',

    xind = 1;

    sp = gobjects([1,0]);
    
% SET color scale max
    maxPfsRate = max(cell2mat(cf(@(p,u) maxRate(p,u,false,'prctile99',0.5),...
                                 [pfs,dfs],repmat({u(2)},[1,numel(pfs)+numel(dfs)]))));

    pfsMaxRatesMean = cell2mat(cf(@(p,u) max(p.maxRate(u,true,'mean')),pfs,repmat({u(2)},[1,numel(pfs)])));
    %dfsMaxRatesMean = cell2mat(cf(@(p,u) max(p.maxRate(u,false,'mean')),dfs,repmat({u(2)},[1,numel(dfs)])));
    dfsMaxRatesMean = mean(dfs{1}.maxRate(u,false,1));

% ACCG 
% $$$     sp(end+1) = axes('Units','centimeters',...
% $$$                      'Position',[xpos(xind),ypos(yind),pwidth,pheight],...
% $$$                      'FontSize', 8);
% $$$     bar(tbins,accg(:,u(2)));
% $$$     axis('tight');
% $$$     sp(end).YTickLabel = {};
% $$$     sp(end).XTickLabel = {};
% $$$     if yind == 3, title(labels{xind});end
% $$$     xind = xind+1;
    

% DRZFIELDS 
    for s = 1:numel(dfs),
        sp(end+1) = axes('Units','centimeters',...
                         'Position',[xpos(xind),ypos(yind),pwidth,pheight],...
                         'FontSize', 8);
        
        dfs{s}.plot(u(2),'mean',false,[0,maxPfsRate],false,0.5,false,[],@jet);
        sp(end).YTickLabel = {};
        sp(end).XTickLabel = {};        
        text(-2,-0.475,num2str(cond_round(dfsMaxRatesMean(s))),'FontSize',8,'Color',[1,1,1]);
        xlim([-2,nonzeros(xlim.*[0,1])-0.2])
        sp(end).Color = [0,0,0];
              
        if yind == 3, title(labels{xind});end        
        xind = xind+1;
    end

    
    for s = 1:numStates,
% PLACEFIELDS MTAApfs
        sp(end+1) = axes('Units','centimeters',...
                         'Position',[xpos(xind),ypos(yind),pwidth,pheight],...
                         'FontSize', 8);

        plot(pfs{s},u(2),'mean',false,[0,maxPfsRate],true,0.5,false,interpPar,@jet);
        sp(end).YTickLabel = {};
        sp(end).XTickLabel = {};
        if yind == 3, title(labels{xind});end        
        xind = xind+1;
        text(-495,-380,num2str(cond_round(pfsMaxRatesMean(s))),'FontSize',8,'Color',[1,1,1]);
    end
    yind = yind+1;

    
end



%% MjgER2016F3C - behavior field erpPCA eigenvectors
% see: req20180123_vis_HBPITCHxBPITCH_erpPCA.m



if ~exist('pfd','var'),
    [pfd,tags,eigVec,eigVar,eigScore,validDims,unitSubsets,unitIntersection,zrmMean,zrmStd] = req20180123_ver5(Trials);
end

numComp = size(eigVec{1},2);
pfindex = 1;


bins = dfs{1}.adata.bins;

% LOAD Behavioral state contours
[stateContourMaps,stateContourHandles] =                           ...
    bhv_contours(sessionListName,                                  ... sessionListName
                 'fet_HB_pitchB',                                  ... featureSet
                 [1,2],                                            ... featureInd
                 {linspace(-2,2,50),linspace(-2,2,50)},            ... featureBin
                 'Ed05-20140529.ont.all',                          ... referenceTrial
                 {'lloc+lpause&theta','hloc+hpause&theta',         ... states
                   'rear&theta'},                                  ...
                 'wcr'                                             ... stateColors
);




fpc  = cell([1,numComp]);
for i = 1:numComp,
    fpc{i} = nan(size(validDims{pfindex}));
    fpc{i}(validDims{pfindex}) = eigVec{pfindex}(:,i);
end

fpcMinMax = [min(cellfun(@min,fpc)),max(cellfun(@max,fpc))];

% PLOT 
yind = yind + 1;
for i = 1:3,
    sp(end+1) = axes('Units','centimeters',...
                     'Position',[xpos(i),ypos(yind)+0.5,pwidth,pheight],...
                     'FontSize', 8);    
    
    imagescnan({bins{:},abs(reshape_eigen_vector(fpc{i},pfd(1,pfindex)))},...
               fpcMinMax,'linear',false,[0,0,0],1,1);                % PRINT eigenvectors
    axis('xy');
    axis('tight');
    hold('on');    
    for s = 1:numel(stateContourHandles),                            % OVERLAY state Contours
        copyobj(stateContourHandles{s},sp(end));
    end
    sp(end).YTickLabel = {};
    sp(end).XTickLabel = {};
    sp(end).Color = [0,0,0];
        
    xlim(pfd{1}.adata.bins{1}([1,end]));
    xlim([-2,nonzeros(xlim.*[0,1])-0.2])            
    ylim(pfd{1}.adata.bins{2}([1,end]));

    %yind = yind + 2;
end



af(@(ax) set(ax,'LineWidth',1), sp);

%% MjgER2016F3D t-SNE mapping of fscores within HPxBP of first 3 eigenvectors
% see: req20180319.m
MjgER2016_load_bhv_erpPCA_scores();
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




% GET auxiliary features
clu =  cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:), pfd(:,pfindex),units');    
si  =  cf(@(p,u) p.data.si(:,ismember(p.data.clu,u),:),  pfd(:,pfindex),units');
tlu =  cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
si  = cat(2, si{:});   
clu = cat(2, clu{:});    
tlu = cat(2, tlu{:});    
clu = [tlu',clu'];


[~,rind] = sortrows(clu);
si  = si(:,rind);


% $$$ dsm = pdist([FSrC(:,1:3),si(unitSubsets{pfindex})']);
% $$$ mapa = mdscale(dsm,2);
% $$$ figure,plot(mapa(:,1),mapa(:,2),'.');


yind = yind - 2;
xind = 4;
sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind),ypos(yind),(pwidth+0.1)*6-0.1,(pheight+0.1)*6-0.1],...
                 'FontSize', 8);    

sigUnits = any(abs(fsrcz(:,1:3))>=1.96,2);

cc = eigScore{pfindex}(:,[2,1,3])+0.75;
cc(~sigUnits,:) = repmat([0.75,0.75,0.75],[sum(~sigUnits),1]);


% $$$ mapa = tsne([FSrC(:,1:3),si(unitSubsets{pfindex})'],[],2,4,225);
mapa = tsne([fsrcz(:,1:3),si(unitSubsets{pfindex})'],[],2,4,25);
% $$$ figure,scatter(mapa(:,1),mapa(:,2),5,cc,'filled');


% $$$ ss = ones([size(cc,1),1])*10;
% $$$ ss(all(abs(fsrcz(:,1:3))>1.96,2)) = 20;
% $$$ ss(all(abs(fsrcz(:,1:3))<1.96,2)) = 5;

cla();
mi = [2,1];
%scatter3(FSrC(:,1),FSrC(:,3),FSrC(:,2),15,cc,'o','filled'); 
scatter(mapa(:,mi(1)),mapa(:,mi(2)),15,cc,'o','filled'); 
sp(end).YTickLabel = {};
sp(end).XTickLabel = {};
box('on');
axis('tight')

cluSessionSubset = cluSessionMap(unitSubsets{pfindex},:);
for u = cluMap'
    uind = find(ismember(cluSessionSubset,u','rows'));
    sigUnits(uind);
end

    
% $$$ sp(end+1)=subplot(359); hold('on');scatter(mapz(:,1)/10,mapz(:,2)/10,10,cc,'o','filled'); grid('on');
% $$$ title('tsne on zscores')

%plot(mapa(:,1),mapa(:,2),'.');
%mapz = tsne([fsrcz(:,1:3),si(unitSubsets{pfindex})'],[],2,4,25);
%figure,scatter(mapz(:,1),mapz(:,2),10,cc,'filled');


yind = yind + 6;
xind = 4;
sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind),ypos(yind),(pwidth+0.1)*6-0.1,(pheight+0.1)*4-0.1],...
                 'FontSize', 8);    

hold('on');
for i = 1:3,
    [F,X] = ecdf(fsrcz(:,i));
    cdfplot(X);
end
legend({'rear','low prone','high prone'},'Location','southeast');
Lines(-1.96,[],'k');
Lines(1.96,[],'k');

Lines(-3.1,[],'k');
Lines(3.1,[],'k');


