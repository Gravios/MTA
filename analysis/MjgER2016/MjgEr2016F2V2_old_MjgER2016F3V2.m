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


% LOAD Data
MjgER2016_load_data;
%  Variables:
%      Trials
%      units
%      cluSessionMap
%      pitchReferenceTrial
%      FigDir
%      sessionListName
%      sessionList
%      states
%      numStatesg
%      interpParPfsp
%      interpParDfs
%
%  Functions:
%      reshape_eigen_vector


tind = 20;
Trial = Trials{tind};
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
dfst = {'HPITCHxBPITCH','HPITCHxBSPEED','BPITCHxBSPEED','BPITCHxHSPEED','HPITCHxRHM'};

% LOAD accgs
[accg,tbins] = autoccg(Trial);

cond_round = @(rate) max([round(rate,0),round(rate,1)].*[rate>=10,rate<10]);


pageWidth  = 29.7;
pageHeight = 42.0;

pwidth = 1.25;
pheight = 1.25;

xpad = 0.3;
ypad = 0.3;

xpos = 3.5:(pwidth+xpad):pageWidth;
ypos = fliplr(0.5:(pheight+xpad):pageHeight-3.7);

% SET figure opts
hfig = figure(666001);
hfig.Units = 'centimeters';
hfig.Position = [1, 1, pageWidth,pageHeight];
hfig.PaperPositionMode = 'auto';




% LOAD data for
sampleRate = 8;
xyz = preproc_xyz(Trial);
ufr = Trial.ufr.copy();
ufr.create(Trial,xyz,'theta-groom-sit',units{tind},0.5);
xyz.filter('ButFilter',3,2.5,'low');
xyz.resample(sampleRate);
xyz.data = sq(xyz(:,'nose',[1,2]));
ufr.resample(xyz);
drzPfs = compute_drz(Trial,units{tind},pft,[],[],[],[],xyz);
feature = fet_HB_pitchB(Trial,sampleRate);
drzBhv = compute_drz(Trial,units{tind},dfs{1},[],[],[],interpParDfs,feature);
tper = [Trial.stc{'theta-groom-sit'}];
tper.cast('TimeSeries');
tper.resample(xyz);
lper = resample(cast([Trial.stc{'loc&gper'}],'TimeSeries'),xyz.sampleRate);
pper = resample(cast([Trial.stc{'pause&gper'}] ,'TimeSeries'),xyz.sampleRate);



clf();

dfsYlim = [-pi/2.6,pi/2.6]+0.5;
dfsXlim = [-pi/2.6,pi/2.6]-0.25;
dfsTextPos = {-1.4,1.45};

pfsTextPos = {-625,450};

%% MjgER2016F3A - behavior and place field examples %%
cluMap = [20,63];
yind = 1;
xind = 1;





% SET color scale max
u = cluMap';
maxPfsRate = max(cell2mat(cf(@(p,u) maxRate(p,u,false,'prctile99',0.5),...
                             [pfs,dfs],repmat({u(2)},[1,numel(pfs)+numel(dfs)]))))-1;
pfsMaxRatesMean = cell2mat(cf(@(p,u) max(p.maxRate(u,true,'mean')),pfs,repmat({u(2)},[1,numel(pfs)])));
dfsMaxRatesMean = cell2mat(cf(@(p,u) max(p.maxRate(u,false,'mean')),dfs,repmat({u(2)},[1,numel(dfs)])));
% $$$ % ACCG 
% $$$ sp(end+1) = axes('Units','centimeters',...
% $$$                  'Position',[xpos(xind),ypos(yind)+pheight+ypad,pwidth,pheight],...
% $$$                  'FontSize', 8);
% $$$ bar(tbins,accg(:,u(2)));
% $$$ axis('tight');


% $$$ xind = 3;
% PLOT theta example    
sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind),ypos(yind),(pwidth+xpad)*2-xpad,(pheight+ypad)*2-ypad],...
                 'FontSize', 8);

plot(pfs{1},u(2),'mean',false,[0,maxPfsRate],true,0.5,false,interpParPfs,@jet,1,[1,1,1]);
text(-500,400,num2str(cond_round(pfsMaxRatesMean(1))),'FontSize',14,'Color',[0,0,0]);
xind = xind+3;
title('Place Field')

sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind),ypos(yind),(pwidth+xpad)*2-xpad,(pheight+ypad)*2-ypad],...
                 'FontSize', 8);
dfs{1}.plot(u(2),'mean',false,[0,maxPfsRate],false,0.5,false,interpParDfs,@jet,1,[1,1,1]);
text(-1.4,1.5,num2str(cond_round(dfsMaxRatesMean(1))),'FontSize',12,'Color',[0,0,0]);
xlim(dfsXlim);
ylim(dfsYlim);
xlabel('Head-Body Pitch (rad)');
ylabel('Body Pitch (radians)');
xind = xind+3;
title('Behavior Field')


sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind),ypos(yind),(pwidth+xpad)*2-xpad,(pheight+ypad)*2-ypad],...
                 'FontSize', 8);
dfs{2}.plot(u(2),'mean',false,[0,maxPfsRate],false,0.5,false,interpParDfs,@jet,1,[1,1,1]);
text(-1.4,1.7,num2str(cond_round(dfsMaxRatesMean(2))),'FontSize',12,'Color',[0,0,0]);
ylabel({'Body Speed','(log10(cm/s))'});
xlabel('Head-Body Pitch (rad)');
xlim(dfsXlim);
ylim([-1.5,2]);
title('Behavior Field')




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

labels = {'ACCG','HPxBP','HPxBS','Theta','Rear','HLoc','HPause','LLoc','LPause'};

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
    dfsMaxRatesMean = cell2mat(cf(@(p,u,ip) max(p.maxRate(u,false,'mean',[],ip)),...
                               dfs(1),repmat({u(2)},[1,numel(1)]),...
                               repmat({interpParDfs},[1,numel(1)])));

% ACCG 
    sp(end+1) = axes('Units','centimeters',...
                     'Position',[xpos(xind),ypos(yind),pwidth,pheight],...
                     'FontSize', 8);
    bar(tbins,accg(:,u(2)));
    axis('tight');
    sp(end).YTickLabel = {};
    sp(end).XTickLabel = {};
    if yind == 3, title(labels{xind});end
    xind = xind+1;
    

% DRZFIELDS 
    for s = 1,%:numel(dfs),
        sp(end+1) = axes('Units','centimeters',...
                         'Position',[xpos(xind),ypos(yind),pwidth,pheight],...
                         'FontSize', 8);
        
        dfs{s}.plot(u(2),'mean',false,[0,maxPfsRate],false,0.5,false,interpParDfs,@jet,1,[1,1,1]);
        sp(end).YTickLabel = {};
        sp(end).XTickLabel = {};        
        if s==1, tpos = dfsTextPos;tpos{1} = tpos{1}+1.45;else,tpos={dfsTextPos{1},1.55}; end
        text(tpos{:},num2str(cond_round(dfsMaxRatesMean(s))),'FontSize',8,'Color',[0,0,0]);        
        if yind == 3, title(labels{xind});end        
        xind = xind+1;
        if s==1, ylim(dfsYlim);else, ylim([-1.5,2]);end
        xlim(dfsXlim);
    end


% SPEED Distributions for locomotion vs pause
    sp(end+1) = axes('Units','centimeters',...
                     'Position',[xpos(xind),ypos(yind),pwidth,pheight],...
                     'FontSize', 8);
    xind = xind+1;    
    uind = u(2)==units{tind};
    dper = MTADepoch([],[],-0.5<drzPfs(:,uind)&drzPfs(:,uind)<0.5&-0.5<drzBhv(:,uind)&drzBhv(:,uind)<0.5,...
                     xyz.sampleRate,...
                     xyz.sync.copy(),...
                     xyz.origin,...
                     'TimeSeries','sts',[],'tdrz','d');
    indLoc = dper&tper.data&lper.data;
    ufrLoc =  ufr(indLoc,uind);
    indPas = dper&tper.data&pper.data;
    ufrPas =  ufr(logical(indPas.data),uind);
    eds = linspace(-1,2,50);
    hold('on');
    hc = histc(log10(ufrLoc),eds);
    %    hax = bar(eds,hc./sum(hc),'histc');
    hax.FaceColor = 'b';
    hax.EdgeColor = 'none';
    hax.FaceAlpha = 0.66;
    hax.EdgeAlpha = 0;
    hc = histc(log10(ufrPas),eds);
    hax = bar(eds,hc./sum(hc),'histc');
    hax.FaceColor = 'r';
    hax.EdgeColor = 'none';
    hax.FaceAlpha = 0.66;
    hax.EdgeAlpha = 0;
    xlim(eds([1,end]));
    sp(end).YTickLabel = {};
    sp(end).XTickLabel = {};        
    
    for s = 1:numStates,
% PLACEFIELDS MTAApfs
        sp(end+1) = axes('Units','centimeters',...
                         'Position',[xpos(xind),ypos(yind),pwidth,pheight],...
                         'FontSize', 8);

        plot(pfs{s},u(2),'mean',false,[0,maxPfsRate],true,0.5,false,interpParPfs,@jet,1,[1,1,1]);
        sp(end).YTickLabel = {};
        sp(end).XTickLabel = {};
        if yind == 3, title(labels{xind});end        
        xind = xind+1;
        text(pfsTextPos{:},num2str(cond_round(pfsMaxRatesMean(s))),'FontSize',8,'Color',[0,0,0]);
        set(gca,'XColor','none','YColor','none');
    end
    yind = yind+1;

end

% $$$ 
% $$$ figure();
% $$$     hold('on');
% $$$     hc = histc(log10(ufrLoc),eds);
% $$$     hax = stairs(eds,hc./sum(hc));    
% $$$     hax.EdgeColor = 'b';    
% $$$     hc = histc(log10(ufrPas),eds);
% $$$     hax = stairs(eds,hc./sum(hc));        
% $$$     hax.EdgeColor = 'r';
% $$$     xlim(eds([1,end]));

%test


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
cf(@(h) set(h,'LineWidth',1), stateContourHandles);



fpc  = cell([1,numComp]);
for i = 1:numComp,
    fpc{i} = nan(size(validDims{pfindex}));
    fpc{i}(validDims{pfindex}) = eigVec{pfindex}(:,i);
end

fpcMinMax = [min(cellfun(@min,fpc)),max(cellfun(@max,fpc))];

% PLOT 
yind = 1;
xind = 10;
for i = 1:3,
    sp(end+1) = axes('Units','centimeters',...
                     'Position',[xpos(xind)+1,ypos(yind),pwidth*2+xpad,pheight*2+ypad],...
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
    xlim(dfsXlim);
    ylim(dfsYlim);
    xind = xind + 2;
end





%% MjgER2016F3D t-SNE mapping of fscores within HPxBP of first 3 eigenvectors
% see: req20180319.m

if ~exist('pfdShuffled','var'),
    pfdShuffled =  cf(@(t,g) MTAApfs(t,'tag',g), Trials, repmat({'HBPITCHxBPITCH_shuffled'},size(Trials)));
end

% GET bhv rate maps
rmaps = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'), pfd(:,pfindex),units');
clu =  cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:), pfd(:,pfindex),units');    
tlu =  cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
rmaps = cat(2, rmaps{:});   
clu = cat(2, clu{:});    
tlu = cat(2, tlu{:});    
clu = [tlu',clu'];
[~,rind] = sortrows(clu);
rmaps = rmaps(:,rind);
rmaps = rmaps(validDims{pfindex},unitSubsets{pfindex});
rmaps(isnan(rmaps)) = 0;

rmapsShuffledMean = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'), pfdShuffled',units');
clu =  cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:), pfdShuffled',units');    
tlu =  cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
rmapsShuffledMean = cat(2, rmapsShuffledMean{:});   
clu = cat(2, clu{:});    
tlu = cat(2, tlu{:});    
clu = [tlu',clu'];
[~,rind] = sortrows(clu);
rmapsShuffledMean = rmapsShuffledMean(:,rind);
rmapsShuffledMean = rmapsShuffledMean(validDims{pfindex},unitSubsets{pfindex});
rmapsShuffledMean(isnan(rmapsShuffledMean)) = 0;

rmapsShuffled = cf(@(p,u) p.data.rateMap(:,ismember(p.data.clu,u),:), pfdShuffled',units');
clu =  cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:), pfdShuffled',units');    
tlu =  cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
rmapsShuffled = cat(2, rmapsShuffled{:});   
clu = cat(2, clu{:});    
tlu = cat(2, tlu{:});    
clu = [tlu',clu'];
[~,rind] = sortrows(clu);
rmapsShuffled = rmapsShuffled(:,rind,:);
rmapsShuffled = rmapsShuffled(validDims{pfindex},unitSubsets{pfindex},:);
rmapsShuffled(isnan(rmapsShuffled)) = 0;

D = cov(rmapsShuffledMean');
LR = eigVec{pfindex};
% compute rotated FS coefficients        
FSCFr = LR * inv(LR' * LR);          % this is pseudo-inverse of LR
% rescale rotated FS coefficients by the corresponding SDs 
rk = size(FSCFr,2);%rank(D,1e-4);       % why not on D? would save on corrcoef(X) computation
FSCFr = FSCFr .* repmat(sqrt(diag(D)),1,rk);
% compute rotated factor scores from the normalized raw data and  the
% corresponding rescaled factor score coefficients
rsMean = mean(rmapsShuffledMean');
rsStd  = std( rmapsShuffledMean');

% MEAN shuffled score
FSrM =  ((rmapsShuffledMean'-rsMean)./rsStd) * FSCFr;

% MEAN normal score
FSrC =  ((rmaps'-rsMean)./rsStd) * FSCFr;

FSrS = [];
for i = 1:pfd{1}.parameters.numIter
    FSrS(:,:,end+1) =  ((rmapsShuffled(:,:,i)'-rsMean)./rsStd) * FSCFr;
end
fsrsMean = mean(FSrS,3);
fsrsStd = std(FSrS,[],3);

fsrcz = (FSrC-fsrsMean)./fsrsStd;



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

mapa = tsne([FSrC(:,1:3),si(unitSubsets{pfindex})'],[],2,3,23);
%figure,plot(mapa(:,1),mapa(:,2),'.');
figure,scatter(mapa(:,mi(1))/10,mapa(:,mi(2))/10,15,cc,'o','filled'); 

D = pdist(FSrC(:,1:3));
mapa = mdscale(D,2);
figure




yind = 8;
xind = 10;

sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind)+1,ypos(yind),(pwidth+xpad)*6-xpad,(pheight+ypad)*6-ypad],...
                 'FontSize', 8);    

sigUnits = any(abs(fsrcz(:,1:3))>=1.96,2);

%cc = eigScore{pfindex}(:,1:3)+0.75;
cc = FSrC(:,[1,2,3])+0.75;
%cc = FSrC(:,[2,1,3])+0.75;
cc(~sigUnits,:) = repmat([0.75,0.75,0.75],[sum(~sigUnits),1]);

% $$$ ss = ones([size(cc,1),1])*10;
% $$$ ss(all(abs(fsrcz(:,1:3))>1.96,2)) = 20;
% $$$ ss(all(abs(fsrcz(:,1:3))<1.96,2)) = 5;


cla();
mi = [2,1];
scatter(-mapa(:,mi(1)),mapa(:,mi(2)),10,cc,'filled');
%scatter(-mapa(:,mi(1))/10,-mapa(:,mi(2))/10,15,cc,'o','filled'); 
xlim([-4,4]);
sp(end).YTickLabel = {};
sp(end).XTickLabel = {};
box('on');
hold('on');
set(gca,'XColor',[0.75,0.75,0.75],'YColor',[0.75,0.75,0.75]);



cluMap = [20,74;...
          20,79;...
          20,104];


cluSessionSubset = cluSessionMap(unitSubsets{pfindex},:);
for u = cluMap'
    uind = find(ismember(cluSessionSubset,u','rows'));
    %plot(-mapa(uind,mi(1))/10,-mapa(uind,mi(2))/10,'+k');
    plot(-mapa(uind,mi(1)),mapa(uind,mi(2)),'+k');    
end

%delete(sp(end));sp(end) = [];    

% $$$ sp(end+1)=subplot(359); hold('on');scatter(mapz(:,1)/10,mapz(:,2)/10,10,cc,'o','filled'); grid('on');
% $$$ title('tsne on zscores')

%plot(mapa(:,1),mapa(:,2),'.');
%mapz = tsne([fsrcz(:,1:3),si(unitSubsets{pfindex})'],[],2,4,25);
%figure,plot(mapz(:,1),mapz(:,2),'.');

%delete(sp(end));sp(end) = [];    

yind = 11;
xind = 10;
sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind)+1,ypos(yind),(pwidth+xpad)*6-xpad,(pheight+ypad)*3-ypad-1],...
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





% $$$ 
% $$$ yind = 14;
% $$$ xind = 10;
% $$$ sp(end+1) = axes('Units','centimeters',...
% $$$                  'Position',[xpos(xind)+1,ypos(yind),(pwidth+xpad)*3-xpad,(pheight+ypad)*3-ypad-1],...
% $$$                  'FontSize', 8);    
% $$$ hold('on');
% $$$ plot(fsrcz(:,2),fsrcz(:,1),'.');
% $$$ circle(0,0,1.96,'-r');
% $$$ daspect([1,1,1])
% $$$ xlim([-8,8]);
% $$$ ylim([-8,8]);
% $$$ 
% $$$ xind = 13;
% $$$ sp(end+1) = axes('Units','centimeters',...
% $$$                  'Position',[xpos(xind)+1,ypos(yind),(pwidth+xpad)*3-xpad,(pheight+ypad)*3-ypad-1],...
% $$$                  'FontSize', 8);    
% $$$ hold('on');
% $$$ plot(fsrcz(:,3),fsrcz(:,1),'.');
% $$$ circle(0,0,1.96,'-r');
% $$$ daspect([1,1,1])
% $$$ xlim([-8,8]);
% $$$ ylim([-8,8]);

D = pdist(FSrC(:,1:3));
Y = cmdscale(D,3);
Y = mdscale(D,3);
Y = FSrC(:,1:3);
figure,
scatter3(Y(:,1),Y(:,2),Y(:,3),10,cc,'filled');
Y = mdscale(D,2);
figure
scatter(-Y(:,2),Y(:,1),10,cc,'filled');

Y = mdscale(D,2);
figure
scatter(-Y(:,2),Y(:,1),10,cc,'filled');
figure,
scatter(sign(-Y(:,2)).*(abs(-Y(:,2))).^(3/4),sign(Y(:,1)).*(abs(Y(:,1))).^(2/3),10,cc,'filled');



cf(@(t) t.load('nq'), Trials);

bRat = cf(@(t,u) t.nq.SpkWidthR(u), Trials,units);
bRat = cf(@(t,u) t.nq.bRat(u), Trials,units);
bRat = cat(1,bRat{:});
bRat = bRat(unitSubsets{1});



figure
scatter(mapa(:,2),mapa(:,1),20,log10(bRat),'filled');
colormap(jet)



figure,plot(bRat,FSrC(:,1),'.')
figure,plot(bRat,FSrC(:,2),'.')
figure,plot(bRat,FSrC(:,3),'.')

figure
scatter(-Y(:,2),Y(:,1),20,log10(bRat+1),'filled');
colormap(jet)

figure,
scatter(fsrcz(:,2),fsrcz(:,1),20,log10(bRat+1),'filled');
colormap('jet')


