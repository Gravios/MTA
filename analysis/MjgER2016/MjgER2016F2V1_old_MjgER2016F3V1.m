
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


sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);

Trials = af(@(s) MTATrial.validate(s), sessionList);
units = cf(@(T)  select_placefields(T),  Trials); 
units = req20180123_remove_bad_units(units);

cluSessionMap = [];
for u = 1:numel(units)
    cluSessionMap = cat(1,cluSessionMap,[u*ones([numel(units{u}),1]),units{u}(:)]);
end

pitchReferenceTrial = 'Ed05-20140529.ont.all';


% SET helper function to reshape eigenvectors
reshape_eigen_vector = @(V,pfd) reshape(V(:,1),pfd{1}.adata.binSizes')';
 

FigDir = create_directory('/storage/gravio/figures/analysis/parts/MjgER2016/');

% LOAD session list
sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);


% SET states to plot
states = {'theta-groom-sit','rear&theta','hloc&theta','lloc&theta',     ...
          'hpause&theta','lpause&theta'};
numStates = numel(states);


interpPar = struct('bins',{{linspace(-500,500,100),linspace(-500,500,100)}},...
                   'nanMaskThreshold', 0,...
                   'methodNanMap',     'linear',...
                   'methodRateMap',    'linear');



Trial = Trials{20};
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


cluMap = [20,63];


yind = 1;
xind = 1;

% SET color scale max
u = cluMap';
maxPfsRate = max(cell2mat(cf(@(p,u) maxRate(p,u,false,'prctile99',0.5),...
                             [pfs,dfs],repmat({u(2)},[1,numel(pfs)+numel(dfs)]))));
pfsMaxRatesMean = cell2mat(cf(@(p,u) max(p.maxRate(u,true,'mean')),pfs,repmat({u(2)},[1,numel(pfs)])));
dfsMaxRatesMean = cell2mat(cf(@(p,u) max(p.maxRate(u,false,'mean')),dfs,repmat({u(2)},[1,numel(dfs)])));



% ACCG 
sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind),ypos(yind)+pheight+ypad,pwidth,pheight],...
                 'FontSize', 8);
bar(tbins,accg(:,u(2)));
axis('tight');


xind = 3;

% PLOT theta example    
sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind),ypos(yind),(pwidth+xpad)*2-xpad,(pheight+ypad)*2-ypad],...
                 'FontSize', 8);

plot(pfs{1},u(2),'mean',false,[0,maxPfsRate],true,0.5,false,interpPar,@jet);
text(-490,-420,num2str(cond_round(pfsMaxRatesMean(1))),'FontSize',14,'Color',[1,1,1]);
xind = xind+3;

sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind),ypos(yind),(pwidth+xpad)*2-xpad,(pheight+ypad)*2-ypad],...
                 'FontSize', 8);
dfs{1}.plot(u(2),'mean',false,[0,maxPfsRate],false,0.5,false,[],@jet);
text(-1.95,-1.65,num2str(cond_round(dfsMaxRatesMean(1))),'FontSize',12,'Color',[1,1,1]);
xlabel('Head-Body Pitch (rad)');
ylabel('Body Pitch (radians)');
xind = xind+3;

sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind),ypos(yind),(pwidth+xpad)*2-xpad,(pheight+ypad)*2-ypad],...
                 'FontSize', 8);
dfs{2}.plot(u(2),'mean',false,[0,maxPfsRate],false,0.5,false,[],@jet);
text(-1.95,-1.65,num2str(cond_round(dfsMaxRatesMean(2))),'FontSize',12,'Color',[1,1,1]);
ylabel({'Body Speed','(log10(cm/s))'});





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




labels = {'ACCG','HPxBP','HPxBS','Theta','Rear','HLoc','LLoc','HPause','LPause'};

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
    dfsMaxRatesMean = cell2mat(cf(@(p,u) max(p.maxRate(u,false,'mean')),dfs,repmat({u(2)},[1,numel(dfs)])));
    

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
    for s = 1:numel(dfs),
        sp(end+1) = axes('Units','centimeters',...
                         'Position',[xpos(xind),ypos(yind),pwidth,pheight],...
                         'FontSize', 8);
        
        dfs{s}.plot(u(2),'mean',false,[0,maxPfsRate],false,0.5,false,[],@jet);
        sp(end).YTickLabel = {};
        sp(end).XTickLabel = {};        
        %title(sprintf('Max Rate: %3.2f',dfsMaxRatesMean(s)));    
        text(-2,-1.5,num2str(cond_round(dfsMaxRatesMean(s))),'FontSize',8,'Color',[1,1,1]);
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
        text(-495,-350,num2str(cond_round(pfsMaxRatesMean(s))),'FontSize',8,'Color',[1,1,1]);
    %title(sprintf('Max Rate: %3.2f',pfsMaxRatesMean(s)));
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
yind = yind + 2;
for i = 1:3,
    sp(end+1) = axes('Units','centimeters',...
                     'Position',[xpos(1),ypos(yind),pwidth*2+0.1,pheight*2+0.1],...
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
    
    xlim([-2,2]);
    ylim([-2,2]);
    yind = yind + 2;
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

% $$$ mapa = tsne([FSrC(:,1:3),si(unitSubsets{pfindex})'],[],2,3,22);
% $$$ figure,plot(mapa(:,1),mapa(:,2),'.');


yind = yind - 2;
xind = 4;
sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind),ypos(yind),(pwidth+0.1)*6-0.1,(pheight+0.1)*6-0.1],...
                 'FontSize', 8);    

sigUnits = any(abs(fsrcz(:,1:3))>=1.96,2);

cc = eigScore{pfindex}(:,1:3)+0.75;
cc(~sigUnits,:) = repmat([0.75,0.75,0.75],[sum(~sigUnits),1]);

% $$$ ss = ones([size(cc,1),1])*10;
% $$$ ss(all(abs(fsrcz(:,1:3))>1.96,2)) = 20;
% $$$ ss(all(abs(fsrcz(:,1:3))<1.96,2)) = 5;

cla();
mi = [1,2];
scatter(mapa(:,mi(1))/10,mapa(:,mi(2))/10,15,cc,'o','filled'); 
xlim([-4.5,4.5]);
sp(end).YTickLabel = {};
sp(end).XTickLabel = {};
box('on');

cluSessionSubset = cluSessionMap(unitSubsets{pfindex},:);
for u = cluMap'
    uind = find(ismember(cluSessionSubset,u','rows'));
    sigUnits(uind)
end

    
% $$$ sp(end+1)=subplot(359); hold('on');scatter(mapz(:,1)/10,mapz(:,2)/10,10,cc,'o','filled'); grid('on');
% $$$ title('tsne on zscores')

%plot(mapa(:,1),mapa(:,2),'.');
%mapz = tsne([fsrcz(:,1:3),si(unitSubsets{pfindex})'],[],2,4,25);
%figure,plot(mapz(:,1),mapz(:,2),'.');


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


