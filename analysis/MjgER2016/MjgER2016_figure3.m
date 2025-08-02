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

%%%<<< LOAD MjgER2016 data
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
%%%>>>

%%%<<< SET MjgER2016 figure 3 global default arguments
MjgER2016_figure3_args();
% MjgER2016_figure3_args:
%
% Default arument overrides:
%     fet_HB_pitchB
%     compute_bhv_ratemaps
%     compute_bhv_ratemaps_shuffled
%     compute_bhv_ratemap_erpPCA
%%%>>>


%%%<<< LOAD data
        % EXAMPLE Trial : jg05-20120312
tind = 20;
Trial = Trials{tind};
stc = Trials{tind}.stc.copy();


% LOAD place restricted behavior fields


bhvMask = false(size(validDims));
bhvMask(validDims) = true;
bhvMask = reshape_eigen_vector(bhvMask)';

bhvLims = [-1.6, 0.6; ...
           -0.5, 1.7];


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


% COMPUTE place field centers for each state
pfssMr = {};
pfssMp = {};
uCount = 0;
for t = 1:numel(Trials),
    [mr,mp] = cf(@(p,u) p.maxRate(units{t},'interpPar',interpParPfs), pfsa{t});
    for s = 1:numel(mr),
        pfssMr{s}(1+uCount:uCount+numel(mr{1}),:) = mr{s}(:,1);
        pfssMp{s}(1+uCount:uCount+size(mp{1},1),:) = mp{s}(:,:);        
    end
    uCount = uCount+numel(mr{1});
end
pfssMrSub = cat(3,pfssMr{:});
pfssMrSub = pfssMrSub(unitSubset,:,:);
pfssMpSub = cat(3,pfssMp{:});
pfssMpSub = pfssMpSub(unitSubset,:,:);


%%%>>>


%%%<<< LOAD Behavioral state contours
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
%%%>>>


%%%<<< LOAD and TRANSFORM(tSNE) erpPCA scores

bfs        = cf(@(t,u)  compute_bhv_ratemaps(t,u),                               Trials, units);
bfsShuff   = cf(@(t,u)  compute_bhv_ratemaps_shuffled(t,u),                      Trials, units);

[fsrcz,FSrC,rmaps,FSCFr,FSrM,FSrS,fsrsMean,fsrsStd,rmapsShuffledMean,rmapsShuffled] = ...
    compute_bhv_ratemaps_erpPCA_scores(Trials,units,bfs,bfsShuff,eigVecs,validDims,unitSubset,true);    
% $$$     compute_bhv_ratemaps_erpPCA_scores('tag','68536b9a930748fd25a91dfc8ca84d39');


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

cluSessionSubset = cluSessionMap(unitSubset,:);
for u = cluMap'
    uind = find(ismember(cluSessionSubset,u','rows'));
    sigUnits(uind);
end

cc = eigScrs(:,[2,1,3])+0.75;
cc(~sigUnits,:) = repmat([0.75,0.75,0.75],[sum(~sigUnits),1]);


if ~exist('mapa','var'),
mapa = tsne([FSrC(:,1:3),(si(unitSubset)'-mean(si(unitSubset))')./std(si(unitSubset))'],[],2,3,40);    
% $$$     mapa = tsne([fsrcz(:,1:3),(si(unitSubset)'-mean(si(unitSubset))')./std(si(unitSubset))'],[],2,4,80);
% $$$     mapa = tsne([fsrcz(:,1:3),si(unitSubset)'],[],2,3,50);
end

%%%>>>


%%%<<< SET Figure Opts
cond_round = @(rate) max([round(rate,0),round(rate,1)].*[rate>=10,rate<10]);
nanColor = [0.25,0.25,0.25];

pageWidth  = 21.;
pageHeight = 29.7;
pwidth = 1.5;
pheight = 1.5;
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

fax = axes('Position',[0,0,1,1],'Visible','off','Units','centimeters');
xlim([0,hfig.Position(3)]);
ylim([0,hfig.Position(4)]);

sp = gobjects([1,0]);

%%%>>>


%%%<<< MjgER2016-F-3A: behavior and place field examples

cluMap = [20,119];
yind = 1;
xind = 1;


u = cluMap';
% SET color scale max
maxPfsRate = max(cell2mat(cf(@(p,u) maxRate(p,u,false,'prctile99',0.5),...
                             {pfs{1},bfsEx},repmat({u(2)},[1,2]))));


xind = 1;
sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind),ypos(yind),pwidth*1.25,pheight*1.25],...
                 'FontSize', 8);
bfs{tind}.plot(u(2),'mean',false,[0,maxPfsRate],false,0.5,false,[],@jet,[],nanColor);
text(-2,-0.45,num2str(cond_round(bfs{tind}.maxRate(u(2),false,1))),'FontSize',10,'Color',[1,1,1]);
xlabel({'HP (rad)'});
ylabel({'BP (rad)'});
xind = xind+2;
xlim([-2,nonzeros(xlim.*[0,1])-0.2])
sp(end).Color = nanColor;
title({'Behavior','Rate Map'});    


sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind),ypos(yind),pwidth*1.25,pheight*1.25],...
                 'FontSize', 8);
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

%%%>>>


%%%<<< MjgER2016-F-3B: behavior and place field examples

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
          20,85;...
          20,83;...
          20,141;...
          20,79;...
          20,103;...
          20,104];
% ...          20,109;...
%          20,59];

% {accg,BFS{HP,BP},BFS{HP,BS},theta,rear,hloc,lloc,hpause,lpause}


labels = {{'Behavior','Rate Map'},{'Theta','Rate Map'},'Rear','H Loc','H Pause','L Loc','L Pause'};

% PLOT 
xind = 1;
yind = 1;
yinit = yind;
for u = cluMap',

    xind = 1;

% SET color scale max
    maxPfsRate = max(cell2mat(cf(@(p,u,m) ...
                                 p.maxRate(u,true&numel(m)==1,1,[],[],[],m),...
                                 [pfs,{bfsEx}],repmat({u(2)},[1,numel(pfs)+1]),[repmat({true},[1,numel(pfs)]),{bhvMask(:)}])));

    pfsMaxRatesMean = cell2mat(cf(@(p,u) ...
                                  max(p.maxRate(u,true,'mean')),...
                                  pfs,repmat({u(2)},[1,numel(pfs)])));
    %bfsMaxRatesMean = cell2mat(cf(@(p,u) max(p.maxRate(u,false,'mean')),bfs,repmat({u(2)},[1,numel(bfs)])));
    bfsMaxRatesMean = max(bfsEx.maxRate(u(2),false,1,[],[],[],bhvMask(:)));

    

% DRZFIELDS 
    sp(end+1) = axes('Units','centimeters',...
                     'Position',[xpos(xind),ypos(yind),pwidth,pheight],...
                     'FontSize', 8);

    bfsEx.plot(u(2),'mean',false,[0,maxPfsRate],false,0.5,false,[],@jet,bhvMask,nanColor);
    sp(end).YTickLabel = {};
    sp(end).XTickLabel = {};        
    %text(-0.1,1.5,num2str(cond_round(bfsMaxRatesMean)),'FontSize',8,'Color',[1,1,1]);
    text(0.55,1.5,num2str(round(bfsMaxRatesMean)),'FontSize',8,'Color',[1,1,1],'HorizontalAlignment','right');
    xlim(bhvLims(1,:));
    ylim(bhvLims(2,:));
    sp(end).Color = nanColor;

    if yind == yinit, 
        title(labels{xind});
    end 
    if yind==round(size(cluMap,1)./2),
        ylabel({'Head Pitch'});
    end
    if yind==size(cluMap,1),
        xlabel({'Body Pitch'});
    end
    
    xind = 2;
    for s = 1:numStates,
% PLACEFIELDS MTAApfs
        sp(end+1) = axes('Units','centimeters',...
                         'Position',[xpos(xind)+0.2+0.2.*double(s>=2),ypos(yind),pwidth,pheight],...
                         'FontSize', 8);

        plot(pfs{s},u(2),'mean',false,[0,maxPfsRate],true,0.5,false,interpParPfs,@jet,[],nanColor);
        sp(end).YTickLabel = {};
        sp(end).XTickLabel = {};
        if yind == yinit, 
            tax = title(labels{xind});
            if s==1,
                %tax.Position(2) = 800;
            end        
        end
        
        xind = xind+1;
        text(490,380,num2str(round(pfsMaxRatesMean(s))),'FontSize',8,'Color',[1,1,1],'HorizontalAlignment','right');
        %text(200,380,num2str(cond_round(pfsMaxRatesMean(s))),'FontSize',8,'Color',[1,1,1]);
        
        if xind ==3
            hold('on');
            rmap = plot(pfs{1},u(2),'mean',false,[0,maxPfsRate],true,0.5,false,interpParPfs);
            [~,mxp] = pft.maxRate(u(2));
            binGrid = cell([1,2]);
            [binGrid{:}] = ndgrid(interpParPfs.bins{:});
            binGrid = cat(3,binGrid{:});
            binDist = sqrt(sum((binGrid-repmat(permute(mxp,[1,3,2]),[cellfun(@numel,interpParPfs.bins),1])).^2,3));
            [cont,conHax] = contour(binGrid(:,:,1),binGrid(:,:,2),(rmap./max(rmap(:)))>0.2&binDist<250,[0.5,0.5],'m','LineWidth',1);

            % REMOVE first coordinates for this example
            cont(:,1) = [];
            
            [contYMax,contYMaxInd] = max(reshape(binGrid(:,:,2).*double((rmap./max(rmap(:)))>0.2&binDist<250),[],1));
            [contYMin,contYMinInd] = min(reshape(binGrid(:,:,2).*double((rmap./max(rmap(:)))>0.2&binDist<250),[],1));
            contMaxInd = zeros([1,2]);
            [contMaxInd(1),contMaxInd(2)] = ind2sub(cellfun(@numel,interpParPfs.bins),contYMaxInd);
            contMinInd = zeros([1,2]);
            [contMinInd(1),contMinInd(2)] = ind2sub(cellfun(@numel,interpParPfs.bins),contYMinInd);
            contMaxPos = sq(binGrid(contMaxInd(1),contMaxInd(2),:));
            contMinPos = sq(binGrid(contMinInd(1),contMinInd(2),:));
            
            [contMaxYord, contMaxInd] = max(cont(2,:));
            contMaxXord = cont(1,contMaxInd);
            [contMinYord, contMinInd] = min(cont(2,:));
            contMinXord = cont(1,contMinInd);
            
            [~,contMinXind] = NearestNeighbour(interpParPfs.bins{1},contMinXord);
            [~,contMinYind] = NearestNeighbour(interpParPfs.bins{2},contMinYord);            
            [~,contMaxXind] = NearestNeighbour(interpParPfs.bins{1},contMaxXord);
            [~,contMaxYind] = NearestNeighbour(interpParPfs.bins{2},contMaxYord);            

            % TOP line
            line(fax,...
                 [sp(end-1).Position(1)+pwidth, sp(end).Position(1)+(contMaxXind./numel(interpParPfs.bins{1})).*pwidth],...
                 [sp(end-1).Position(2)+pheight,sp(end).Position(2)+(contMaxYind./numel(interpParPfs.bins{2})).*pheight],...
                 'Color','m','LineWidth',1);
            % BOTTOM line            
            line(fax,...
                 [sp(end-1).Position(1)+pwidth, sp(end).Position(1)+(contMinXind./numel(interpParPfs.bins{1})).*pwidth],...
                 [sp(end-1).Position(2),sp(end).Position(2)+(contMinYind./numel(interpParPfs.bins{2})).*pheight],...
                 'Color','m','LineWidth',1);
% $$$             % TOP line
% $$$             line(fax,...
% $$$                  [sp(end-1).Position(1)+pwidth,sp(end).Position(1)+(contMaxInd(1)./numel(interpParPfs.bins{1})).*pwidth],...
% $$$                  [sp(end-1).Position(2)+pheight,sp(end).Position(2)+(contMaxInd(2)./numel(interpParPfs.bins{2})).*pheight],...
% $$$                  'Color','m','LineWidth',1);
% $$$             % BOTTOM line            
% $$$             line(fax,...
% $$$                  [sp(end-1).Position(1)+pwidth,sp(end).Position(1)+(contMinInd(1)./numel(interpParPfs.bins{1})).*pwidth],...
% $$$                  [sp(end-1).Position(2),sp(end).Position(2) + (contMinInd(2)./numel(interpParPfs.bins{2})).*pheight], ...
% $$$                  'Color','m','LineWidth',1);
            
            rectangle(fax,'Position',sp(end-1).Position+[-0.02,-0.02,0.04,0.04],'EdgeColor','m');
% $$$             % 
% $$$             line(fax,                                                                             ...
% $$$                  [sp(end-1).Position(1),sp(end-1).Position(1)+pwidth],                      ...
% $$$                  [sp(end-1).Position(2),sp(end-1).Position(2)],                                   ...
% $$$                  'Color','m','LineWidth',1);      
% $$$             line(fax,                                                                             ...
% $$$                  [sp(end-1).Position(1),sp(end-1).Position(1)+pwidth],                      ...
% $$$                  [sp(end-1).Position(2)+pheight,sp(end-1).Position(2)+pheight],       ...
% $$$                  'Color','m','LineWidth',1);      
% $$$             line(fax,                                                                             ...
% $$$                  [sp(end-1).Position(1)+pwidth,sp(end-1).Position(1)+pwidth],         ...
% $$$                  [sp(end-1).Position(2),sp(end-1).Position(2)+pheight],                     ...
% $$$                  'Color','m','LineWidth',1);      
% $$$             line(fax,                                                                             ...
% $$$                  [sp(end-1).Position(1),sp(end-1).Position(1)],                                   ...
% $$$                  [sp(end-1).Position(2),sp(end-1).Position(2)+pheight],                     ...
% $$$                  'Color','m','LineWidth',1);      
            uistack(fax,'top');
            end
    end
    yind = yind+1;

    
end



% DRAW line and arrows for theta to behavioral state decomposition
s = 7*size(cluMap,1)-2;
%sp(end-s:end-s+5)
% HORIZONTAL line
line(fax,...
     [sp(end-s).Position(1)+pwidth-0.2,sp(end-(s-5)).Position(1)+pwidth/2+0.07],...
      repmat(sp(end-s).Position(2)+sp(end-s).Position(4)+0.575,[1,2]),...
     'Color','k',...
     'LineWidth',1);
for a = 1:5,
    patch(fax,                                                                                    ...
          repmat(sp(end-(s-a)).Position(1)+pwidth/2,[1,3])+[-0.07,0.07,0],                          ...
          repmat(sp(end-s).Position(2)+sp(end-s).Position(4)+0.475,[1,3])+[0.1,0.1,-0.09],               ...
          'k');
end


%%%>>> 

%%%<<< MjgER2016-F-3B-2: plot placefield centers
csteps = floor(2*pi*480);
grid   = [1:csteps];
angle  = grid*2*pi/csteps;
xind = 2;
for s = 1:6,
        sp(end+1) = axes('Units','centimeters',...
                         'Position',[xpos(xind)+0.2+0.2.*double(s>=2),ypos(yind),pwidth,pheight],...
                         'FontSize', 8,...
                         'Color', nanColor);
        hold('on');
        % DRAW disk
        pax = patch( sin(angle) * 480, cos(angle) * 480, [0.8,0.8,0.8]);        
        sind = pfssMrSub(:,1,s)>3;
        plot(pfssMpSub(sind,1,s)+randn([sum(sind),1])*10,...
             pfssMpSub(sind,2,s)+randn([sum(sind),1])*10,'.b');
        xind = xind +1;
        sp(end).YTickLabel = {};
        sp(end).XTickLabel = {};        

end


%%%>>>


%%%<<< MjgER2016-F-3C: behavior field erpPCA eigenvectors

% PLOT 
yind = size(cluMap,1)+3;
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
    if i == 1,
        title({'erpPCA','Factors'});
    end
end

af(@(ax) set(ax,'LineWidth',1), sp);

%%%>>>


%%%<<< MjgER2016-F-3D: t-SNE mapping of fscores within HPxBP of first 3 eigenvectors

% DATA SOURCE:  compute_bhv_ratemaps_erpPCA_scores

yind = yind - 1;
xind = 2;
sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind)+0.2,ypos(yind)+0.5,(pwidth+0.1)*3-0.1,(pheight+0.1)*3-0.1],...
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
scatter(mapa(:,mi(1)),mapa(:,mi(2)),10,cc,'o','filled'); 
sp(end).YTickLabel = {};
sp(end).XTickLabel = {};
box('on');
axis('tight')
xlim(xlim.*1.05);
ylim(ylim.*1.05);

set(gca(),'Color',nanColor)
title({'t-SNE Transformed','erpPCA F-Scores'});
    
% $$$ sp(end+1)=subplot(359); hold('on');scatter(mapz(:,1)/10,mapz(:,2)/10,10,cc,'o','filled'); grid('on');
% $$$ title('tsne on zscores')

%plot(mapa(:,1),mapa(:,2),'.');
%mapz = tsne([fsrcz(:,1:3),si(unitSubsets{pfindex})'],[],2,4,25);
%figure,scatter(mapz(:,1),mapz(:,2),10,cc,'filled');

%%%>>>


%%%<<< MjgER2016-F-3E: ECDF of eigenvector significances
yind = yind-1;
xind = 5;
sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind)+0.4,ypos(yind)+0.5+pheight/2,(pwidth+0.1)*3-0.1,(pheight+0.1)*1.5],...
                 'FontSize', 8);    

hold('on');
sclr = 'grb';
for i = [2,1,3],
    [F,X] = ecdf(fsrcz(:,i));
    cax = cdfplot(X);
    cax.Color = sclr(i);
end
set(sp(end),'YAxisLocation','right')
lax = legend({'rear','high','low'},'Location','southoutside');
lax.Units = 'centimeters';
lax.Position(1) = sp(end).Position(1);
drawnow();
lax.Position(2) = sp(end).Position([2])-1.75;
Lines(-1.96,[],'k');
Lines(1.96,[],'k');

Lines(-3.1,[],'k');
Lines(3.1,[],'k');
xlim([-20,20])

%%%>>>

