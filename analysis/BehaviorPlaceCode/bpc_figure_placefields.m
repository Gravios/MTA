% BehaviorPlaceCode BhvPlacefields
%
% HP:= Head Pitch
% BP:= Body Pitch
% PFD:= behavior field restricted to DRZ[-0.5,0.5]
%
% Goals:
% 1. Show the behavioral dependence of place cell firing rates
%  
%  Given:
%
%   -  The firing rate of place cells change depending on the environment and the cues placed within and around
%      the environment. 
%   -  T
%
%  Hypothesis:
%      The firing rate of place cells change dependent upon the active behavioral state.
% 
%  Measure:
%      Mean Firing rate within field 
%        1. definition of field - Contiguous patch whos boarders are half the maximum firing rate
%
%  Results:
%   -  Mean peak firing rate of cells within the primary theta patch
%   -  
%   -  The firing rate of several place cells change significantly depending upon behavioral state
%   ?  Definitions of signifcance.
%      -  within patch state permutations of mean rate
%         1. equal sample size 
%         2. 3 state pairs 
%            a. rear vs high
%            b. rear vs low
%            c. high vs low
%         3. 
% 
% 2. Intra field firing maps onto the
% 
% Subplots:
%    A. Place field examples
%        1. ratemaps maze theta
%        2. ratemaps maze rear & theta
%        3. ratemaps maze hloc & theta
%        4. ratemaps maze lloc & theta
%        5. ratemaps maze hpause & theta
%        6. ratemaps maze lpause & theta
%        7. ratemaps posture HPxBP (first patch)
%        8. ratemaps posture HPxBP (second patch)
%
%    --- Examination of the neurons' behavior-space-information distribution.
%    ??? How many neurons have significant behavioral information ???
%    ??? should the erpPCA only include those with siginificant behavioral information ???
%
%    B. Eigenvectors of PFD decomposition of HPxBP space
%    C. t-SNE mapping of fscores within HPxBP of first 3 eigenvectors
%    D. Empirical CDF of zscores
%
% Suplementary plots
%    Place field centers for each state where the max rate is greater than 3 Hz.
%    For each cell get max rate in each state
%    Plot position of place field
%        rear vs low
%        rear vs high
%        high vs low
%    
%    
%    
% TODO 
%       Better differentiate behavior space rate maps from maze space rate maps
%  


% >>> LOAD BehaviorPlaceCode data >>> -----------------------------------------

%  BehaviorPlaceCode_load_data:
%
%  Variables:
%      Trials      Units    cluSessionMap   pitchReferenceTrial  FigDir
%      sessionListName      sessionList     states               numStates      
%      interpParPfsp        interpParDfs
%      
%  Functions:
%      reshape_eigen_vector
bpc_load_data();

% <<< LOAD BehaviorPlaceCode data <<< -----------------------------------------

% >>> GLOBAL ARGS function input override >>> ---------------------------------
% SET global overide pramaters for the listed functions.
% - fet_HB_pitchB
% - compute_bhv_ratemaps
% - compute_bhv_ratemaps_shuffled
% - compute_ratemaps
% - compute_bhv_ratemap_erpPCA
configure_default_args();
% <<< GLOBAL ARGS function input override <<< ---------------------------------

% >>> EXAMPLE Trial : jg05-20120312 >>> ---------------------------------------
% EXAMPLE Trial : jg05-20120312
tind = 23;
Trial = Trials{tind};
stc = Trials{tind}.stc.copy();
exunit = 18;
% <<< EXAMPLE Trial : jg05-20120312 <<< ---------------------------------------

% >>> LOAD rate maps >>> ------------------------------------------------------
% LOAD rate maps --------------------------------------------------------------
% LOAD place restricted behavior fields
Bfs      = cf(@(T,U)  compute_bhv_ratemaps(T,U),          Trials, Units);
bfsEx    = bfs{tind};
BfsShuff = cf(@(T,U)  compute_bhv_ratemaps_shuffled(T,U), Trials, Units);

% LOAD all place fields
pfst = cf(@(T,U)  pfs_2d_theta(T,U),    Trials,Units);
pfss = cf(@(T,U)  pfs_2d_states(T,U),   Trials,Units);
pfsr = cf(@(T,U)  pfs_2d_states(T,U,'', ...
                                {'rear&theta','hbhv&theta','lbhv&theta'}),  ...
          Trials,Units);
pfsa = cf(@(s,t)  cat(2,{t},s),       pfss,pfst);

% <<< LOAD rate maps <<< ------------------------------------------------------

% >>> COMPUTE bfs erpPCA >>> --------------------------------------------------

[eigVecs, eigScrs, eigVars, unitSubset, validDims, zrmMean, zrmStd] = ...
                    compute_bhv_ratemaps_erpPCA(Bfs, Units, [], [], true);
numComp = size(eigVecs,2);
fpc  = cell([1,numComp]);
for i = 1:numComp,
    fpc{i} = nan(size(validDims));
    fpc{i}(validDims) = eigVecs(:,i);
end
fpcLims = [min(cellfun(@min,fpc)),max(cellfun(@max,fpc))];
% RESHAPE eigenvectors
bhvMask = false(size(validDims));
bhvMask(validDims) = true;
bhvMask = reshape_eigen_vector(bhvMask,bfs)';
% SET behavior space visual limits
bhvLims = [-1.7, 0.5; ...
           -0.5, 1.7];

% <<< COMPUTE bfs erpPCA <<< --------------------------------------------------

% >>> SELECT example rate maps >>> --------------------------------------------
pft = pfst{tind};
pft.parameters.states = 'theta-groom-sit';
% LOAD behavior place field
pfs = pfss{tind};
% SORT place field states to match states
pfStates = cf(@(s) ['^',s,'$'],cf(@(p) p.parameters.states,pfs));
for s = 1:numel(pfStates),
    pfStates{s} = strrep(pfStates{s},'&','[&]');
    %pfStates{s} = strrep(pfStates{s},'-','[-]');
end    
for s = 1:numStates,
    psi(s) = find(~cellfun(@isempty,...
                           regexp( repmat(states(s), size(pfStates)),...
                                  pfStates)));
end
pfs = pfs(psi);
% <<< SELECT example rate maps <<< --------------------------------------------

t = 23;
figure,plot(pfsr{t}{3},18);


% >>> ACCUMULATE place fields >>> ---------------------------------------------
accumulate_patch_stats();
% <<< ACCUMULATE place fields <<< ---------------------------------------------

t = 23;
figure,
for u = Units{t};
    uid = find(ismember(cluSessionMap,[t,u],'rows'));
    mrate = max(apfstats.peakFR(uid,:));
    clf();
    for s = 1:3
        subplot(1,3,s);
        plot(pfsr{t}{s},u,1,'colorbar',[0,mrate]);
        hold('on');
        for p = 1:2
            if apfstats.patchPFR(uid,s,p)>rthresh
                rinds = sq(apfstats.patchRateInd(uid,s,p,:,:))';
                rinds = rinds(nniz(rinds),:);
                if size(rinds,1)>athresh/400
                    plot(pfsr{t}{s}.adata.bins{1}(rinds(:,1)),...
                         pfsr{t}{s}.adata.bins{2}(rinds(:,2)),...
                         'm*');                    
                end        
            end
        end
        for p = 1:sum(nniz(patchCntrF(uid,:,1)'))
            circle(patchCntrF(uid,p,1),patchCntrF(uid,p,2),150);
        end
        title(num2str(u))
    end
     waitforbuttonpress();
end




%% Now that I have patches
% get patch max firing rate across states {'rear','high','low'}
% Why?
gridPoints = cell([1,2]);
[gridPoints{:}] = ndgrid(pfsr{1}{1}.adata.bins{:});
gridPoints = cat(3,gridPoints{:});
patchPFRF = nan([numUnits,numStates,numPatches]);
patchMFRF = nan([numUnits,numStates,numPatches]);
patchCNTF = nan([numUnits,1]);
for uid = 1:numUnits
    t = cluSessionMap(uid,1);
    u = cluSessionMap(uid,2);
    for p = 1:numPatches
        if ~isnan(patchCntrF(uid,p,1))
            gridInds = find(sqrt(sum((gridPoints-repmat(patchCntrF(uid,p,:),...
                                                   [size(gridPoints,1),size(gridPoints,2),1])).^2,3)) < 150);
            patchPFRF(uid,:,p) = cell2array(cf(@(pfs) max(pfs.data.rateMap(gridInds,pfs.data.clu==u)), pfsr{t}));
            patchMFRF(uid,:,p) = cell2array(cf(@(pfs) mean(pfs.data.rateMap(gridInds,pfs.data.clu==u),'omitnan'), pfsr{t}));
        end
    end
    patchCNTF(uid) = sum(~isnan(patchCntrF(uid,:,1)));
end




uid = find(ismember(cluSessionMap,[tind,exunit],'rows'));            
% $$$ reshape(sq(tempPatchCenter(uid,:,:,:,:)),[],2),'rows')

figure,imagesc(sq(cdist(uid,:,:,1)))
figure,imagesc(sq(cdist(uid,:,:,1)))            

% COMPUTE place field centers for each state
pfssMr = {};
pfssMp = {};
uCount = 0;
for t = 1:numel(Trials),
    [mr,mp] = cf(@(p,u) p.maxRate(Units{t},'interpPar',interpParPfs), pfsa{t});
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


%%%<<< Supplementary figure: inter-state patch peak rates ------------------------------------------
figure,
subplot(131);
    hold('on');
    plot(log10(apfstats.patchPFR(unitSubset,1,1)),log10(apfstats.patchPFR(unitSubset,2,1)),'.');
    line(log10([0.001,30]),log10([0.001,30]),'Color','r');
    line(log10([0.001,30]),log10([0.001,15]),'Color','b');    
    line(log10([0.001,15]),log10([0.001,30]),'Color','b');        
    line(log10([0.001,10]),log10([0.001,30]),'Color','b');            
    line(log10([0.001,30]),log10([0.001,10]),'Color','b');                
    title('rear x high'),
    xlim([-0.75,1.5]);
    ylim([-0.75,1.5]);   
subplot(132);
    hold('on');
    plot(log10(apfstats.patchPFR(unitSubset,1,1)),log10(apfstats.patchPFR(unitSubset,3,1)),'.');
    line(log10([0.001,30]),log10([0.001,30]),'Color','r');
    line(log10([0.001,30]),log10([0.001,15]),'Color','b');    
    line(log10([0.001,15]),log10([0.001,30]),'Color','b');            
    line(log10([0.001,10]),log10([0.001,30]),'Color','b');            
    line(log10([0.001,30]),log10([0.001,10]),'Color','b');                
    title('rear x low');
    xlim([-0.75,1.5]);
    ylim([-0.75,1.5]);   
subplot(133);
    hold('on');
    plot(log10(apfstats.patchPFR(unitSubset,2,1)),...
         log10(apfstats.patchPFR(unitSubset,3,1)),'.');
    line(log10([0.001,30]),log10([0.001,30]),'Color','r');
    line(log10([0.001,30]),log10([0.001,15]),'Color','b');    
    line(log10([0.001,15]),log10([0.001,30]),'Color','b');        
    line(log10([0.001,10]),log10([0.001,30]),'Color','b');            
    line(log10([0.001,30]),log10([0.001,10]),'Color','b');                
    title('high x low')
    xlim([-0.75,1.5]);
    ylim([-0.75,1.5]);   

figure,
    hold('on');
    plot(apfstats.patchPFR(unitSubset,1,1),apfstats.patchPFR(unitSubset,2,1),'.k');
    line([0.1,30],[0.1,30],'Color','r');
    line([0.1,30],[0.1,15],'Color','b');    
    line([0.1,15],[0.1,30],'Color','b');        
    line([0.1,10],[0.1,30],'Color','b');            
    line([0.1,30],[0.1,10],'Color','b');                
    plot(apfstats.patchPFR(unitSubset,1,1),apfstats.patchPFR(unitSubset,3,1),'rs');
    plot(apfstats.patchPFR(unitSubset,2,1),apfstats.patchPFR(unitSubset,3,1),'g^');
% $$$ figure;
% $$$ subplot(131);
% $$$     hist(log2(apfstats.patchPFR(unitSubset,1,1)./apfstats.patchPFR(unitSubset,2,1)),linspace(-6,6,26));    
% $$$ subplot(132);
% $$$     hist(log2(apfstats.patchPFR(unitSubset,1,1)./apfstats.patchPFR(unitSubset,3,1)),linspace(-6,6,26));    
% $$$ subplot(133);
% $$$     hist(log2(apfstats.patchPFR(unitSubset,2,1)./apfstats.patchPFR(unitSubset,3,1)),linspace(-6,6,26));    

% $$$ figure();
% $$$ hold('on');
% $$$ plot(pfsr{t}{1},17,1,'colorbar');
% $$$ rinds = sq(pfstats(uid,1).patchRateInd(1,1,1,:,:))';
% $$$ rinds = rinds(nniz(rinds),:);
% $$$     plot(pfsr{t}{s}.adata.bins{1}(rinds(:,1)),...
% $$$          pfsr{t}{s}.adata.bins{2}(rinds(:,2)),...
% $$$          'm*');
%circle(apfstats.patchCenters(uid,s,p,1),apfstats.patchCenters(uid,s,p,2),150);

%%%>>>


%%%<<< LOAD Behavioral state contours --------------------------------------------------------------
[stateContourMaps,stateContourHandles] =                           ...
    bhv_contours('MjgER2016',                                  ... sessionListName
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

%%%<<< LOAD and TRANSFORM(tSNE) erpPCA scores ------------------------------------------------------

[fsrcz,FSrC,rmaps,FSCFr,FSrM,FSrS,fsrsMean,fsrsStd,rmapsShuffledMean,rmapsShuffled] = ...
    compute_bhv_ratemaps_erpPCA_scores(Trials,Units,Bfs,BfsShuff,eigVecs,validDims,unitSubset,true);    

% GET auxiliary features
clu =  cf(@(B,U) B.data.clu(:,ismember(B.data.clu,U),:), Bfs,Units);    
si  =  cf(@(B,U) B.data.si(:,ismember(B.data.clu,U),:),  Bfs,Units);
tlu =  cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(Units),1,ones([1,numel(Units)])),Units);
si  = cat(2, si{:});   
clu = cat(2, clu{:});    
tlu = cat(2, tlu{:});    
clu = [tlu',clu'];

% REORDER clusters
[~,rind] = sortrows(clu);
si  = si(:,rind);

bsi = sum(1/size(rmaps,1).*bsxfun(@rdivide,rmaps,mean(rmaps)) ...
          .*log2(bsxfun(@rdivide,rmaps,mean(rmaps))),'omitnan')';


sigUnits = any(abs(fsrcz(:,1:3))>=1.96,2);
sigUnits = any(abs(fsrcz(:,1:3))>=2.58,2);
sigUnitsSB = sum(fsrcz(:,1:3)>=2.58,2)==1;
sigUnitsDB = sum(fsrcz(:,1:3)>=2.58,2)==2;
sigUnitsTB = sum(fsrcz(:,1:3)>=2.58,2)==3;
sigUnitsNB = sum(abs(fsrcz(:,1:3))>=2.58,2)==0;

sum(fsrcz(:,2)<=-2.58 &fsrcz(:,1)>=2.58)./numel(unitSubset)
sum(fsrcz(:,3)<=-2.58 &fsrcz(:,1)>=2.58)./numel(unitSubset)
sum(fsrcz(:,3)<=-2.58 & fsrcz(:,2)<=-2.58 &fsrcz(:,1)>=2.58)./size(fsrcz,1)
sum(fsrcz(:,1)<=-2.58 & fsrcz(:,2)<=-2.58 &fsrcz(:,3)>=2.58)./size(fsrcz,1)
sum(fsrcz(:,1)<=-2.58 & fsrcz(:,3)<=-2.58 &fsrcz(:,2)>=2.58)./size(fsrcz,1)

find(fsrcz(:,1)<=-2.58 & fsrcz(:,3)<=-2.58 &fsrcz(:,2)>=2.58)
find(fsrcz(:,1)<=-2.58 & fsrcz(:,2)<=-2.58 &fsrcz(:,3)>=2.58)

sigUnitsAB = sum(fsrcz(:,1:3)<=-2.58,2)==1;
sum(sigUnitsAB)./numel(sigUnitsNB)

sum(sigUnitsSB)./numel(sigUnitsNB)
sum(sigUnitsDB)./numel(sigUnitsNB)
sum(sigUnitsTB)./numel(sigUnitsNB)
(sum(sigUnitsSB) + sum(sigUnitsDB) + sum(sigUnitsNB))./numel(sigUnitsNB) 
(sum(sigUnitsDB))./numel(sigUnitsNB) 

find(sigUnitsTB)
u = 281;
figure,plot(bfs{cluSessionMap(unitSubset(u),1)},cluSessionMap(unitSubset(u),2),1,'text',[],false);

cluSessionSubset = cluSessionMap(unitSubset,:);
% $$$ for u = cluMap'
% $$$     uind = find(ismember(cluSessionSubset,u','rows'));
% $$$     sigUnits(uind);
% $$$ end
colorOrder = [1,2,3];

cc = eigScrs(:,colorOrder)+0.75;
cc(~sigUnits,:) = repmat([0.75,0.75,0.75],[sum(~sigUnits),1]);

if ~exist('mapa','var'),
    %mapa = tsne([FSrC(:,1:3),(si(unitSubset)'-mean(si(unitSubset))')./std(si(unitSubset))'],[],2,3,40);    
    mapn = tsne([FSrC(:,1:3)],[],2,3,15);    
    %figure,scatter(mapn(:,1),mapn(:,2),10,cc,'Filled')    
    mapa = tsne([FSrC(:,1:3)],[],2,3,40);
    mapa = tsne([FSrC(:,1:3)],[],2,3,15);

figure();
mi = [1,2];
%scatter3(FSrC(:,1),FSrC(:,3),FSrC(:,2),5,cc,'o','filled'); 
scatter(mapa(:,mi(1)),mapa(:,mi(2)),10,cc,'o','filled'); 
    
% $$$     mapa = tsne([fsrcz(:,1:3),(si(unitSubset)'-mean(si(unitSubset))')./std(si(unitSubset))'],[],2,4,80);
% $$$     mapa = tsne([fsrcz(:,1:3),si(unitSubset)'],[],2,3,50);
end
%%%>>>



%%%<<< SET Figure Opts
cond_round = @(rate) max([round(rate,0),round(rate,1)].*[rate>=10,rate<10]);
nanColor = [0.25,0.25,0.25];

% SET figure opts
[hfig,fig,fax,sax] = set_figure_layout(figure(666007),'A4','portrait',[],1.5,1.5,0.05,0.05);
% $$$ pageWidth  = 21.;
% $$$ pageHeight = 29.7;
% $$$ pwidth = 1.5;
% $$$ pheight = 1.5;
% $$$ xpad = 0.05;
% $$$ ypad = 0.05;
% $$$ xpos = 3.5:(pwidth+xpad):pageWidth;
% $$$ ypos = fliplr(0.5:(pheight+xpad):pageHeight-3.7);
%%%>>>

% $$$ %%%<<< BehaviorPlaceCode-F2-Sup: behavior and place field examples -> move to supfig
% $$$ 
% $$$ cluMap = [20,119];
% $$$ yind = 1;
% $$$ xind = 1;
% $$$ 
% $$$ u = cluMap';
% $$$ % SET color scale max
% $$$ maxPfsRate = max(cell2mat(cf(@(p,u) maxRate(p,u,false,'prctile99',0.5),...
% $$$                              {pfs{1},bfsEx},repmat({u(2)},[1,2]))));
% $$$ 
% $$$ 
% $$$ % ADJUST subplot coordinates
% $$$ [yind, yOffSet, xind, xOffSet] = deal(1,0, 1, 0);
% $$$ % CREATE subplot axes
% $$$ sax(end+1) = axes('Units','centimeters',                                        ...
% $$$                   'Position',[fig.page.xpos(xind)+xOffSet,                      ...
% $$$                               fig.page.ypos(yind)+yOffSet,                      ...
% $$$                               fig.subplot.width.*1.25,                          ...
% $$$                               fig.subplot.height.*1.25],                        ...
% $$$                   'FontSize', 8,                                                ...
% $$$                   'LineWidth',1);
% $$$ hold(sax(end),'on');
% $$$ bfs{tind}.plot(u(2),'mean',false,[0,maxPfsRate],true,0.5,false,[],@jet,bhvMask,nanColor);
% $$$ text(-2,-0.45,num2str(cond_round(bfs{tind}.maxRate(u(2),false,1))),'FontSize',10,'Color',[1,1,1]);
% $$$ xlabel({'HP (rad)'});
% $$$ ylabel({'BP (rad)'});
% $$$ xind = xind+2;
% $$$ xlim([-2,nonzeros(xlim.*[0,1])-0.2])
% $$$ sax(end).Color = nanColor;
% $$$ title({'Behavior','Rate Map'});    
% $$$ 
% $$$ [yind, yOffSet, xind, xOffSet] = deal(1,0, 3, 0);
% $$$ % CREATE subplot axes
% $$$ sax(end+1) = axes('Units','centimeters',                                        ...
% $$$                   'Position',[fig.page.xpos(xind)+xOffSet,                      ...
% $$$                               fig.page.ypos(yind)+yOffSet,                      ...
% $$$                               fig.subplot.width.*1.25,                          ...
% $$$                               fig.subplot.height.*1.25],                        ...
% $$$                   'FontSize', 8,                                                ...
% $$$                   'LineWidth',1);
% $$$ hold(sax(end),'on');
% $$$ plot(pfs{1},u(2),'mean',false,[0,maxPfsRate],true,0.5,false,interpParPfs,@jet,[],nanColor);
% $$$ text(-490,-380,num2str(cond_round(maxPfsRate)),'FontSize',10,'Color',[1,1,1]);
% $$$ 
% $$$ %%%<<< ADD contours to ratemap indicating the selected space

hold('on');
rmap = plot(pfs{1},u(2),'mean',false, [0, maxPfsRate], true, 0.5, false, ...
            interpParPfs);
[~,mxp] = Pft.maxRate(u(2));
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
     [sax(end-1).Position(1)+fig.subplot.width.*1.25,sax(end).Position(1)+(contMaxInd(1)./numel(interpParPfs.bins{1})).*fig.subplot.width.*1.25],...
     [sax(end-1).Position(2)+fig.subplot.height.*1.25,sax(end).Position(2)+(contMaxInd(2)./numel(interpParPfs.bins{2})).*fig.subplot.height.*1.25],...
     'Color','m','LineWidth',1);
line(fax,...
     [sax(end-1).Position(1)+fig.subplot.width.*1.25,sax(end).Position(1)+(contMinInd(1)./numel(interpParPfs.bins{1})).*fig.subplot.width.*1.25],...
     [sax(end-1).Position(2),sax(end).Position(2)+(contMinInd(2)./numel(interpParPfs.bins{2})).*fig.subplot.height.*1.25],...
     'Color','m','LineWidth',1);
line(fax,...
     [sax(end-1).Position(1),sax(end-1).Position(1)+fig.subplot.width.*1.25],...
     [sax(end-1).Position(2),sax(end-1).Position(2)],...
     'Color','m','LineWidth',1);      
line(fax,...
     [sax(end-1).Position(1),sax(end-1).Position(1)+fig.subplot.width.*1.25],...
     [sax(end-1).Position(2)+fig.subplot.height.*1.25,sax(end-1).Position(2)+fig.subplot.height.*1.25],...
     'Color','m','LineWidth',1);      
line(fax,...
     [sax(end-1).Position(1)+fig.subplot.width.*1.25,sax(end-1).Position(1)+fig.subplot.width.*1.25],...
     [sax(end-1).Position(2),sax(end-1).Position(2)+fig.subplot.height.*1.25],...
     'Color','m','LineWidth',1);      
line(fax,...
     [sax(end-1).Position(1),sax(end-1).Position(1)],...
     [sax(end-1).Position(2),sax(end-1).Position(2)+fig.subplot.height.*1.25],...
     'Color','m','LineWidth',1);      
uistack(fax,'top');

%%%>>>
% $$$ 
% $$$ title({'Spatial','Rate Map'});   
% $$$ chax = colorbar(sax(end));
% $$$ chax.LineWidth = 1;
% $$$ chax.Units = 'centimeters';
% $$$ colormap(sax(end),'jet');
% $$$ caxis([0,maxPfsRate]);
% $$$ chax.Position(1) = sax(end).Position(1)+fig.subplot.width*1.25+0.1;
% $$$ 
% $$$ sax(end).YTickLabel = {};
% $$$ sax(end).XTickLabel = {};
% $$$ line(fax,fig.page.xpos(xind)+[0,fig.subplot.width*1.25/2],fig.page.ypos(yind).*[1,1]-0.125,'Color','k','LineWidth',2);
% $$$ text(fax,fig.page.xpos(xind),fig.page.ypos(yind)-0.35,'50cm','FontSize',8,'Color',[0,0,0]);
% $$$ 
% $$$ %%%>>>


%%%<<< BehaviorPlaceCode-F2-A: behavior and place field examples

[hfig,fig,fax,sax] = set_figure_layout(figure(666007),'A4','portrait',[],1.75,1.75,0.05,0.05);
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

% two patches [63, 89, 133


cluMap = [23,42;... 74;...
          23,29;... %20,73;...
          23,32;... double
          23,79;... %20,85;...
          23,106;...
          23,137;...
          23,79;...
          23,101];
%20,104];
% ...          20,109;...
%          20,59];

% {accg,BFS{HP,BP},BFS{HP,BS},theta,rear,hloc,lloc,hpause,lpause}


%labels = {{'Behavior','Rate Map'},{'Theta','Rate Map'},'Rear','H Loc','H Pause','L Loc','L Pause'};
labels = {{'Theta','    '},'Rear','H Loc','H Pause','L Loc','L Pause',{'Posture'}};

numStates = numel(states);
% PLOT 
xind = 1;
yind = 1;
yinit = yind;
for u = cluMap',
    uid = find(ismember(cluSessionMap,u','rows'));
    xind = 7;
% SET color scale max for each row -----------------------------------------------------------------
    maxPfsRate = max([cell2mat(cf(@(p,u,m)                                                        ...
                                 p.maxRate(u,true&numel(m)==1,1,[],[],[],m),                     ...
                                 [pfs,{bfsEx}],                                                  ...
                                 repmat({u(2)},[1,numel(pfs)+1]),                                ...
                                 [repmat({true},[1,numel(pfs)]),{bhvMask(:)}])                   ...
                              ),max(nonzeros(rmapP(validDims,uid,:)))]);
    pfsMaxRatesMean = cell2mat(cf(@(p,u) ...
                                  max(p.maxRate(u,true,'mean')),...
                                  pfs,repmat({u(2)},[1,numel(pfs)])));
    bfsMaxRatesMean = max(bfsEx.maxRate(u(2),false,1,[],[],[],bhvMask(:)));


% BHV FIELDS ---------------------------------------------------------------------------------------

    patchInds = find(~isnan(patchCntrF(uid,:,1)));
    for patchId = patchInds
        if xind > 8
            break;
        end
        sax(end+1) = axes('Units','centimeters',                                                 ...
                          'Position',[fig.page.xpos(xind)+0.2*(numStates-3),                     ...
                                      fig.page.ypos(yind),                                       ...
                                      fig.subplot.width,                                         ...
                                      fig.subplot.height],                                       ...
                          'FontSize', 8);
        %bfsEx.plot(u(2),'mean',false,[0,maxPfsRate],false,0.5,false,[],@jet,bhvMask,nanColor);
        rmapEx = rmapP(:,ismember(cluSessionMap,u','rows'),patchId);
        rmapEx(~validDims) = nan;
        imagescnan({bfsEx.adata.bins{:},reshape(rmapEx,bfsEx.adata.binSizes')'},                 ...
                   [0,maxPfsRate],                                                               ...
                   'linear',                                                                     ...
                   false,                                                                        ...
                   nanColor,                                                                     ...
                   [],[],@jet);
        axis('xy');
        sax(end).YTickLabel = {};
        sax(end).XTickLabel = {};        
        % ANNOTATE plot with max rate
        text(0.45,1.475,num2str(round(max(rmapEx))),'FontSize',8,'Color',[1,1,1],'HorizontalAlignment','right');
        xlim(bhvLims(1,:));
        ylim(bhvLims(2,:));
        sax(end).Color = nanColor;
        
        patchColors = 'mw';        
        rectangle(fax,'Position',sax(end).Position+[0.02,0.02,-0.04,-0.04],'EdgeColor',patchColors(patchId));
        uistack(fax,'top');                  
        xind = xind + 1;        

        if yind == yinit && patchId == 1
            line([-1.6,-1.6],...
                 [1.1,1.6],...
                 'LineWidth',2,...
                 'Color','w');
            line([-1.6,-1.1],...
                 [1.6,1.6],...
                 'LineWidth',2,...
                 'Color','w');
        end
    end%for patchId 

% $$$     if yind == yinit, 
% $$$         title(labels{7});
% $$$     end 
    if yind==round(size(cluMap,1)),
        ylabel({'Body Pitch'});
        sax(end).YAxisLocation = 'right';
    end

% DRAW the scale bars for placefields (spatial ratemaps)
    if yind==size(cluMap,1),
        xlabel({'Head-Body','Pitch'});
% $$$         % y-axis scale bar
% $$$         line(fax,...
% $$$              [sax(end).Position([1])].*[1,1]-0.15,...     
% $$$              [sax(end).Position([2]),sum(sax(end).Position([2,4]).*[1,0.5])],...     
% $$$              'LineWidth',1,...
% $$$              'Color',[0,0,0]);
% $$$         % y-axis scale bar label        
% $$$         text(fax,...
% $$$              sax(end).Position(1)-0.35,...
% $$$              sum(sax(end).Position([2,4]).*[1,0.25]),...
% $$$              '1 rad',...
% $$$              'HorizontalAlignment','center',...
% $$$              'VerticalAlignment','middle',...
% $$$              'Rotation',90,...
% $$$              'FontSize',8);
    end
    
    xind = 1;
    for s = 1:numStates
% PLACEFIELDS MTAApfs ------------------------------------------------------------------------------
        sax(end+1) = axes('Units','centimeters',                                                 ...
                         'Position',[fig.page.xpos(xind)+0.2+0.2.*double(s>=2),                  ...
                                     fig.page.ypos(yind),                                        ...
                                     fig.subplot.width,                                          ...
                                     fig.subplot.height],                                        ...
                         'FontSize', 8);
        plot(pfs{s},u(2),'mean',false,[0,maxPfsRate],true,0.5,false,interpParPfs,@jet,[],nanColor);
        sax(end).YTickLabel = {};
        sax(end).XTickLabel = {};
        if yind == yinit, 
            tax = title(labels{xind});
            if s==1,
                %tax.Position(2) = 800;
            end        
        end
                
        xind = xind+1;
        text(490,380,num2str(round(pfsMaxRatesMean(s))),'FontSize',8,'Color',[1,1,1],'HorizontalAlignment','right');
        %text(200,380,num2str(cond_round(pfsMaxRatesMean(s))),'FontSize',8,'Color',[1,1,1]);
        if yind == yinit && s == 1
            line([-480,-480],...
                 [280,480],...
                 'LineWidth',2,...
                 'Color','w');
            line([-480,-280],...
                 [480,480],...
                 'LineWidth',2,...
                 'Color','w');
        end
        
        
        if xind ==2
            hold(sax(end),'on');
            % Circles
            uid = find(ismember(cluSessionMap,u','rows'));
            gpatches = find(~isnan(patchCntrF(uid,:,1)));
            patchColors = 'mw';
            if numel(gpatches)>2,
                gpatches = gpatches(1:2);
            end
            for gp = gpatches
                circle(patchCntrF(uid,gp,1),patchCntrF(uid,gp,2),100,['-',patchColors(gp)]);
            end
% CONTURES 
% $$$             rmap = plot(pfs{1},u(2),'mean',false,[0,maxPfsRate],true,0.5,false,interpParPfs);
% $$$             [~,mxp] = pft.maxRate(u(2));
% $$$             binGrid = cell([1,2]);
% $$$             [binGrid{:}] = ndgrid(interpParPfs.bins{:});
% $$$             binGrid = cat(3,binGrid{:});
% $$$             binDist = sqrt(sum((binGrid-repmat(permute(mxp,[1,3,2]),[cellfun(@numel,interpParPfs.bins),1])).^2,3));
% $$$             [cont,conHax] = contour(binGrid(:,:,1),binGrid(:,:,2),(rmap./max(rmap(:)))>0.2&binDist<250,[0.5,0.5],'m','LineWidth',1);
% $$$ 
% $$$             % REMOVE first coordinates for this example
% $$$             cont(:,1) = [];
% $$$             
% $$$             [contYMax,contYMaxInd] = max(reshape(binGrid(:,:,2).*double((rmap./max(rmap(:)))>0.2&binDist<250),[],1));
% $$$             [contYMin,contYMinInd] = min(reshape(binGrid(:,:,2).*double((rmap./max(rmap(:)))>0.2&binDist<250),[],1));
% $$$             contMaxInd = zeros([1,2]);
% $$$             [contMaxInd(1),contMaxInd(2)] = ind2sub(cellfun(@numel,interpParPfs.bins),contYMaxInd);
% $$$             contMinInd = zeros([1,2]);
% $$$             [contMinInd(1),contMinInd(2)] = ind2sub(cellfun(@numel,interpParPfs.bins),contYMinInd);
% $$$             contMaxPos = sq(binGrid(contMaxInd(1),contMaxInd(2),:));
% $$$             contMinPos = sq(binGrid(contMinInd(1),contMinInd(2),:));
% $$$             
% $$$             [contMaxYord, contMaxInd] = max(cont(2,:));
% $$$             contMaxXord = cont(1,contMaxInd);
% $$$             [contMinYord, contMinInd] = min(cont(2,:));
% $$$             contMinXord = cont(1,contMinInd);
% $$$             
% $$$             [~,contMinXind] = NearestNeighbour(interpParPfs.bins{1},contMinXord);
% $$$             [~,contMinYind] = NearestNeighbour(interpParPfs.bins{2},contMinYord);            
% $$$             [~,contMaxXind] = NearestNeighbour(interpParPfs.bins{1},contMaxXord);
% $$$             [~,contMaxYind] = NearestNeighbour(interpParPfs.bins{2},contMaxYord);            

% $$$             if yind == 1,
% $$$              % TOP line                
% $$$              line(fax,...
% $$$                   [sax(end-1).Position(1)+fig.subplot.width, sax(end).Position(1)+(contMaxXind./numel(interpParPfs.bins{1})).*fig.subplot.width],...
% $$$                   [sax(end-1).Position(2)+fig.subplot.height,sax(end).Position(2)+(contMaxYind./numel(interpParPfs.bins{2})).*fig.subplot.height],...
% $$$                   'Color','m','LineWidth',1);
% $$$             % BOTTOM line            
% $$$              line(fax,...
% $$$                   [sax(end-1).Position(1)+fig.subplot.width, sax(end).Position(1)+(contMinXind./numel(interpParPfs.bins{1})).*fig.subplot.width],...
% $$$                   [sax(end-1).Position(2),sax(end).Position(2)+(contMinYind./numel(interpParPfs.bins{2})).*fig.subplot.height],...
% $$$                   'Color','m','LineWidth',1);
% $$$              rectangle(fax,'Position',sax(end-1).Position+[-0.02,-0.02,0.04,0.04],'EdgeColor','m');
% $$$             end% if yind==1
% $$$             uistack(fax,'top');
        end% if xind==2
% $$$         if  yind == 1 && s == numStates,
% $$$             cax = colorbar(sax(end));
% $$$             colormap(cax,'jet');
% $$$             cax.Units = 'centimeters';
% $$$             cax.Position = [sum(sax(end).Position([1,3]))+0.15,sax(end).Position(2),0.15,sax(end).Position(4)];
% $$$             cax.XTick = [0,1];
% $$$             cax.XTickLabel = {'0','Max'};
% $$$             cax.Label.String = 'Hz';
% $$$             cax.Label.Units = 'centimeters';
% $$$             cax.Label.FontSize = 8;
% $$$             cax.Label.Position(1) = 0.3;
% $$$         end
    end% for s
    yind = yind+1;
end


% DRAW line and arrows for theta to behavioral state decomposition
s = 7*size(cluMap,1)-3;
%sax(end-s:end-s+5)
% HORIZONTAL line
line(fax,                                                                                         ...
     [sax(3).Position(1)+fig.subplot.width-0.2,sax(8).Position(1)+fig.subplot.width/2+0.07],                  ...
      repmat(sax(3).Position(2)+sax(3).Position(4)+0.575,[1,2]),                            ...
     'Color','k',                                                                                 ...
     'LineWidth',1);
for a = 1:5,
    patch(fax,                                                                                    ...
          repmat(sax(3+a).Position(1)+fig.subplot.width/2,[1,3])+[-0.07,0.07,0],                        ...
          repmat(sax(3+a).Position(2)+sax(3+a).Position(4)+0.475,[1,3])+[0.1,0.1,-0.09],        ...
          'k');
end

% $$$ % DRAW the scale bars for placefields (spatial ratemaps)
% $$$ % x-axis scale bar
% $$$ line(fax,...
% $$$      [sum(sax(end-numStates+1).Position([1,3]).*[1,0.5]),sum(sax(end-numStates+1).Position([1,3]))],...
% $$$      [sax(end-numStates+1).Position([2]),sax(end-numStates+1).Position([2])]-0.15,...
% $$$      'LineWidth',1,...
% $$$      'Color',[0,0,0]);
% $$$ text(fax,...
% $$$      sum(sax(end-numStates+1).Position([1,3]).*[1,0.75]),...
% $$$      sax(end-numStates+1).Position([2])-0.35,...
% $$$      '50cm',...
% $$$      'HorizontalAlignment','center',...
% $$$      'VerticalAlignment','middle',...
% $$$      'FontSize',8);
% $$$ 
% $$$ % y-axis scale bar
% $$$ line(fax,...
% $$$      [sax(end-numStates+1).Position([1]),sax(end-numStates+1).Position([1])]-0.15,...     
% $$$      [sax(end-numStates+1).Position([2]),sum(sax(end-numStates+1).Position([2,4]).*[1,0.5])],...     
% $$$      'LineWidth',1,...
% $$$      'Color',[0,0,0]);
% $$$ text(fax,...
% $$$      sum(sax(end-numStates+1).Position([1]))-0.35,...
% $$$      sum(sax(end-numStates+1).Position([2,4]).*[1,0.25]),...
% $$$      '50cm',...
% $$$      'HorizontalAlignment','center',...
% $$$      'VerticalAlignment','middle',...
% $$$      'Rotation',90,...
% $$$      'FontSize',8);



%%%>>> 


%%%<<< BehaviorPlaceCode-F2-Sup: plot placefield centers -> move to supfigure
% $$$ csteps = floor(2*pi*480);
% $$$ grid   = [1:csteps];
% $$$ angle  = grid*2*pi/csteps;
% $$$ xind = 2;
% $$$ for s = 1:6,
% $$$         sax(end+1) = axes('Units','centimeters',...
% $$$                          'Position',[fig.page.xpos(xind)+0.2+0.2.*double(s>=2),fig.page.ypos(yind)-0.2,fig.subplot.width,fig.subplot.height],...
% $$$                          'FontSize', 8,...
% $$$                          'Color', nanColor);
% $$$         hold('on');
% $$$         % DRAW disk
% $$$         pax = patch( sin(angle) * 480, cos(angle) * 480, [0.8,0.8,0.8]);        
% $$$         sind = pfssMrSub(:,1,s)>3;
% $$$         plot(pfssMpSub(sind,1,s)+randn([sum(sind),1])*10,...
% $$$              pfssMpSub(sind,2,s)+randn([sum(sind),1])*10,'.b');
% $$$         xind = xind +1;
% $$$         sax(end).YTickLabel = {};
% $$$         sax(end).XTickLabel = {};        
% $$$ end



%%%<<< BehaviorPlaceCode-BhvPlaceFields: patch cntyind = size(cluMap,1)+4;
yind = size(cluMap,1)+2;
sax(end+1) = axes('Units','centimeters',                     ...
                  'Position',[fig.page.xpos(1),              ...
                              fig.page.ypos(yind),           ...
                              fig.subplot.width*1.5,             ...
                              fig.subplot.height*1.5],           ...
                  'FontSize', 8);    
patchCntF = sum(~isnan(patchCntrF(unitSubset,:,1)),2);

%out = histcounts(patchCntF,0.5:4.5,'Normalization','probability');
anatCA1Id = [4:6,14:15,21:28,30];
anatCA3Id = [1:3,7:13,16:19,29];

outCA1 = histcounts(patchCntF(ismember(cluSessionSubset(:,1),anatCA1Id)),0.5:4.5,'Normalization','count');
outCA3 = histcounts(patchCntF(ismember(cluSessionSubset(:,1),anatCA3Id)),0.5:4.5,'Normalization','count');


bar(1:4,[outCA1;outCA3],'stacked')
legend({'CA1','CA3'});
xlabel('# of Patches');
ylabel('Count');
xlim([0.35,4.65]);

% $$$ unitSubset
% $$$ pieh = pie(out);
% $$$ pieLbls = pieh(2:2:end);
% $$$ extraS = ' sss';
% $$$ for pl = 1:numel(pieLbls)-1,
% $$$     %pieLbls(pl).String = [num2str(pl),' field',extraS(pl),': ',pieLbls(pl).String];
% $$$     pieLbls(pl).String = [num2str(pl)];
% $$$     pieLbls(pl).Position = pieLbls(pl).Position./3+[0,0.05,0];
% $$$     pieLbls(pl).Color = [1,1,1];
% $$$     pieLbls(pl).FontSize = 10;
% $$$     pieLbls(pl).FontWeight = 'bold';
% $$$ end
% $$$ pieLbls(pl+1).String = [num2str(pl+1)];
% $$$ pieLbls(pl+1).FontWeight = 'bold';

%%%>>>

%%%<<< BehaviorPlaceCode-BhvPlaceFields: inter pfs patch correlation
yind = size(cluMap,1)+2;
sax(end+1) = axes('Units','centimeters',                     ...
                  'Position',[fig.page.xpos(3)+0.5,              ...
                              fig.page.ypos(yind),           ...
                              fig.subplot.width*1.5,             ...
                              fig.subplot.height*1.5],           ...
                  'FontSize', 8);    
plot(rmapIPDist/10,rmapIPCorr,'.')
xlabel('Patch Distance (cm)')
ylabel('Bhv Ratemap Corr');
%%%>>>


%%%<<< BehaviorPlaceCode-BhvPlaceFields: inter pfs patch correlation
yind = size(cluMap,1)+2;
sax(end+1) = axes('Units','centimeters',                     ...
                  'Position',[fig.page.xpos(5)+0.5,              ...
                              fig.page.ypos(yind),           ...
                              fig.subplot.width*1.5,             ...
                              fig.subplot.height*1.5],           ...
                  'FontSize', 8);    
hold('on');
ph = patch([10.1,10.1,20.0,20.0],[-1,1,1,-1],[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',0.3);
pd = reshape(patchDistF(unitSubset,:,:),[],1)/10;
rcorr = reshape(rmapCorr(1,unitSubset,:,:),[],1);
nind = nniz(pd) & nniz(rcorr) & pd>10;
plot(pd(nind),rcorr(nind),'.');
box('on');
xlabel('Patch Distance (cm)')
ylabel('');
xlim([10,75])
sax(end).YTickLabel = {};
%%%>>>


%%%<<< BehaviorPlaceCode-F2-B: behavior field erpPCA eigenvectors
% PLOT 
yind = size(cluMap,1)+4;
for i = 1:3,
    sax(end+1) = axes('Units','centimeters',                    ...
                     'Position',[fig.page.xpos(1),              ...
                                 fig.page.ypos(yind),           ...
                                 fig.subplot.width,             ...
                                 fig.subplot.height],           ...
                     'FontSize', 8);    
    imagescnan({bfsEx.adata.bins{:},abs(reshape_eigen_vector(fpc{i},{bfsEx}))},...
               fpcLims,'linear',false,nanColor,1,1);                % PRINT eigenvectors
    axis('xy');
    axis('tight');
    hold('on');    
    for s = 1:numel(stateContourHandles),                            % OVERLAY state Contours
        copyobj(stateContourHandles{s},sax(end));
    end
    sax(end).YTickLabel = {};
    sax(end).XTickLabel = {};
    sax(end).Color = nanColor;
    ylabel(['F',num2str(i)])
    xlim(bfsEx.adata.bins{1}([1,end]));
    xlim([-2,nonzeros(xlim.*[0,1])-0.2])            
    ylim(bfsEx.adata.bins{2}([1,end]));
    yind = yind + 1;
    if i == 1,
        title({'erpPCA','Factors'});
    end
end
af(@(ax) set(ax,'LineWidth',1), sax);

%%%>>>

%%%<<< BehaviorPlaceCode-F2-D: t-SNE mapping of fscores within HPxBP of first 3 eigenvectors

% DATA SOURCE:  compute_bhv_ratemaps_erpPCA_scores
gca();
yind = yind-1;
xind = 2;
sax(end+1) = axes('Units','centimeters',...
                 'Position',[fig.page.xpos(xind)+0.2,fig.page.ypos(yind),(fig.subplot.width+0.1)*3-0.1,(fig.subplot.height+0.1)*3-0.1],...
                 'FontSize', 8);    
mi = [1,2];
scatter(mapa(:,mi(1)),mapa(:,mi(2)),10,cc,'o','filled'); 
sax(end).YTickLabel = {};
sax(end).XTickLabel = {};
box('on');
axis('tight')
xlim(xlim.*1.05);
ylim(ylim.*1.05);

set(gca(),'Color',nanColor)
title({'t-SNE Transformed','erpPCA F-Scores'});
%%%>>>




%%%<<< BehaviorPlaceCode-F2-E: ECDF of eigenvector significances
yind = yind-1;
xind = 5;
sax(end+1) = axes('Units','centimeters',...
                 'Position',[fig.page.xpos(xind)+0.4,...
                             fig.page.ypos(yind)+(fig.subplot.height+0.1)*0.5,...
                            (fig.subplot.width+0.1)*3-0.1,...
                            (fig.subplot.height+0.1)*1.5],...
                 'FontSize', 8);    

hold('on');
sclr = 'rgb';
sclr = sclr(colorOrder);
for i = [1,2,3],
    [F,X] = ecdf(fsrcz(:,i));
    cax = cdfplot(X);
    cax.Color = sclr(i);
end
set(sax(end),'YAxisLocation','right')
Lines(-1.96,[],'k');
Lines(1.96,[],'k');
Lines(-3.1,[],'k');
Lines(3.1,[],'k');
legendLabels = {'rear','high','low'};
lax = legend(legendLabels(colorOrder),'Location','southoutside');
lax.Units = 'centimeters';
xlim([-20,20])
lax.Position(1) = sax(end).Position(1);
drawnow();
lax.Position(2) = sax(end).Position([2])-1.75;
sax(end).YTick  = [0,0.2,0.4,0.6,0.8,1];
%%%>>>



% END FIGURE ---------------------------------------------------------------
% $$$ 
% $$$ RESULTS
% $$$ 
% $$$ FACTOR ANALYSIS OF THE RATE MAPS OF THE POSTURE SPACE RESULTED IN # PRIMARY 
% $$$ FACTORS WHICH CORRESPOND WELL WITH REAR< HIGH< AND LOW BEHAVIORS


% low loc cells
sum(fsrcz(:,3)>2.58 & fsrcz(:,2)<2.58 & fsrcz(:,1)<2.58)/size(fsrcz,1)
% 12.4 percent 61 of 492
% high loc cells
sum(fsrcz(:,2)>2.58 & fsrcz(:,3)<2.58 & fsrcz(:,1)<2.58)/size(fsrcz,1)
sum(fsrcz(:,1)>2.58 &  fsrcz(:,3)<2.58)/size(fsrcz,1)

1-(1-0.99).^(1/size(fsrcz,1))


u = find((fsrcz(:,3)>2.58 & fsrcz(:,2)<2.58 & fsrcz(:,1)<2.58));

anatomicalLocation = [0,0,0,                        ... er01
                      1,1,1,0,0,                    ... ER06
                      0,0,0,0,0,                    ... Ed10
                      1,1,0,0,0,0,                  ... jg04                      
                      1,1,1,1,1,1,1,1,1,0           ... jg05
                      1];% new units - FS03


anatomicalLocation = [0,0,                          ... er01
                      1,1,1,                        ... ER06
                      0,0,                          ... Ed10
                      0,0,0,0,0,1,1,1,1             ... jg04                      
                      1,1,1,1,1,1,1,                ... jg05
                      1,1,0,0,0];% new units - jg05, jg05, ER06, Ed10, er01                    


inCA1 = anatomicalLocation(cluSessionSubset(:,1))==1;
inCA3 = ~inCA1;


% 3.6% (14 of 380) CA1 cells were restricted to rearing exploration
sum(fsrcz(inCA1,1)>2.58 & fsrcz(inCA1,2)<-2.58 & fsrcz(inCA1,2)<-2.58)/sum(inCA1)
% 4.2% (16 of 380) CA1 cells were 
sum(fsrcz(inCA1,3)>2.58 & fsrcz(inCA1,2)<-2.58 & fsrcz(inCA1,1)<-2.58)/sum(inCA1)
% 5.3% (20 of 380) CA1 cells were 
sum(fsrcz(inCA1,2)>2.58 & fsrcz(inCA1,3)<-2.58 & fsrcz(inCA1,1)<-2.58)/sum(inCA1)

sum(find(fsrcz(inCA1,1)<-2.58))/sum(inCA1)


figure,bar(linspace(0,40,40),histc(sqrt(sum(fsrcz(:,1:3).^2,2)),linspace(0,40,40)),'histc');


% 7.1% (27 of 380) CA1 cells were selective for rearing exploration
sum(fsrcz(inCA1,1)>2.58 & fsrcz(inCA1,2)<2.58 & fsrcz(inCA1,3)<2.58)/sum(inCA1)
% 22.1% (84 of 380) CA1 cells were selective high exploration
sum(fsrcz(inCA1,2)>2.58 & fsrcz(inCA1,1)<2.58 & fsrcz(inCA1,1)<2.58)/sum(inCA1)
% 12.4% (47 of 380) CA1 cells were selective for low exploration
sum(fsrcz(inCA1,3)>2.58 & fsrcz(inCA1,2)<2.58 & fsrcz(inCA1,1)<2.58)/sum(inCA1)


sum(fsrcz(inCA1,1)<-2.58)/sum(inCA1)

% 4.5% (5 of 112) CA1 cells were selective for rearing exploration
sum(fsrcz(inCA3,1)>2.58 & fsrcz(inCA3,2)<2.58 & fsrcz(inCA3,3)<2.58)/sum(inCA3)
% 17.9% (20 of 112) CA3 cells are dominantly active during high exploration
sum(fsrcz(inCA3,2)>2.58 & fsrcz(inCA3,1)<2.58 & fsrcz(inCA3,1)<2.58)/sum(inCA3)
% 12.5% (14 of 112) CA3 cells are dominantly active during low exploration
sum(fsrcz(inCA3,3)>2.58 & fsrcz(inCA3,2)<2.58 & fsrcz(inCA3,1)<2.58)/sum(inCA3)



u = find(fsrcz(:,1)<-2.58 & fsrcz(:,2)<2.58 & fsrcz(:,3)<2.58);
u = 407;
figure();
for i = u',
    plot(bfs{cluSessionSubset(i,1)},cluSessionSubset(i,2),'mean','text',[],true,'mazeMask',bhvMask);
    title(num2str(cluSessionSubset(i,:)));
    waitforbuttonpress();
end

cluSessionSubset(i,:)
fsrcz(find(ismember(cluSessionSubset,cluSessionSubset(i,:),'row')),:)

sum(ismember(cluSessionSubset(:,1),[9:12]))

% CA1
% jg04(15), ER06(103), jg05(255)
% 
% CA3
% er01(26), Ed10(36), jg04(22), ER06(33)

bsi = sum(1/size(rmaps,1).*(rmaps(:,225)./mean(rmaps(:,225))).*log2(rmaps(:,225)./mean(rmaps(:,225))),'omitnan');

bsio = sum(1/size(rmaps,1).*bsxfun(@rdivide,rmaps,mean(rmaps)) ...
          .*log2(bsxfun(@rdivide,rmaps,mean(rmaps))),'omitnan');

brmaps = rmaps;
brmaps(brmaps==0)=nan;

bsi = sum(repmat(1./sum(~isnan(brmaps(:,:,1))),[size(brmaps,1),1]).*bsxfun(@rdivide,brmaps,mean(brmaps,'omitnan')) ...
          .*log2(bsxfun(@rdivide,rmaps,mean(brmaps,'omitnan'))),'omitnan');


figure();
plot3(fsrcz(:,1),fsrcz(:,2),fsrcz(:,3),'.');

plot3(FSrC(:,1),FSrC(:,2),FSrC(:,3),'.');


figure();
plot(log10(max(fsrcz(:,1:3),[],2)),bsi,'.');
plot(log2(sqrt(sum(fsrcz(:,1:3).^2,2))),bsi,'.');

find(sqrt(sum(fsrcz(:,1:3).^2,2))<1&bsi'>1.4)

bxyv = fet_bref(Trial);

figure,
bar(linspace(0,40,40),...
    histc(sqrt(sum(fsrcz(:,1:3).^2,2)),linspace(0,40,40)),'histc');

rmapsShuffled = decapsulate_and_concatenate_mtaapfs(bfsShuff,units);
rmapsShuffled = rmapsShuffled(validDims,unitSubset,:);
rmapsShuffled(~nniz(rmapsShuffled(:))) = 0;    
bsiShuffled = sq(sum(1/size(rmapsShuffled,1).*bsxfun(@rdivide,rmapsShuffled,mean(rmapsShuffled)) ...
          .*log2(bsxfun(@rdivide,rmapsShuffled,mean(rmapsShuffled))),'omitnan'));

bsiShuffled = sq(sum(repmat(1./sum(~isnan(rmapsShuffled(:,:,1))),[size(rmapsShuffled,1),1,size(rmapsShuffled,3)]) ...
                     .*bsxfun(@rdivide,rmapsShuffled,mean(rmapsShuffled,'omitnan')) ...
          .*log2(bsxfun(@rdivide,rmapsShuffled,mean(rmapsShuffled,'omitnan'))),'omitnan'));

bsiShuffled = sq(sum(1/size(rmapsShuffled,1),'omitnan').*bsxfun(@rdivide,rmapsShuffled,mean(rmapsShuffled,'omitnan')) ...
          .*log2(bsxfun(@rdivide,rmapsShuffled,mean(rmapsShuffled))),'omitnan'));

bsiZ = (log(bsi) - mean(log(bsiShuffled),2))./std(log(bsiShuffled),[],2);
bsiZ = (log2(bsi') - mean(log2(bsiShuffled),2))./std(log2(bsiShuffled),[],2);
bsiZ = ((bsi) - mean((bsiShuffled),2))./std((bsiShuffled),[],2);


(log2(bsi(i)) - mean(log2(bsiShuffled(i,:)),2))./std(log2(bsiShuffled(i,:)),[],2)
 
u = 1:numel(unitSubset);
figure();
for i = u,
    subplot2(2,3,1,1);
        plot(bfs{cluSessionSubset(i,1)},...
             cluSessionSubset(i,2),...
             'mean','colorbar',[],true,...
             'mazeMask',bhvMask,...
             'colorMap',@jet);
        title(num2str(cluSessionSubset(i,:)));
    subplot2(2,3,1,2);
        hist((bsiShuffled(i,:)),100);
        Lines((bsi(i)),[],'r');
    subplot2(2,3,1,3);
        hist(log2(bsiShuffled(i,:)),100);
        Lines(log2(bsi(i)),[],'r');
        title(bsiZ(i));    
    waitforbuttonpress();
end
 



