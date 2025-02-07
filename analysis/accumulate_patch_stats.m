% SCRIPT accumulate_patch_stats 
% DEPS MjgER2016_figure_BhvPlacefields.m

% ACCUMULATE place fields
% COMPUTE center point of all good patches
% ACCUMULATE independent patches 
% MERGE close patches
% RECOMPUTE patch distance after merger
% MERGE close patches again
% RECOMPUTE patch distance after merger
% SELECT final patches
%% SECOND attempt
% get patch max firing rate across states {'rear','high','low'}
% compute_permuted_patch_ratemap
% COMPUTE the patch bhv rate maps
% COMPUTE inter placecell bhv ratemap correlations
%         only compare units on separate electrodes
cf(@(T) T.load('nq'),Trials);

% ACCUMULATE place fields
t = 20;
[pfstats] = PlaceFieldStats(Trials{t},pfsr{t}{1},20,2,false,true);
for t = 1:30,
    disp(['[INFO] accumulating placefield stats: ' Trials{t}.filebase]);
    for u = units{t},
        [pfstats(end+1,1)] = PlaceFieldStats(Trials{t},pfsr{t}{1},u,2,false,true);
        [pfstats(end  ,2)] = PlaceFieldStats(Trials{t},pfsr{t}{2},u,2,false,true);
        [pfstats(end  ,3)] = PlaceFieldStats(Trials{t},pfsr{t}{3},u,2,false,true);
    end
end
pfstats(1,:) = [];

% $$$ figure();
% $$$ subplot(131);plot([pfstats(:,1).peakFR],[pfstats(:,2).peakFR],'.');xlim([0,30]);ylim([0,30]);line([0,30],[0,30],'Color','k');
% $$$ subplot(132);plot([pfstats(:,1).peakFR],[pfstats(:,3).peakFR],'.');xlim([0,30]);ylim([0,30]);line([0,30],[0,30],'Color','k');
% $$$ subplot(133);plot([pfstats(:,2).peakFR],[pfstats(:,3).peakFR],'.');xlim([0,30]);ylim([0,30]);line([0,30],[0,30],'Color','k');
% $$$ 
% $$$ rpatches = sq([pfstats(:,3).patchMFR]);
% $$$ figure,plot(rpatches(:,1),rpatches(:,2),'.')
% $$$ line([0,30],[0,30],'Color','k');


apfstats = CatStruct(pfstats(1,:),[],2);
for field = fieldnames(pfstats)'
    field = field{1};
    if numel(apfstats.(field)) > 3
        apfstats.(field) = permute(apfstats.(field),[5,1,2,3,4]);
    end
    for uid = 2:size(cluSessionMap,1)
        apfstats.(field) = cat(1,apfstats.(field),cat(2,pfstats(uid,:).(field)));
    end
end

% COMPUTE center point of all good patches
apfstats.patchCenters = nan([size(cluSessionMap,1),3,2,2]);
for uid = 1:size(cluSessionMap,1),
    t = cluSessionMap(uid,1);
    for s = 1:3
        for p = 1:2
            if apfstats.patchPFR(uid,s,p)>1
                rinds = sq(apfstats.patchRateInd(uid,s,p,:,:))';
                rinds = rinds(nniz(rinds),:);
                apfstats.patchCenters(uid,s,p,:) = [mean(pfsr{t}{s}.adata.bins{1}(rinds(:,1))),...
                                    mean(pfsr{t}{s}.adata.bins{2}(rinds(:,2)))];
            end
        end
    end
end


% $$$ % DIAGNOSTIC figure
% $$$ t = 20
% $$$ figure,
% $$$ for u = units{t};
% $$$     uid = find(ismember(cluSessionMap,[t,u],'rows'));
% $$$     mrate = max(apfstats.peakFR(uid,:));
% $$$     clf();
% $$$     for s = 1:3
% $$$         subplot(1,3,s);
% $$$         plot(pfsr{t}{s},u,1,'colorbar',[0,mrate]);
% $$$         hold('on');
% $$$         for p = 1:2
% $$$             if apfstats.patchPFR(uid,s,p)>1
% $$$                 rinds = sq(apfstats.patchRateInd(uid,s,p,:,:))';
% $$$                 rinds = rinds(nniz(rinds),:);
% $$$                 if size(rinds,1)>40
% $$$                     plot(pfsr{t}{s}.adata.bins{1}(rinds(:,1)),...
% $$$                          pfsr{t}{s}.adata.bins{2}(rinds(:,2)),...
% $$$                          'm*');
% $$$                     circle(apfstats.patchCenters(uid,s,p,1),apfstats.patchCenters(uid,s,p,2),150);
% $$$                 end        
% $$$             end
% $$$         end
% $$$         title(num2str(u))
% $$$     end
% $$$     waitforbuttonpress();
% $$$ end


% ACCUMULATE independent patches 
%    IF two patchCenters are closer than 150 mm then count patch with the 
%       final patch -> mean center and radius = 150 + dist(center1,center2)


% COMPUTE center point of all good patches
apfstats.patchCenters = nan([size(cluSessionMap,1),3,2,2]);
for uid = 1:size(cluSessionMap,1),
    t = cluSessionMap(uid,1);
    for s = 1:3
        for p = 1:2
            if apfstats.patchPFR(uid,s,p)>1
                rinds = sq(apfstats.patchRateInd(uid,s,p,:,:))';
                rinds = rinds(nniz(rinds),:);
                apfstats.patchCenters(uid,s,p,:) = [mean(pfsr{t}{s}.adata.bins{1}(rinds(:,1))),...
                                    mean(pfsr{t}{s}.adata.bins{2}(rinds(:,2)))];
            end
        end
    end
end

cdist = nan([size(cluSessionMap,1),3,3,2,2]);
for uid = 1:size(cluSessionMap,1),
    t = cluSessionMap(uid,1);
    for s1 = 1:3
        for s2 = 1:3
            for p1  = 1:2
                for p2 = 1:2
                    cdist(uid,s1,s2,p1,p2) = sqrt(sum(sq(apfstats.patchCenters(uid,s1,p1,:) - ...
                                                         apfstats.patchCenters(uid,s2,p2,:)).^2));
                end
            end
        end
    end
end
tempPatchCenter = nan([size(cluSessionMap,1),3,3,2,2,2]);
tempPatchRadius = nan([size(cluSessionMap,1),3,3,2,2]);
% MERGE close patches
for uid = 1:size(cluSessionMap,1),
    for s1 = 1:3
        for s2 = 1:3
            for p1 = 1:2
                for p2 = 1:2                
                    if cdist(uid,s1,s2,p1,p2) < 150           ...
                            && apfstats.patchPFR(uid,s1,p1)>1 ...
                            && apfstats.patchPFR(uid,s2,p2)>1 ...
                            && apfstats.patchPFR(uid,s1,p2)>1 ...
                            && apfstats.patchPFR(uid,s2,p1)>1 ...
                            && sum(nniz(apfstats.patchRateMap(uid,s1,p1)))>40 ...
                            && sum(nniz(apfstats.patchRateMap(uid,s1,p2)))>40 ...
                            && sum(nniz(apfstats.patchRateMap(uid,s2,p1)))>40 ...
                            && sum(nniz(apfstats.patchRateMap(uid,s2,p2)))>40 
                        % 18 comparisons
                        % don't look at same state patch centers 
                        tempPatchCenter(uid,s1,s2,p1,p2,:) = (apfstats.patchCenters(uid,s1,p1,:)+apfstats.patchCenters(uid,s2,p2,:))./2;
                        tempPatchRadius(uid,s1,s2,p1,p2) = 150 + cdist(uid,s1,s2,p1,p2);
                        cdist(uid,s1,s2,p1,p2) = 0;
                        apfstats.patchCenters(uid,s1,p1,:) = tempPatchCenter(uid,s1,s2,p1,p2,:);
                        apfstats.patchCenters(uid,s2,p2,:) = tempPatchCenter(uid,s1,s2,p1,p2,:);
                    end
                end
            end
        end
    end
end

% RECOMPUTE patch distance after merger
cdist = nan([size(cluSessionMap,1),3,3,2,2]);
for uid = 1:size(cluSessionMap,1),
    t = cluSessionMap(uid,1);
    for s1 = 1:3
        for s2 = 1:3
            for p1  = 1:2
                for p2 = 1:2
                    cdist(uid,s1,s2,p1,p2) = sqrt(sum(sq(apfstats.patchCenters(uid,s1,p1,:) - ...
                                                         apfstats.patchCenters(uid,s2,p2,:)).^2));
                end
            end
        end
    end
end
% MERGE close patches again
tempPatchCenter = nan([size(cluSessionMap,1),3,3,2,2,2]);
tempPatchRadius = nan([size(cluSessionMap,1),3,3,2,2]);
for uid = 1:size(cluSessionMap,1),
    for s1 = 1:3
        for s2 = 1:3
            for p1 = 1:2
                for p2 = 1:2                
                    if cdist(uid,s1,s2,p1,p2) < 150           
                        % 18 comparisons
                        % don't look at same state patch centers 
                        tempPatchCenter(uid,s1,s2,p1,p2,:) = (apfstats.patchCenters(uid,s1,p1,:)+apfstats.patchCenters(uid,s2,p2,:))./2;
                        tempPatchRadius(uid,s1,s2,p1,p2) = 150 + cdist(uid,s1,s2,p1,p2);
                        cdist(uid,s1,s2,p1,p2) = 0;
                        apfstats.patchCenters(uid,s1,p1,:) = tempPatchCenter(uid,s1,s2,p1,p2,:);
                        apfstats.patchCenters(uid,s2,p2,:) = tempPatchCenter(uid,s1,s2,p1,p2,:);
                    end
                end
            end
        end
    end
end

% RECOMPUTE patch distance after merger
cdist = nan([size(cluSessionMap,1),3,3,2,2]);
for uid = 1:size(cluSessionMap,1),
    t = cluSessionMap(uid,1);
    for s1 = 1:3
        for s2 = 1:3
            for p1  = 1:2
                for p2 = 1:2
                    cdist(uid,s1,s2,p1,p2) = sqrt(sum(sq(apfstats.patchCenters(uid,s1,p1,:) - ...
                                                         apfstats.patchCenters(uid,s2,p2,:)).^2));
                end
            end
        end
    end
end
% SELECT final patches
tempPatchCenter = nan([size(cluSessionMap,1),3,3,2,2,2]);
tempPatchRadius = nan([size(cluSessionMap,1),3,3,2,2]);
for uid = 1:size(cluSessionMap,1),
    for s1 = 1:3
        for s2 = 1:3
            for p1 = 1:2
                for p2 = 1:2                
                    if cdist(uid,s1,s2,p1,p2) < 150 
                        % 18 comparisons
                        % don't look at same state patch centers 
                        tempPatchCenter(uid,s1,s2,p1,p2,:) = (apfstats.patchCenters(uid,s1,p1,:)+apfstats.patchCenters(uid,s2,p2,:))./2;
                        tempPatchRadius(uid,s1,s2,p1,p2) = 150 + cdist(uid,s1,s2,p1,p2);
                    else
                        tempPatchCenter(uid,s1,s2,p1,1,:) = apfstats.patchCenters(uid,s1,p1,:);
                        tempPatchCenter(uid,s1,s2,p1,2,:) = apfstats.patchCenters(uid,s2,p2,:);
                        tempPatchRadius(uid,s1,s2,p1,p2) = 150 + cdist(uid,s1,s2,p1,p2);
                    end
                end
            end
        end
    end
end

apfstats.finalPatchCenters = nan([size(cluSessionMap,1),10,2]);
for uid = 1:size(cluSessionMap,1)
    t = cluSessionMap(uid,1);
    u = cluSessionMap(uid,2);
    [ patches, patchInd ] = unique(reshape(sq(tempPatchCenter(uid,:,:,:,:)),[],2),'rows');
    patches = patches(nniz(patches),:);
    for p1 = 1:size(patches,1)
        for p2 = p1+1:size(patches,1)
            patchDist = sqrt(sum((patches(p1,:)-patches(p2,:)).^2));
            if patchDist < 150
                patches(p1,:) = (patches(p1,:)+patches(p2,:))/2;
                patches(p2,:) = nan([1,2]);
            end
        end
    end
    
    for p =  1:size(patches,1)
        if  isnan(patches(p,1))
            continue
        end
        [~,mapindx] = NearestNeighbour(pfsr{t}{s}.adata.bins{1},patches(p,1));
        [~,mapindy] = NearestNeighbour(pfsr{t}{s}.adata.bins{2},patches(p,2));        
        rmap = cf(@(pp) plot(pp,u), pfsr{t});
        if ~any(cellfun(@(rr) rr(mapindx,mapindy), rmap)>1)
            patches(p,:) = nan;
        end
    end
    patches(~nniz(patches),:) = [];
    apfstats.finalPatchCenters(uid,1:size(patches,1),:) = patches;
end


% $$$ 
% $$$ t = 4
% $$$ figure,
% $$$ for u = units{t};
% $$$     uid = find(ismember(cluSessionMap,[t,u],'rows'));
% $$$     mrate = max(apfstats.peakFR(uid,:));
% $$$     clf();
% $$$     for s = 1:3
% $$$         subplot(1,3,s);
% $$$         plot(pfsr{t}{s},u,1,'colorbar',[0,mrate]);
% $$$         hold('on');
% $$$         for p = 1:2
% $$$             if apfstats.patchPFR(uid,s,p)>1
% $$$                 rinds = sq(apfstats.patchRateInd(uid,s,p,:,:))';
% $$$                 rinds = rinds(nniz(rinds),:);
% $$$                 if size(rinds,1)>40
% $$$                     plot(pfsr{t}{s}.adata.bins{1}(rinds(:,1)),...
% $$$                          pfsr{t}{s}.adata.bins{2}(rinds(:,2)),...
% $$$                          'm*');                    
% $$$                 end        
% $$$             end
% $$$         end
% $$$         for p = 1:sum(nniz(apfstats.finalPatchCenters(uid,:,1)'))
% $$$             circle(apfstats.finalPatchCenters(uid,p,1),apfstats.finalPatchCenters(uid,p,2),175);
% $$$         end
% $$$         title(num2str(u))
% $$$     end
% $$$     waitforbuttonpress();
% $$$ end

%% SECOND attempt
numUnits = size(cluSessionMap,1);
numStates = numel(pfsr{1});
maxNumPatches = 2;
rthresh = 1.5;
dthresh = 100;
athresh = 9600;%16000
numPatches =numel(pfsr{1})*maxNumPatches;
patchArea = reshape(apfstats.patchArea,[],numPatches);
patchPFR  = reshape(apfstats.patchPFR,[],numPatches);
patchInds = reshape(apfstats.patchRateInd,[],numPatches,2,550);
patchCntr = nan([numUnits,numPatches,2]);
for uid = 1:numUnits
    for p = 1:numPatches
        if patchArea(uid,p) < athresh  |  patchPFR(uid,p) < rthresh
            patchInds(uid,p,:,:) = nan;
            patchArea(uid,p) = nan;
            patchPFR(uid,p) = nan;
        end
    end
end

for uid = 1:numUnits
    t = cluSessionMap(uid,1);
    for p = 1:numPatches
        if ~isnan( patchArea(uid,p) )
            rinds = sq(patchInds(uid,p,:,:))';
            rinds = rinds(nniz(rinds),:);
            patchCntr(uid,p,:) = [mean(pfsr{t}{1}.adata.bins{1}(rinds(:,1))),...
                                  mean(pfsr{t}{1}.adata.bins{2}(rinds(:,2)))];
        end
    end
end

patchDist = nan([numUnits,numPatches,numPatches]);
for uid = 1:numUnits
    for p1 = 1:numPatches    
        for p2 = 1:numPatches
            patchDist(uid,p1,p2) = sqrt(sum((patchCntr(uid,p1,:)-patchCntr(uid,p2,:)).^2));
        end
    end
end



patchCntrF = nan([size(patchCntr)]);
patchCntrI = nan([size(patchCntr)]);
for uid = 1:numUnits
    [cnt,sind] = sort(sum(sq(patchDist(uid,:,:))<dthresh  &  sq(patchDist(uid,:,:))~=0),'descend');
    for p = 1:numel(cnt),
        pweights = permute(repmat(sq(patchPFR(uid,sq(patchDist(uid,sind(p),:))<dthresh)./sum(patchPFR(uid,sq(patchDist(uid,sind(p),:))<dthresh))),...
                          [2,1]),...
                           [3,2,1]);
        patchCntrI(uid,p,:) = sum(patchCntr(uid,sq(patchDist(uid,sind(p),:))<dthresh,:).*pweights,2);
    end
    patchCntrUnique = unique(sq(patchCntrI(uid,:,:)),'rows');
    patchCntrUnique(~nniz(patchCntrUnique),:) = [];
    patchCntrF(uid,1:size(patchCntrUnique,1),:) = patchCntrUnique;
end


patchCntr = patchCntrF;
patchCntrF = nan([size(patchCntr)]);
patchCntrI = nan([size(patchCntr)]);
patchDist = nan([numUnits,numPatches,numPatches]);
for uid = 1:numUnits
    for p1 = 1:numPatches    
        for p2 = 1:numPatches
            patchDist(uid,p1,p2) = sqrt(sum((patchCntr(uid,p1,:)-patchCntr(uid,p2,:)).^2));
        end
    end
end
for uid = 1:numUnits
    [cnt,sind] = sort(sum(sq(patchDist(uid,:,:))<dthresh  &  sq(patchDist(uid,:,:))~=0),'descend');
    for p = 1:numel(cnt),
        patchCntrI(uid,p,:) = mean(patchCntr(uid,sq(patchDist(uid,sind(p),:))<dthresh,:),2);
    end
    patchCntrUnique = unique(sq(patchCntrI(uid,:,:)),'rows');
    patchCntrUnique(~nniz(patchCntrUnique),:) = [];
    patchCntrF(uid,1:size(patchCntrUnique,1),:) = patchCntrUnique;
end

patchDistF = nan([numUnits,numPatches,numPatches]);
for uid = 1:numUnits
    for p1 = 1:numPatches    
        for p2 = 1:numPatches
            patchDistF(uid,p1,p2) = sqrt(sum((patchCntrF(uid,p1,:)-patchCntrF(uid,p2,:)).^2));
        end
    end
end


% $$$ t = 20;
% $$$ figure,
% $$$ for u = units{t};
% $$$     uid = find(ismember(cluSessionMap,[t,u],'rows'));
% $$$     mrate = max(apfstats.peakFR(uid,:));
% $$$     clf();
% $$$     for s = 1:3
% $$$         subplot(1,3,s);
% $$$         plot(pfsr{t}{s},u,1,'colorbar',[0,mrate]);
% $$$         hold('on');
% $$$         for p = 1:2
% $$$             if apfstats.patchPFR(uid,s,p)>rthresh
% $$$                 rinds = sq(apfstats.patchRateInd(uid,s,p,:,:))';
% $$$                 rinds = rinds(nniz(rinds),:);
% $$$                 if size(rinds,1)>athresh/400
% $$$                     plot(pfsr{t}{s}.adata.bins{1}(rinds(:,1)),...
% $$$                          pfsr{t}{s}.adata.bins{2}(rinds(:,2)),...
% $$$                          'm*');                    
% $$$                 end        
% $$$             end
% $$$         end
% $$$         for p = 1:sum(nniz(patchCntrF(uid,:,1)'))
% $$$             circle(patchCntrF(uid,p,1),patchCntrF(uid,p,2),150);
% $$$         end
% $$$         title(num2str(u))
% $$$     end
% $$$      waitforbuttonpress();
% $$$ end




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




uid = find(ismember(cluSessionMap,[20,20],'rows'));            
% $$$ reshape(sq(tempPatchCenter(uid,:,:,:,:)),[],2),'rows')
% $$$ 
% $$$ figure,imagesc(sq(cdist(uid,:,:,1)))
% $$$ figure,imagesc(sq(cdist(uid,:,:,1)))            

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


%figure,hist(nonzeros(patchArea)/100,50)
figure,hist(sqrt(nonzeros(patchArea)),50)

%% Permutations 


% center1 center2
% area around center1 as state

overwrite = false;
threshDist = 150;
sampleRate = 16;
pfsArgs = struct('states',           'theta-groom-sit',  ...
                 'binDims',          [0.1,0.1],          ...
                 'SmoothingWeights', [1.8,1.8],          ...
                 'numIter',          200,                ...                
                 'boundaryLimits',   [-2,0.8;-0.8,2],    ...
                 'mask',             validDims,          ...
                 'halfsample',       false);

% $$$ unit = 34;
%% compute_permuted_patch_ratemap
%uids = find(ismember(cluSessionMap(:,1),tid));

rmapA = zeros([784,numUnits,6,6]);
rmapB = zeros([784,numUnits,6,6]);
rmapZscr = zeros([784,numUnits,6,6]);
rmapCorr = zeros([200,numUnits,6,6]);
tid = 0;
for uid = 1:numUnits
    if tid ~= cluSessionMap(uid,1) | uid==1
        tid = cluSessionMap(uid,1);
        Trial = Trials{tid};
        xyz = preproc_xyz(Trial,'trb',sampleRate);
        fet = fet_HB_pitchB(Trial,sampleRate);
        spk  = create(copy(Trial.spk),Trial,sampleRate,pfsArgs.states,units{tid},'');
    end
        
    unit = cluSessionMap(uid,2);
    numPatches = sum(~isnan(patchCntrF(uid,:,1)));
    tic
    if numPatches > 1
        for p1 = 1:numPatches-1
            for p2 = p1+1:numPatches
                tag = ['pfsPerm_',num2str([unit,p1,p2],'%d_%d-%d')];
% $$$                 [rmapDiff,ratemapA,ratemapB] = ...
                [rmapZscr(:,uid,p1,p2),                                         ...
                 rmapA(:,uid,p1,p2),rmapB(:,uid,p1,p2),                         ...
                 rmapCorr(:,uid,p1,p2)] =                                       ...
                    compute_permuted_patch_ratemap(Trial,                       ...
                                                   unit,                        ...
                                                   [],                          ... fetset
                                                   sampleRate,                  ...
                                                   patchCntrF(uid,[p1,p2],:),   ...
                                                   'hcom',                      ... marker
                                                   pfsArgs,                     ... 
                                                   [],                          ... threshRate
                                                   threshDist,                  ...
                                                   xyz,                         ...
                                                   fet,                         ...
                                                   spk,                         ...
                                                   tag,                         ...
                                                   overwrite);

            end
        end
    end
    toc
end

% Hmm this is where I need to get the code for the intra patch angle corr... ah copy the p1 p2 part


rmapAngle = zeros([numUnits,6,6]);
rmapDist = zeros([numUnits,6,6]);;
for uid = 1:numUnits
    unit = cluSessionMap(uid,2);
    numPatches = sum(~isnan(patchCntrF(uid,:,1)));
    if numPatches > 1
        for p1 = 1:numPatches-1
            for p2 = p1+1:numPatches
                    a = acos(dot(patchPFRF(uid,:,p1),patchPFRF(uid,:,p2)) ./ ...
                             (sqrt(sum(patchPFRF(uid,:,p1).^2)).*sqrt(sum(patchPFRF(uid,:,p2).^2))));
                    rmapAngle(uid,p1,p2) = a;
                    rmapDist(uid,p1,p2) = patchDistF(uid,p1,p2);
            end
        end
    end
end

rangle = abs(reshape(rmapAngle(unitSubset,:,:),[],1));
pd = reshape(patchDistF(unitSubset,:,:),[],1)/10;
% $$$ nind = nniz(pd) & nniz(rangle) & nniz(rcorr) & pd>10;
% $$$ figure,plot(pd(nind),rcorr(nind),'.');
% $$$ figure,plot(pd(nind),rangle(nind),'.');
% $$$ figure,plot(rcorr(nind),rangle(nind),'.');


% now how do I compare this to the rmapCorr?
                


%%%<<< COMPUTE the patch bhv rate maps -------------------------------------------------------------
rmapP = zeros([784,numUnits,6]);
tid = 0;
for uid = 1:numUnits
    if tid ~= cluSessionMap(uid,1)
        tid = cluSessionMap(uid,1);
        Trial = Trials{tid};
        xyz = preproc_xyz(Trial,'trb',sampleRate);
        fet = fet_HB_pitchB(Trial,sampleRate);
        spk  = create(copy(Trial.spk),Trial,sampleRate,pfsArgs.states,units{tid},'');
    end
    unit = cluSessionMap(uid,2);
    numPatches = sum(~isnan(patchCntrF(uid,:,1)));
    tic
    for p = 1:numPatches
        [rmapP(:,uid,p)] =                                         ...
                compute_patch_ratemap(Trial,                       ...
                                      unit,                        ...
                                      [],                          ... fetset
                                      sampleRate,                  ...
                                      patchCntrF(uid,p,:),         ...
                                      'hcom',                      ... marker
                                      pfsArgs,                     ... 
                                      [],                          ... threshRate
                                      threshDist,                  ...
                                      xyz,                         ...
                                      fet,                         ...
                                      spk,                         ...
                                      tag,                         ...
                                      overwrite);
    end
    toc
end
%%%>>>


% $$$ figure,
% $$$ for p = 1:6,
% $$$ subplot(1,6,p);
% $$$ imagesc(reshape(rmapP(:,ismember(cluSessionMap,[20,31],'rows'),p).*validDims,[28,28])');axis('xy');colormap('jet');colorbar();
% $$$ end


% $$$ bins = bfs{1}.adata.bins
% $$$ figure();
% $$$ stateLabels = {'Rear','High','Low'};
% $$$ tid = 20;
% $$$ u = 34;
% $$$ p1 = 1; p2 = 2;
% $$$ %for u = units{tid};
% $$$     uid = find(ismember(cluSessionMap,[tid,u],'rows'));
% $$$     mrate = max(apfstats.peakFR(uid,:));
% $$$     clf();
% $$$     for s = 1:3
% $$$         subplot2(2,3,1,s);
% $$$         plot(pfsr{tid}{s},u,1,'colorbar',[0,mrate],'colorMap',@parula);
% $$$         hold('on');
% $$$         for p = 1:2
% $$$             if apfstats.patchPFR(uid,s,p)>rthresh
% $$$                 rinds = sq(apfstats.patchRateInd(uid,s,p,:,:))';
% $$$                 rinds = rinds(nniz(rinds),:);
% $$$                 if size(rinds,1)>athresh/400
% $$$                     plot(pfsr{tid}{s}.adata.bins{1}(rinds(:,1)),...
% $$$                          pfsr{tid}{s}.adata.bins{2}(rinds(:,2)),...
% $$$                          'm.');                    
% $$$                 end        
% $$$             end
% $$$         end
% $$$         for p = 1:sum(nniz(patchCntrF(uid,:,1)'))
% $$$             circle(patchCntrF(uid,p,1),patchCntrF(uid,p,2),150);
% $$$             text(patchCntrF(uid,p,1),patchCntrF(uid,p,2),num2str(p));
% $$$         end
% $$$         title({['Unit: ',num2str(u)],[stateLabels{s}, ' State Ratemap']});
% $$$     end
% $$$     subplot2(2,3,2,1);
% $$$         imagesc(bins{:},reshape(rmapA(:,uid,p1,p2).*validDims,cellfun(@numel,bins))');
% $$$         title(['Patch ',num2str(p1),' bhv ratemap']);
% $$$         axis('xy');
% $$$         cax = colorbar();
% $$$         ylabel(cax,'Hz');                
% $$$         colormap(gca(),'jet');
% $$$         xlabel('Head-Body pitch (rad)');
% $$$         ylabel('Body pitch (rad)');    
% $$$     subplot2(2,3,2,2);
% $$$         imagesc(bins{:},reshape(rmapB(:,uid,p1,p2).*validDims,cellfun(@numel,bins))');
% $$$         title(['Patch ',num2str(p2),' bhv ratemap']);
% $$$         xlabel('Head-Body pitch (rad)');
% $$$         ylabel('Body pitch (rad)');
% $$$         axis('xy');
% $$$         cax = colorbar();
% $$$         ylabel(cax,'Hz');                
% $$$         colormap(gca(),'jet');
% $$$     subplot2(2,3,2,3);
% $$$         imagesc(bins{:},reshape(rmapZscr(:,uid,p1,p2).*validDims,cellfun(@numel,bins))');
% $$$         axis('xy');
% $$$         cax = colorbar();
% $$$         ylabel(cax,'z-score');
% $$$         colormap(gca(),'jet');    
% $$$         title(num2str([p1,p2],'%d <-> %d'));
% $$$         xlabel('Head-Body pitch (rad)');
% $$$         ylabel('Body pitch (rad)');
% $$$         title('Permuted Patch Ratemap Z-Score');
% $$$         
% $$$     
% $$$ % $$$     waitforbuttonpress();
% $$$ % $$$ end



%%%<<< COMPUTE inter placecell bhv ratemap correlations
% only compare units on separate electrodes
rmapIPCorr = [];
rmapIPDist = [];
rmapIPBcnt = [];
for tid = numel(Trials),
    for u1 = units{tid}
        uid1 = find(ismember(cluSessionMap,[tid,u1],'rows'));
        for u2 = units{tid}
            if u1 == u2 ...
               || (Trials{tid}.nq.ElNum(u1)==Trials{tid}.nq.ElNum(u2) ...
                   && abs(Trials{tid}.nq.maxAmpChan(u1)-Trials{tid}.nq.maxAmpChan(u2))<=2)
                continue
            end
            uid2 = find(ismember(cluSessionMap,[tid,u2],'rows'));
            for gptch1 = find(~isnan(patchCntrF(uid1,:,2)))
                for gptch2 = find(~isnan(patchCntrF(uid2,:,2)))
                    rmapPC = [rmapP(validDims,uid1,gptch1),rmapP(validDims,uid2,gptch2)];
                    if any(sum(isnan(rmapPC)) > 50)
                        continue
                    end
                    c = corr(rmapPC(~isnan(rmapPC(:,1))&~isnan(rmapPC(:,2)),:),'type','Spearman');
                    rmapIPCorr(end+1) = c(2);
                    rmapIPDist(end+1) = sqrt(sum((patchCntrF(uid1,gptch1,:)-patchCntrF(uid2,gptch2,:)).^2,3));
                    rmapIPBcnt(end+1) = sum(double(~isnan(rmapPC(:,1))&~isnan(rmapPC(:,2))));
                end
            end
        end
    end
end

rmapPCA = nan([784,1]);
rmapPCB = nan([784,1]);
rmapPCA(validDims) = rmapPC(:,1);
rmapPCB(validDims) = rmapPC(:,2);
figure,
subplot(121);
pcolor(reshape(rmapPCA,[28,28])');axis('xy');
subplot(122);
pcolor(reshape(rmapPCB,[28,28])');axis('xy');

% $$$ figure,
% $$$ uid = find(ismember(cluSessionMap,[20,20],'rows'));
% $$$ imagesc(reshape(rmapP(:,uid,1),[28,28])');
% $$$ axis('xy');
% $$$ 
figure,
plot(rmapIPDist/10,rmapIPCorr,'.')
xlabel('Inter Patch Distance (cm)')
ylabel('Bhv Ratemap Correlation');

figure,hist(rmapIPCorr,20)

% $$$ 
% $$$ patchCntF = sum(~isnan(patchCntrF(:,:,1)),2);
% $$$ out = histcounts(patchCntF,0.5:5.5);
% $$$ figure,pie(out)
%%%>>>
% $$$ 
% $$$ % rmapIPCorr should mirror rmapIPAngle
% $$$ rmapIPAngle = []; 
% $$$ for tid = numel(Trials),
% $$$     for u1 = units{tid}
% $$$         uid1 = find(ismember(cluSessionMap,[tid,u1],'rows'));
% $$$         for u2 = units{tid}
% $$$             if u1 == u2 ...
% $$$                || (Trials{tid}.nq.ElNum(u1)==Trials{tid}.nq.ElNum(u2) ...
% $$$                    && abs(Trials{tid}.nq.maxAmpChan(u1)-Trials{tid}.nq.maxAmpChan(u2))<=2)
% $$$                 continue
% $$$             end
% $$$             uid2 = find(ismember(cluSessionMap,[tid,u2],'rows'));
% $$$             for gptch1 = find(~isnan(patchCntrF(uid1,:,2)))
% $$$                 for gptch2 = find(~isnan(patchCntrF(uid2,:,2)))
% $$$                     % have to keep this to make comparison to bhvfield corr fair 
% $$$                     rmapPC = [rmapP(validDims,uid1,gptch1),rmapP(validDims,uid2,gptch2)];
% $$$                     if any(sum(isnan(rmapPC)) > 50)
% $$$                         continue
% $$$                     end
% $$$                     a = acos(dot(patchPFRF(uid1,:,gptch1),patchPFRF(uid2,:,gptch2)) ./ ...
% $$$                              (sqrt(sum(patchPFRF(uid1,:,gptch1).^2)).*sqrt(sum(patchPFRF(uid2,:,gptch2).^2))));
% $$$                     rmapIPAngle(end+1) = a;
% $$$                 end
% $$$             end
% $$$         end
% $$$     end
% $$$ end

% $$$ figure();
% $$$ subplot(131);
% $$$     plot(rmapIPDist,rmapIPAngle,'.');
% $$$     xlabel('between place cell inter patch distance (cm)');
% $$$     ylabel('angle between state firing rate vectors (rad)');
% $$$ subplot(132);
% $$$     plot(rmapIPDist,rmapIPCorr,'.');
% $$$     xlabel('between place cell inter patch distance (cm)');
% $$$     ylabel('correlation between bhv firing rate maps (A.U.)');
% $$$ subplot(133);
% $$$     plot(rmapIPCorr,rmapIPAngle,'.');
% $$$     xlabel('correlation between bhv firing rate maps (A.U.)');
% $$$     ylabel('angle between state firing rate vectors (rad)');    
% $$$ % end supfig



        
% $$$ stateLabels = {'Rear','High','Low'};
% $$$ tid = 20;
% $$$ u = 25;
% $$$ p1 = 1; p2 = 2;
% $$$ %for u = units{tid};
% $$$     uid = find(ismember(cluSessionMap,[tid,u],'rows'));
% $$$     mrate = max(apfstats.peakFR(uid,:));
% $$$     clf();
% $$$     for s = 1:3
% $$$         subplot2(2,3,1,s);
% $$$         plot(pfsr{tid}{s},u,1,'colorbar',[0,mrate],'colorMap',@parula);
% $$$         hold('on');
% $$$         for p = 1:2
% $$$             if apfstats.patchPFR(uid,s,p)>rthresh
% $$$                 rinds = sq(apfstats.patchRateInd(uid,s,p,:,:))';
% $$$                 rinds = rinds(nniz(rinds),:);
% $$$                 if size(rinds,1)>athresh/400
% $$$                     plot(pfsr{tid}{s}.adata.bins{1}(rinds(:,1)),...
% $$$                          pfsr{tid}{s}.adata.bins{2}(rinds(:,2)),...
% $$$                          'm.');                    
% $$$                 end        
% $$$             end
% $$$         end
% $$$         for p = 1:sum(nniz(patchCntrF(uid,:,1)'))
% $$$             circle(patchCntrF(uid,p,1),patchCntrF(uid,p,2),150);
% $$$             text(patchCntrF(uid,p,1),patchCntrF(uid,p,2),num2str(p));
% $$$         end
% $$$         title({['Unit: ',num2str(u)],[stateLabels{s}, ' State Ratemap']});
% $$$     end
% $$$         
% $$$ set(gcf(),'PaperOrientation','landscape');
% $$$ set(gcf(),'PaperType','A4');
% $$$         
% $$$ nrA = ratemapA(validDims,uid,p1,p2);
% $$$ nrB = ratemapB(validDims,uid,p1,p2);
% $$$ nind = nniz(nrA) & nniz(nrB);
% $$$ nrA = nrA(nind);
% $$$ nrB = nrB(nind);
% $$$ 
% $$$ corr(nrA,nrB,'type','Spearman')
% $$$ 
% $$$ uid = find(ismember(cluSessionMap,[20,44],'rows'));
% $$$ patchCenter = patchCntrF(uid,[1,2],:);

% $$$ % REPORT intra cell patch distance vs correlation 
% $$$ figure
% $$$ hold('on');
% $$$ ph = patch([10.1,10.1,30.0,30.0],[-1,1,1,-1],[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',0.3);
% $$$ pd = reshape(patchDistF(unitSubset,:,:),[],1)/10;
% $$$ rcorr = reshape(rmapCorr(1,unitSubset,:,:),[],1);
% $$$ nind = nniz(pd) & nniz(rcorr) & pd>10;
% $$$ plot(pd(nind),rcorr(nind),'.');
% $$$ box('on');
% $$$ xlabel('Patch Distance (cm)')
% $$$ ylabel('Bhv Correlation');
% $$$ 
% $$$ 
% $$$ 
% $$$ 

% $$$ polyfit(pd(nind),rcorr(nind),'.');
% $$$ 
% $$$ ind = unitSubset;
% $$$ ind = ':';
% $$$ figure,
% $$$ subplot(131);
% $$$ plot(reshape(patchPFRF(ind,2,:),[],1),reshape(patchPFRF(ind,1,:),[],1),'.');xlim([0,25]);ylim([0,25]);
% $$$ line([0,25],[0,25],'Color','k');
% $$$ line([0,12.5],[0,25],'Color','r');
% $$$ line([0,25],[0,12.5],'Color','r');
% $$$ xlabel('Peak Patch Rate Hz (high)');
% $$$ ylabel('Peak Patch Rate Hz (rear)');
% $$$ subplot(132);
% $$$ plot(reshape(patchMFRF(ind,3,:),[],1),reshape(patchMFRF(ind,1,:),[],1),'.');xlim([0,25]);ylim([0,25]);
% $$$ line([0,25],[0,25],'Color','k');
% $$$ line([0,12.5],[0,25],'Color','r');
% $$$ line([0,25],[0,12.5],'Color','r');
% $$$ xlabel('Peak Patch Rate Hz (low)');
% $$$ ylabel('Peak Patch Rate Hz (rear)');
% $$$ subplot(133);
% $$$ plot(reshape(patchPFRF(ind,3,:),[],1),reshape(patchPFRF(ind,2,:),[],1),'.');xlim([0,25]);ylim([0,25]);
% $$$ line([0,25],[0,25],'Color','k');
% $$$ line([0,12.5],[0,25],'Color','r');
% $$$ line([0,25],[0,12.5],'Color','r');
% $$$ xlabel('Peak Patch Rate Hz (low)');
% $$$ ylabel('Peak Patch Rate Hz (high)');

% $$$ 
% $$$ ind = unitSubset;
% $$$ %ind = ':';
% $$$ figure,
% $$$ subplot(131);
% $$$ plot(log10(reshape(patchPFRF(ind,2,:),[],1)),log10(reshape(patchPFRF(ind,1,:),[],1)),'.');
% $$$ xlim([-1.5,1.5]);ylim([-1.5,1.5]);
% $$$ xlabel('Peak Patch Rate Hz (high)');
% $$$ ylabel('Peak Patch Rate Hz (rear)');
% $$$ grid('on');
% $$$ subplot(132);
% $$$ plot(log10(reshape(patchPFRF(ind,3,:),[],1)),log10(reshape(patchPFRF(ind,1,:),[],1)),'.');
% $$$ xlim([-1.5,1.5]);ylim([-1.5,1.5]);
% $$$ xlabel('Peak Patch Rate Hz (low)');
% $$$ ylabel('Peak Patch Rate Hz (rear)');
% $$$ grid('on');
% $$$ subplot(133);
% $$$ plot(log10(reshape(patchPFRF(ind,3,:),[],1)),log10(reshape(patchPFRF(ind,2,:),[],1)),'.');
% $$$ xlim([-1.5,1.5]);ylim([-1.5,1.5]);
% $$$ xlabel('Peak Patch Rate Hz (low)');
% $$$ ylabel('Peak Patch Rate Hz (high)');
% $$$ grid('on');


% $$$ figure();
% $$$ subplot(131);
% $$$ plot(log2(reshape(patchPFRF(ind,1,:),[],1)./reshape(patchPFRF(ind,2,:),[],1)),...
% $$$      log2(reshape(patchPFRF(ind,1,:),[],1)./reshape(patchPFRF(ind,3,:),[],1)),...
% $$$      '.');
% $$$ subplot(132);
% $$$ plot(log2(reshape(patchPFRF(ind,2,:),[],1)./reshape(patchPFRF(ind,1,:),[],1)),...
% $$$      log2(reshape(patchPFRF(ind,2,:),[],1)./reshape(patchPFRF(ind,3,:),[],1)),...
% $$$      '.');
% $$$ subplot(133);
% $$$ plot(log2(reshape(patchPFRF(ind,3,:),[],1)./reshape(patchPFRF(ind,1,:),[],1)),...
% $$$      log2(reshape(patchPFRF(ind,3,:),[],1)./reshape(patchPFRF(ind,2,:),[],1)),...
% $$$      '.');

% SHOULD i compute the same patch PFR ration for shuffeled state id?
% PRETTY sure I've done this in the past ...



% $$$ 
% $$$ figure();
% $$$ subplot(131);
% $$$ violin((reshape(patchPFRF(ind,1,:),[],1)-reshape(patchPFRF(ind,2,:),[],1))./(reshape(patchPFRF(ind,1,:),[],1)+reshape(patchPFRF(ind,2,:),[],1)));
% $$$ subplot(132);
% $$$ violin((reshape(patchPFRF(ind,1,:),[],1)-reshape(patchPFRF(ind,3,:),[],1))./(reshape(patchPFRF(ind,1,:),[],1)+reshape(patchPFRF(ind,3,:),[],1)));
% $$$ subplot(133);
% $$$ violin((reshape(patchPFRF(ind,3,:),[],1)-reshape(patchPFRF(ind,2,:),[],1))./(reshape(patchPFRF(ind,3,:),[],1)+reshape(patchPFRF(ind,2,:),[],1)));
% $$$ 
% $$$ figure
% $$$ histogram(patchCNTF(ind),0.5:6.5);
% $$$  