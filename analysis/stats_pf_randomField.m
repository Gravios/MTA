function [pfstats,pfshuff] = stats_pf_randomField(Trial,pf,unit)

map = sq(pf.data.rateMap(:,unit==pf.data.clu,:));
pkfr = max(map(:,1));
initRateThreshold = pkfr;

ratemap = pf.plot(unit);

%keyboard
% Maximum firing rate found within the normal place field

% $$$ pfstats.peakFR = pkfr;
% $$$ pfstats.rateThreshold = rateThreshold;
% $$$ pfstats.spatialCoherence = pf.spatialCoherence(unit);
% $$$ pfstats.patchArea = nan([1,1,maxNumPatches]);
% $$$ pfstats.patchCOM  = nan([1,1,maxNumPatches,2]);
% $$$ pfstats.patchPFR = nan([1,1,maxNumPatches]);
% $$$ pfstats.patchMFR = nan([1,1,maxNumPatches]);
% $$$ pfstats.patchCnt = nan([1,1,maxNumPatches]);
% $$$ pfstats.patchRateInd = nan([1,1,maxNumPatches,2,round(prod(pf.adata.binSizes)*0.1)]);
% $$$ pfstats.patchRateMap = nan([1,1,maxNumPatches,round(prod(pf.adata.binSizes)*0.1)]);
rsteps = 20;

patchIndex = [];
patchSize = [];
patchCount = [];
rateThresh = [];

mPatchCount = zeros([size(pf.data.rateMap,3),rsteps]);
mPatchSize = zeros([size(pf.data.rateMap,3),rsteps]);
mRateThresh = zeros([size(pf.data.rateMap,3),rsteps]);

rateThreshSteps = linspace(0,initRateThreshold,rsteps);

for i = 1:size(pf.data.rateMap,3),
    for j = rateThreshSteps,
        ratemap = pf.plot(unit,i);
        [B,~] = bwboundaries(ratemap>j,'noholes');
        if ~isempty(B)            
            pc = numel(B);
            %patchIndex = cat(1,patchIndex,repmat(i,[pc,1]));
            %patchSize = cat(1,patchSize,cellfun(@size,B,repmat({1},[pc,1])));
            %patchCount = cat(1,patchCount,repmat(pc,[pc,1]));
            %rateThresh = cat(1,rateThresh,repmat(j,[pc,1]));
            mPatchCount(i,j==rateThreshSteps) = pc;
            mPatchSize(i,j==rateThreshSteps) =  max(cellfun(@size,B,repmat({1},[pc,1])));

        end
        mRateThresh(i,j==rateThreshSteps) = j;
    end
end

[B,L] = bwboundaries(ratemap>rateThreshold,'noholes');

if ~isempty(B),  
    for i = 1:numel(B)
        [nB{i}(:,1),nB{i}(:,2)] = ind2sub(size(L),find(L==i));
        npix = size(nB{i},1);

        patchArea(1,1,i) = prod(pf.parameters.binDims)*npix;
        patchRateMap(1,1,i,1:npix) = ratemap(sub2ind(pf.adata.binSizes',nB{i}(:,1),nB{i}(:,2)));

        patchRateInd(1,1,i,:,1:npix) = nB{i}';
        pxy = [pf.adata.bins{1}(nB{i}(:,1)),pf.adata.bins{1}(nB{i}(:,2))];
        patchCOM(1,1,i,:)  = [sq(patchRateMap(1,1,i,1:npix))'*pxy(:,1),sq(patchRateMap(1,1,i,1:npix))'*pxy(:,2)]/nansum(patchRateMap(1,1,i,1:npix));
        patchPFR(1,1,i)    = max    (patchRateMap(1,1,i,1:npix),[],4);
        patchMFR(1,1,i)    = nanmean(patchRateMap(1,1,i,1:npix),4);

    end
    
    %keyboard


    [~,pind] = sort(sq(patchArea),'descend');
    if numel(nB)>=maxNumPatches,
        for i = 1:maxNumPatches,
            nNotNan = sum(~isnan(sq(patchRateMap(1,1,pind(i),:))));
            pfstats.patchRateInd(1,1,i,:,1:nNotNan) = patchRateInd(1,1,pind(i),:,:);
            pfstats.patchRateMap(1,1,i,1:nNotNan) = patchRateMap(1,1,pind(i),:);
        end
        pfstats.patchCOM(1,1,:,:)     = patchCOM(1,1,pind(1:maxNumPatches),:);
        pfstats.patchArea(1,1,:)    = patchArea(1,1,pind(1:maxNumPatches));
        pfstats.patchPFR(1,1,:)     = patchPFR(1,1,pind(1:maxNumPatches));
        pfstats.patchMFR(1,1,:)     = patchMFR(1,1,pind(1:maxNumPatches));
    else
        for i = 1:maxNumPatches,
            nNotNan = sum(~isnan(sq(patchRateMap(1,1,1,:))));
            pfstats.patchRateInd(1,1,1,:,1:nNotNan) = patchRateInd(1,1,1,:,:);
            pfstats.patchRateMap(1,1,1,1:nNotNan) = patchRateMap(1,1,1,:);
        end
        pfstats.patchCOM(1,1,1,:,:)   = patchCOM (1,1,1,:);
        pfstats.patchArea(1,1,1)      = patchArea(1,1,1);
        pfstats.patchPFR(1,1,1)       = patchPFR (1,1,1);
        pfstats.patchMFR(1,1,1)       = patchMFR (1,1,1);
        
    end

end


if nargout>1,        
    pfshuff.peakFR = zeros(pf.parameters.numIter,1);
    pfshuff.patchArea = zeros(pf.parameters.numIter,1);
    pfshuff.patchPFR = zeros(pf.parameters.numIter,1);
    pfshuff.patchMFR = zeros(pf.parameters.numIter,1);
    pfshuff.patchCnt = zeros(pf.parameters.numIter,1);;
    %% Patch stats
    bsB = cell(pf.parameters.numIter,1);
    bsPFR = zeros(pf.parameters.numIter,1);
    bsMFR = zeros(pf.parameters.numIter,1);

    for i = 1:pf.parameters.numIter,
        bsB{i} = bwboundaries(reshape(map(:,i),pf.adata.binSizes')'>pfstats.rateThreshold);
        if ~isempty(bsB{i}),
            [~,tmpi] = max(cellfun(@numel,bsB{i}));
            bsB{i} = bsB{i}{tmpi};

            bsMFR(i) = mean(map(sub2ind(pf.adata.binSizes',bsB{i}(:,2),bsB{i}(:,1)),i));
            bsPFR(i) = max(map(sub2ind(pf.adata.binSizes',bsB{i}(:,2),bsB{i}(:,1)),i));
        end
    end

    pfshuff.peakFR = sq(max(map))';
    pfshuff.patchArea = cellfun(@size,bsB,repmat({1},pf.parameters.numIter,1));
    pfshuff.patchPFR = bsPFR;
    pfshuff.patchMFR = bsMFR;
end
