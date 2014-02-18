function pfstats = PlaceFieldStats(Trial,pf,unit)

map = sq(pf.data.rateMap(:,unit==pf.data.clu,:));
pkfr = max(map(:,1));
rateThreshold = prctile(sq(max(map(:,:))),95);

ratemap = reshape(map(:,1),pf.adata.binSizes')';

ratemap = pf.plot(unit);

% Maximum firing rate found within the normal place field
pfstats.peakFR = pkfr;
pfstats.patchArea = 0;
pfstats.patchCOM = zeros(1,2);
pfstats.patchPFR = 0;
pfstats.patchMFR = 0;
pfstats.shuffledPFR = zeros(pf.parameters.numIter,1);
pfstats.shuffledPatchArea = zeros(pf.parameters.numIter,1);
pfstats.shuffledPatchPFR = zeros(pf.parameters.numIter,1);
pfstats.shuffledPatchMFR = zeros(pf.parameters.numIter,1);




B = bwboundaries(ratemap>rateThreshold);

if ~isempty(B),
    [mp,ind] = sort(cellfun(@numel,B)/2,'descend');
    B = B(ind);

    pfstats.patchArea = prod(pf.parameters.binDims)*(numel(B{1})/2);

    % center of mass for the largest patch
    pxy = [pf.adata.bins{1}(B{1}(:,1)),pf.adata.bins{1}(B{1}(:,2))];
    patchRates = ratemap(sub2ind(pf.adata.binSizes',B{1}(:,1),B{1}(:,2)));
    pfstats.patchCOM  = [patchRates'*pxy(:,1),patchRates'*pxy(:,2)]/sum(patchRates);

    % patch peak firing rate
    pfstats.patchPFR = max(patchRates);
    % patch mean firing rate
    pfstats.patchMFR = mean(patchRates);

    %% Patch stats
    bsB = cell(pf.parameters.numIter,1);
    bsPFR = zeros(pf.parameters.numIter,1);
    bsMFR = zeros(pf.parameters.numIter,1);
    parfor i = 1:pf.parameters.numIter,
        bsB{i} = bwboundaries(reshape(map(:,i),pf.adata.binSizes')'>pfstats.peakFR/2);
        if ~isempty(bsB{i}),
            [~,tmpi] = max(cellfun(@numel,bsB{i}));
            bsB{i} = bsB{i}{tmpi};

            bsMFR(i) = mean(map(sub2ind(pf.adata.binSizes',bsB{i}(:,2),bsB{i}(:,1)),i));
            bsPFR(i) = max(map(sub2ind(pf.adata.binSizes',bsB{i}(:,2),bsB{i}(:,1)),i));
        end
    end

    pfstats.shuffledPFR = sq(max(map));
    pfstats.shuffledPatchArea = cellfun(@size,bsB,repmat({1},pf.parameters.numIter,1));
    pfstats.shuffledPatchPFR = bsPFR;
    pfstats.shuffledPatchMFR = bsMFR;

end