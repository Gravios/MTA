function [pfstats,pfshuff] = PlaceFieldStats(Trial,pf,unit)

map = sq(pf.data.rateMap(:,unit==pf.data.clu,:));
pkfr = max(map(:,1));
rateThreshold = prctile(sq(max(map(:,:))),95);


ratemap = pf.plot(unit);

% Maximum firing rate found within the normal place field
pfstats.peakFR = pkfr;
pfstats.rateThreshold = rateThreshold;
pfstats.patchArea = 0;
pfstats.patchCOM = zeros(1,2);
pfstats.patchPFR = 0;
pfstats.patchMFR = 0;
pfstats.patchCnt = 0;




B = bwboundaries(ratemap>rateThreshold);

if ~isempty(B),
    [mp,ind] = sort(cellfun(@numel,B)/2,'descend');
    B = B(ind);

    pfstats.patchArea = prod(pf.parameters.binDims)*(numel(B{1})/2);
    % Center of mass for the largest patch
    pxy = [pf.adata.bins{1}(B{1}(:,1)),pf.adata.bins{1}(B{1}(:,2))];
    patchRates = ratemap(sub2ind(pf.adata.binSizes',B{1}(:,1),B{1}(:,2)));
    pfstats.patchCOM  = ([patchRates'*pxy(:,1),patchRates'*pxy(:,2)]/sum(patchRates))';
    % patch peak firing rate
    pfstats.patchPFR = max(patchRates);
    % patch mean firing rate
    pfstats.patchMFR = mean(patchRates);


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
end