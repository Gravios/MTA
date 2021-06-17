function expectedRate = get_expected_rate(Pfs,units,featureSet,state);


eds = {};
for e = 1:numel(Pfs.parameters.binDims),
   eds{e} = [Pfs.adata.bins{e}-Pfs.parameters.binDims(e)/2;Pfs.adata.bins{e}(end)+Pfs.parameters.binDims(e)/2];
end

posInd = {};
for d = 1:size(featureSet,2),
    posInd{d} = discretize(featureSet(state,d),eds{d});
end

pfsRateIndex = sub2ind(Pfs.adata.binSizes',posInd{:});

expectedRate = nan([numel(posInd{1}),numel(units)]);
for u = 1:numel(units)
    rmap = Pfs.data.rateMap(:,Pfs.data.clu==units(u));
    npri = nniz(pfsRateIndex);
    expectedRate(npri,u) = rmap(pfsRateIndex(npri));
end
