function binEdges = get_bin_edges(pfs)
binEdges = {};
for e = 1:numel(pfs.parameters.binDims),
   binEdges{e} = [pfs.adata.bins{e}-pfs.parameters.binDims(e)/2;pfs.adata.bins{e}(end)+pfs.parameters.binDims(e)/2];
end