function ratemap = collate_ratemaps(pfs,units,mask,rateOffset)
%function ratemap = collate_ratemaps(pfs,units,mask,rateOffset)
ratemap = zeros([prod(pfs.adata.binSizes),0]);
for u = 1:numel(units),
    trm = pfs.plot(units(u),1,false,[],false,0.25,false);
    ratemap = cat(2,ratemap,trm(:));
end
ratemap(~mask) = nan;
ratemap = ratemap+rateOffset;    

% $$$ ratemap(isnan(ratemap)) = 0;
% $$$ ratemap(~mask) = nan;
% $$$ ratemap = ratemap+rateOffset;    
