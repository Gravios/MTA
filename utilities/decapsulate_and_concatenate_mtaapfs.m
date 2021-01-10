function [rmaps,clu] = decapsulate_and_concatenate_mtaapfs(pfs,units)
%function [rmaps] = decapsulate_and_concatenate_mtaapfs(pfs,units)
%
% 
if isempty(units)
    rmaps = {};
    clu = [];
    return;
end

tlu   = cf(@(i,u) repmat(i,[1,numel(u)]), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);

exu = ~cellfun(@isempty,units);
pfs = pfs(exu);
units = units(exu);

rmaps = cf(@(p,u) p.data.rateMap(:,ismember(p.data.clu,u),:),    pfs,units);
clu   = cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:),        pfs,units);

rmaps = cat(2, rmaps{:});
clu   = cat(2, clu{:});
tlu   = cat(2, tlu{:});
clu   = [tlu',clu'];
[clu,rind] = sortrows(clu);
rmaps = rmaps(:,rind,:);

