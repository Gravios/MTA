function [rmaps,clu] = decapsulate_and_concatenate_mtaapfs(pfs,units)
%function [rmaps] = decapsulate_and_concatenate_mtaapfs(pfs,units)
%
% 

rmaps = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'),    pfs,units);
clu   = cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:),                          pfs,units);
tlu   = cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
rmaps = cat(2, rmaps{:});
clu   = cat(2, clu{:});
tlu   = cat(2, tlu{:});
clu   = [tlu',clu'];
[clu,rind] = sortrows(clu);
rmaps = rmaps(:,rind,:);

