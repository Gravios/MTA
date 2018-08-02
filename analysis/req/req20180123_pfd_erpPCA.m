function [LR,FSr,VT,unitSubset,vDims,zrmMean,zrmStd] = req20180123_pfd_erpPCA(pfd,units,range,pfindex,numComp,overwrite)

global MTA_PROJECT_PATH


filePath = fullfile(MTA_PROJECT_PATH,'analysis',['req20180123_pfd_erpPCA-',pfd{1,pfindex}.tag,'.mat']);


if ~exist(filePath,'file')||overwrite,
    assert(range(1)<range(2));

% GET bhvfield maps
    rmaps = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'), pfd(:,pfindex),units');
    clu =  cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:), pfd(:,pfindex),units');    
    tlu =  cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
    rmaps = cat(2, rmaps{:});   
    clu = cat(2, clu{:});    
    tlu = cat(2, tlu{:});    
    clu = [tlu',clu'];
    
    [~,rind] = sortrows(clu);
    rmaps = rmaps(:,rind);
    
% SUM total non nan values per map
    smnd  = sum(isnan(rmaps),1);
    unitSubset = find(...                               % SELECT units with good sampling
                      range(1)<smnd&smnd<range(2)&...   %   SELECT unit subset by range of valid element count 
                      max(rmaps)>4);                    %   SELECT unit subset by max rate
    zrmaps = rmaps(:,unitSubset);                       % SELECT a subset of rate maps by units
    zdims = size(zrmaps);                               % NOTE the dimensions of the original vector space
    zrmaps(isnan(zrmaps)) = 0;                          % SET all nan valued elements to zeros
    vDims = sum(zrmaps==0,2)<zdims(2)/3;                % SELECT subspace {2/3 non-zero samples} 
    zrmaps = zrmaps(vDims,:);                           % REDUCE to selected subspace
                                                        % DECOMPOSE rate maps
    
    %zrmaps = bsxfun(@rdivide,zrmaps,mean(prctile(zrmaps,90)));
    [~,LR,FSr,VT] = erpPCA(zrmaps',numComp);            % COMPUTE covariance-based PCA with Varimax rotation

    zrmMean = mean(zrmaps,2);
    zrmStd  = std(zrmaps,[],2);
    
    save(filePath,'LR','FSr','VT','unitSubset','vDims','zrmMean','zrmStd');
else
    load(filePath);
end

