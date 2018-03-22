function [LR,FSr,VT,unitSubset,validDims] = req20180123_pfd_erpPCA(pfd,units,range,pfindex,numComp,overwrite)

global MTA_PROJECT_PATH



filePath = fullfile(MTA_PROJECT_PATH,'analysis',['req20180123_pfd_erpPCA-',pfd{pfindex}.tag,'.mat']);


if ~exist(filePath,'file')||overwrite,
    assert(range(1)<range(2));

    rmaps = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'), pfd(:,pfindex),units');
    rmaps = cat(2, rmaps{:});
    [smnd,sind] = sort(sum(isnan(rmaps),1));
    unitSubset = sind(range(1)<smnd&smnd<range(2));
    zrmaps = rmaps(:,unitSubset);
    zdims = size(zrmaps);                               % NOTE the dimensions of the original vector space
    zrmaps(isnan(zrmaps)) = 0;                          % SET all nan valued elements to zeros
    validDims = sum(zrmaps==0,2)<zdims(2)/3;            % SELECT subspace {2/3 non-zero samples} 
    zrmaps = zrmaps(validDims,:);                       % REDUCE to selected subspace
                                                        % DECOMPOSE rate maps
    [~,LR,FSr,VT] = erpPCA(zrmaps',numComp);            % COMPUTE covariance-based PCA with Varimax rotation

    save(filePath,'LR','FSr','VT','unitSubset','validDims');
else
    load(filePath);
end

