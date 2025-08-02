function [fsrcz,FSrC,rmaps,FSCFr,FSrM,FSrS,fsrsMean,fsrsStd,rmapsShuffledMean,rmapsShuffled] = ...
    compute_bhv_ratemaps_nnmf_scores(Trials,units,bfrm,bfrmShuffled,eigVecs,validDims,unitSubset)
%function [fsrcz,FSrC,rmaps,FSCFr,FSrM,FSrS,fsrsMean,fsrsStd,rmapsShuffledMean,rmapsShuffled] = ...
%    compute_bhv_ratemaps_erpPCA_scores(Trials,Units,bfrm,bfrmShuffled,validDims,unitSubset)
%
%
% drzrbhv field: directional rate zone restriced behavioral tuning curves
% DIMENSIONS 
% p := number of drz restricted bhv field types 
% T := number of Trials
% U := cumulative number of units
% V := number of factors (eigenvectors)
% D := number of valid elements in validDims{ind} 
% S := number of units in unitSubset{ind}
%
% OUTPUT Vars
%     rmaps    - matrix[D x S](numeric); rate maps corresponding to the valid eigenvector dims
%     FSrC     - matrix[U x V](Numeric); fscores
%     fsrcz    - matrix[U x V](Numeric); normalized fscores


global MTA_PROJECT_PATH


% TAGS ---------------------------------------------------------------------------------------------
pfsTag = DataHash(cf(@(p,ps)  [p.filename,ps.filename],  bfrm,bfrmShuffled));
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------

filePath = fullfile(MTA_PROJECT_PATH,'analysis',[mfilename,'-',pfsTag,'.mat']);


if exist(filePath,'file'),
    load(filePath);
else,
    pfindex = 1;

    % [rmaps,clu] = concatenate_set_pfs(pfs,units);
    % GET bhv rate maps
    rmaps = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'),  bfrm,units);
    clu   = cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:),                        bfrm,units);
    tlu   = cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
    rmaps = cat(2, rmaps{:});
    clu   = cat(2, clu{:});
    tlu   = cat(2, tlu{:});
    clu   = [tlu',clu'];
    [clu,rind] = sortrows(clu);
    clu   = clu(unitSubset,:);
    rmaps = rmaps(:,rind);
    rmaps = rmaps(validDims,unitSubset);
    rmaps(isnan(rmaps)) = 0;

    rmapsShuffledMean = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'), bfrmShuffled,units);
    clu =  cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:), bfrmShuffled,units);    
    tlu =  cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
    rmapsShuffledMean = cat(2, rmapsShuffledMean{:});   
    clu = cat(2, clu{:});    
    tlu = cat(2, tlu{:});    
    clu = [tlu',clu'];
    [~,rind] = sortrows(clu);
    rmapsShuffledMean = rmapsShuffledMean(:,rind);
    rmapsShuffledMean = rmapsShuffledMean(validDims,unitSubset);
    rmapsShuffledMean(isnan(rmapsShuffledMean)) = 0;

    rmapsShuffled = cf(@(p,u) p.data.rateMap(:,ismember(p.data.clu,u),:), bfrmShuffled,units);
    clu =  cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:), bfrmShuffled,units);    
    tlu =  cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
    rmapsShuffled = cat(2, rmapsShuffled{:});   
    clu = cat(2, clu{:});    
    tlu = cat(2, tlu{:});    
    clu = [tlu',clu'];
    [~,rind] = sortrows(clu);
    rmapsShuffled = rmapsShuffled(:,rind,:);
    rmapsShuffled = rmapsShuffled(validDims,unitSubset,:);
    rmapsShuffled(isnan(rmapsShuffled)) = 0;

    D = cov(rmapsShuffledMean');
    LR = eigVecs;
% COMPUTE rotated FS coefficients        
    FSCFr = LR * inv(LR' * LR);          % this is pseudo-inverse of LR
                                         % rescale rotated FS coefficients by the corresponding SDs 
    rk = size(FSCFr,2);%rank(D,1e-4);       % why not on D? would save on corrcoef(X) computation
    FSCFr = FSCFr .* repmat(sqrt(diag(D)),1,rk); % compute rotated factor scores from the normalized raw
                                                 % data and  the corresponding rescaled factor score coefficients
    rsMean = mean(rmapsShuffledMean');
    rsStd  = std( rmapsShuffledMean');

% MEAN shuffled score
    FSrM =  ((rmapsShuffledMean'-rsMean)./rsStd) * FSCFr;
% MEAN normal score
    FSrC =  ((rmaps'-rsMean)./rsStd) * FSCFr;

    FSrS = [];
    for i = 1:bfrm{1}.parameters.numIter
        FSrS(:,:,end+1) =  ((rmapsShuffled(:,:,i)'-rsMean)./rsStd) * FSCFr;
    end
    fsrsMean = mean(FSrS,3);
    fsrsStd = std(FSrS,[],3);

    fsrcz = (FSrC-fsrsMean)./fsrsStd;
     
    save(filePath,...   
         'rmapsShuffledMean',...
         'rmapsShuffled',...
         'FSCFr',...
         'FSrM',...
         'FSrS',...
         'fsrsMean',...
         'fsrsStd',...
         'rmaps',...
         'FSrC',...
         'fsrcz',...
         '-v7.3'...
         );
end