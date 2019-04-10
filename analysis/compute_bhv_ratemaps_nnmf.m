function [W,H,D,unitSubset,validDims,zrmMean,zrmStd] = compute_bhv_ratemaps_nnmf(pfs,units,varargin)
% function [LR,FSr,VT,unitSubset,validDims,zrmMean,zrmStd] = ...
%     comput_bhv_ratemap_erpPCA(pfs,units,numComp,overwrite)
%
% ARGS : 
%    pfs - cellarray(MTAApfs), all ratemap objects required 
%

global MTA_PROJECT_PATH

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('numComp',                         3,                                           ...
                 'range',                           [],                                          ...
                 'overwrite',                       false                                        ...
);
[numComp,range,overwrite] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% TAGS ---------------------------------------------------------------------------------------------
tag = DataHash(cf(@(p)  p.filename,  pfs));
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------

filePath = fullfile(MTA_PROJECT_PATH,'analysis',[mfilename,'-',tag,'.mat']);

if ~exist(filePath,'file') || overwrite,
    assert(range(1)<range(2));

% GET bhvfield maps
    rmaps = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'), pfs,units);
    clu =  cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:), pfs,units);    
    tlu =  cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);

    rmaps = cat(2, rmaps{:});   
    clu = cat(2, clu{:});    
    tlu = cat(2, tlu{:});    
    clu = [tlu',clu'];

    if pfs{1}.parameters.numIter>1
        snr = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan')...
                 ./std(p.data.rateMap(:,ismember(p.data.clu,u),:),[],3,'omitnan'),  pfs,units);
        snr = cat(2, snr{:});               
        rmaps(snr>75|snr<1) = nan;
    end
    
    [~,rind] = sortrows(clu);
    rmaps = rmaps(:,rind);

    
% SUM total nan values per map
    smnd  = sum(isnan(rmaps),1);
    unitSubset = find(...                               % SELECT units with good sampling
                      range(1)<smnd&smnd<range(2)&...   %   SELECT unit subset by range of valid element count 
                      max(rmaps)>2);                    %   SELECT unit subset by max rate
    zrmaps = rmaps(:,unitSubset);                       % SELECT a subset of rate maps by units
    zdims = size(zrmaps);                               % NOTE the dimensions of the original vector space
    zrmaps(isnan(zrmaps)) = 0;                          % SET all nan valued elements to zeros
    validDims = sum(zrmaps==0,2)<zdims(2)/3;                % SELECT subspace {1/2 non-zero samples} 
    zrmaps = zrmaps(validDims,:);                           % REDUCE to selected subspace
                                                        % DECOMPOSE rate maps
    
    %zrmaps = bsxfun(@rdivide,zrmaps,mean(prctile(zrmaps,90)));
    [W,H,D] = nnmf(zrmaps,k);
% $$$     W = {};H = {}; D = [];
% $$$     for k = 1:numComp
% $$$         [W{k},H{k},D(k)] = nnmf(zrmaps,k);
% $$$     end
    %[~,LR,FSr,VT] = erpPCA(zrmaps',numComp);            % COMPUTE covariance-based PCA with Varimax rotation

    zrmMean = mean(zrmaps,2);
    zrmStd  = std(zrmaps,[],2);
    
    save(filePath,'W','H','D','unitSubset','validDims','zrmMean','zrmStd');
else
    load(filePath);
end

% END MAIN -----------------------------------------------------------------------------------------


% DIAGNOSTIC plots ---------------------------------------------------------------------------------
% $$$ figure
% $$$ numComp =       3;
% $$$ for i = 1:numComp,
% $$$     evmap = zeros(pfs{1}.adata.binSizes');
% $$$     evmap(validDims) = W{numComp}(:,i);
% $$$     subplot(1,numComp+1,i);
% $$$     imagesc(pfs{1}.adata.bins{:},evmap');axis('xy')
% $$$ end
% $$$ subplot(1,numComp+1,numComp+1);
% $$$ %plot(VT(1:numComp,4),'*')

% $$$ figure
% $$$ for i = 1:numComp,
% $$$     evmap = zeros(pfs{1}.adata.binSizes');
% $$$     evmap(validDims) = W(:,i);
% $$$     subplot(1,numComp+1,i);
% $$$     imagesc(pfs{1}.adata.bins{:},evmap');axis('xy')
% $$$ end
% $$$ subplot(1,numComp+1,numComp+1);
% $$$ plot(VT(1:numComp,4),'*')

% $$$ figure
% $$$ for i = 1:numel(unitSubset),,
% $$$     evmap = zeros(pfs{1}.adata.binSizes');
% $$$     evmap(validDims) = zrmaps(:,i);
% $$$ clf();
% $$$     imagesc(pfs{1}.adata.bins{:},evmap');axis('xy')
% $$$     title(num2str(clu(unitSubset(i),:)))
% $$$     waitforbuttonpress();
% $$$ end
% $$$ figure,plot3(FSr(:,1),FSr(:,2),FSr(:,3),'.')
% $$$ ind = FSr(:,1)<0;
% $$$ figure,plot3(FSr(ind,2),FSr(ind,3),FSr(ind,4),'.')


% $$$ figure
% $$$ for i = 1:numel(unitSubset),
% $$$     clf();
% $$$     evmap = nan(pfd{1}.adata.binSizes');
% $$$     evmap(vDims) = zrmaps(:,i);
% $$$     imagescnan({pfd{1}.adata.bins{:},evmap'},[],[],true);axis('xy')
% $$$     waitforbuttonpress();
% $$$ end
% $$$ 
