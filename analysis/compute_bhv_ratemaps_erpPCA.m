function [LR,FSr,VT,unitSubset,validDims,zrmMean,zrmStd] = compute_bhv_ratemaps_erpPCA(pfs,units,varargin)
% function [LR,FSr,VT,unitSubset,validDims,zrmMean,zrmStd] = ...
%     comput_bhv_ratemap_erpPCA(pfs,units,numComp,overwrite)
%
% Compute Factors (varimax) over behavior ratemaps
%
%
% VARARGIN : 
%    pfs       - cellarray{MTAApfs}, Ratemap objects  
%    units     - cellarray{Array},   unit sets 
%    numComp   - Integer,            number of expected components 
%    range     - Integer,            number of NANs accepted
%    overwrite - Logical,            FLAG, overwrite save file
%
global MTA_PROJECT_PATH

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('numComp',                         5,                                           ...
                 'range',                           [],                                          ...
                 'overwrite',                       false                                        ...
);
[numComp,range,overwrite] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% TAGS ---------------------------------------------------------------------------------------------
pfsTag = DataHash(cf(@(p)  p.filename,  pfs));
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------

filePath = fullfile(MTA_PROJECT_PATH,'analysis',[mfilename,'-',pfsTag,'.mat']);

if ~exist(filePath,'file') || overwrite,
    assert(range(1)<range(2));

% GET bhvfield maps
    rmaps = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'), pfs,units);
    clu =  cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:), pfs,units);
    tlu =  cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);

    rmaps = cat(2, rmaps{:});   
    clu = cat(2, clu{:});
    tlu = cat(2, tlu{~cellfun(@isempty,tlu)});    
    clu = [tlu',clu'];

    if pfs{1}.parameters.numIter>1
        snr = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan')...
                 ./std(p.data.rateMap(:,ismember(p.data.clu,u),:),[],3,'omitnan'),  pfs,units');
        snr = cat(2, snr{:});               
        rmaps(snr>75|snr<1) = nan;
    end
    
    [~,rind] = sortrows(clu);
    rmaps = rmaps(:,rind);

    
% SUM total nan values per map
    smnd  = sum(~isnan(rmaps),1);
% SELECT units with good sampling    
% ------ unit subset by range of valid element count 
% ------ unit subset by max rate
    unitSubset = find(...                               
                      range(1)<smnd&smnd<range(2)&...   
                      max(rmaps)>4);                    
% SELECT a subset of rate maps by units    
    zrmaps = rmaps(:,unitSubset);                       
% NOTE the dimensions of the original vector space    
    zdims = size(zrmaps);                               
% SET all nan valued elements to zeros    
    zrmaps(isnan(zrmaps)) = 0;                          
% SELECT subspace {1/2 non-zero samples}     
    validDims = sum(zrmaps==0,2)<zdims(2)/10;               
% REDUCE to selected subspace    
    zrmaps = zrmaps(validDims,:);                           
% DECOMPOSE rate maps
    zrind = find(prctile(zrmaps,90)>2);
    unitSubset = unitSubset(zrind);
    zrmaps = zrmaps(:,zrind);

    [~,LR,FSr,VT] = erpPCA(zrmaps',numComp);            % COMPUTE covariance-based PCA with Varimax rotation

    zrmMean = mean(zrmaps,2);
    zrmStd  = std(zrmaps,[],2);
    
    save(filePath,'LR','FSr','VT','unitSubset','validDims','zrmMean','zrmStd');
else
    load(filePath);
end

% END MAIN -----------------------------------------------------------------------------------------

% DIAGNOSTIC plots ---------------------------------------------------------------------------------
% $$$ figure
% $$$ for i = 1:numComp,
% $$$     evmap = zeros(pfs{1}.adata.binSizes');
% $$$     evmap(validDims) = LR(:,i);
% $$$     subplot(1,numComp+1,i);
% $$$     imagesc(pfs{1}.adata.bins{:},evmap');axis('xy')
% $$$ end
% $$$ subplot(1,numComp+1,numComp+1);
% $$$ plot(VT(1:numComp,4),'*')
% $$$ 
% $$$ figure
% $$$ for i = 1:numel(unitSubset),,
% $$$     evmap = zeros(pfs{1}.adata.binSizes');
% $$$     evmap(validDims) = zrmaps(:,i);
% $$$ clf();
% $$$     imagesc(pfs{1}.adata.bins{:},evmap');axis('xy')    
% $$$     title(num2str([clu(unitSubset(i),:),max(evmap(:))]))
% $$$     waitforbuttonpress();
% $$$ end
% $$$ 
% $$$ figure,plot3(FSr(:,1),FSr(:,2),FSr(:,3),'.')
% $$$ ind = FSr(:,1)<0;
% $$$ figure,plot3(FSr(ind,2),FSr(ind,3),FSr(ind,4),'.')
% $$$ 
% $$$ figure
% $$$ for i = 1:numel(unitSubset),
% $$$     clf();
% $$$     evmap = nan(pfd{1}.adata.binSizes');
% $$$     evmap(vDims) = zrmaps(:,i);
% $$$     imagescnan({pfd{1}.adata.bins{:},evmap'},[],[],true);axis('xy')
% $$$     waitforbuttonpress();
% $$$ end

