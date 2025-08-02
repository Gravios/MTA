function [fsrcz,FSrC,rmaps,FSCFr,FSrM,FSrS,fsrsMean,fsrsStd,rmapsShuffledMean,rmapsShuffled] = ...
    compute_bhv_ratemaps_erpPCA_scores(varargin)
%function [fsrcz,FSrC,rmaps,FSCFr,FSrM,FSrS,fsrsMean,fsrsStd,rmapsShuffledMean,rmapsShuffled] = ...
%    compute_bhv_ratemaps_erpPCA_scores(Trials,Units,bfrm,bfrmShuffled,validDims,unitSubset)
% 
% INPUTS:
%    Trials,
%    units,
%    bfrm,
%    bfrmShuffled,
%    eigVecs,
%    validDims,
%    unitSubset,
%    overwrite,
%    tag
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


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('Trials',                        [],                                            ...
                 'units',                         [],                                            ...
                 'bfrm',                          [],                                            ...
                 'bfrmShuffled',                  [],                                            ...
                 'eigVecs',                       [],                                            ...
                 'validDims',                     [],                                            ...
                 'unitSubset',                    [],                                            ...
                 'overwrite',                     false,                                         ...
                 'tag',                           []                                             ...
);
[Trials, units, bfrm, bfrmShuffled, eigVecs, validDims, unitSubset, overwrite, tag] =            ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------

if isempty(tag)
    tag = DataHash(cf(@(p,ps)  [p.filename,ps.filename],  bfrm,bfrmShuffled));
end

filePath = fullfile(MTA_PROJECT_PATH,'analysis',[mfilename(),'-',tag,'.mat']);

if exist(filePath,'file') && ~overwrite,
    load(filePath);
else,
    % GET bhv rate maps
    
    [rmaps,clu]   = decapsulate_and_concatenate_mtaapfs(bfrm,units);    
    rmaps         = mean(rmaps,3,'omitnan');
    clu           = clu(unitSubset,:);        
    rmaps         = rmaps(validDims,unitSubset);
    rmaps(isnan(rmaps)) = 0;
    
    rmapsShuffled = decapsulate_and_concatenate_mtaapfs(bfrmShuffled,units);
    rmapsShuffled = rmapsShuffled(validDims,unitSubset,:);

    rmapsShuffled(~nniz(rmapsShuffled(:))) = 0;    
    rmapsShuffledMean = mean(rmapsShuffled,3,'omitnan');

    
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