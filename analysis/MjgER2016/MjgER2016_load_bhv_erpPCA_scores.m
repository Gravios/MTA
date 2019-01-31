% MjgER2016_load_bhv_erpPCA_scores.m
% DEFINITIONS
% drzrbhv field: directional rate zone restriced behavioral tuning curves

% DIMENSIONS 
% p := number of drz restricted bhv field types 
% T := number of Trials
% U := cumulative number of units
% V := number of factors (eigenvectors)
% D := number of valid elements in validDims{ind} 
% S := number of units in unitSubsets{ind}

% REQUIRED Vars
% pfindex     - Numeric[1x1]                   index of the drzbhv field
% pfd         - cellarray[T x p]{MTAApfs}      drz restricted bhv field
% units       - cellarray[1 x T]{NumericArray} list session of units
% unitSubsets - cellarray[1 x p]{NumericArray} concatenated list of units
% eigVec      - cellarray[1 x p]{NumericArray} eigenvectors
% validDims   - cellarray[1 x p]{NumericArray} valid eigenvector indicies


% CREATED Vars
% 
% NON - ESSENTIAL 
%     clu    - cluster ids
%     tlu    - trial ids
%     rind   - maps cluster order from MTAApfs object to Units order
%     D      - covariance matrix
%     LR     - eigenvectors
%     FSCFr  - pseudo-inverse of LR
%     rsMean   - rmaps shuffled mean
%     rsStd    - rmaps shuffled std
% 
% 
% AUXILLARAY vars
%     pfdShuffled
%     rmapsShuffledMean
%     rmapsShuffled
%     FSrM     - mean of 'shuffled' scores
%     FSrS     - std of 'shuffled' score
%     fsrsMean - 
%     fsrsStd
%     fsrsMean - 
%     fsrsStd
%
%
% OUTPUT Vars
%     rmaps    - matrix[D x S](numeric); rate maps corresponding to the valid eigenvector dims
%     FSrC     - matrix[U x V](Numeric); fscores
%     fsrcz    - matrix[U x V](Numeric); normalized fscores


version = '13';

%pfdParameters
% $$$ tags{end+1}        = ['HBPITCHxBPITCH_shuffled_v',version];
% $$$ fetSets{end+1}     = 'fet_HB_pitchB';
% $$$ fetInds{end+1}     = [1,2];
% $$$ boundaryLimits     = [-2.25,1;-0.5,2];
% $$$ binDims            = [0.2,0.2];
% $$$ SmoothingWeights   = [1.1,1.1];
% $$$ states{end+1}      = 'theta-groom-sit';
% $$$ ranges{end+1}      = [0,110];%[1100,1550];%1450];




filename = fullfile(MTASession([]).path.data,'analysis',['MjgER2016_erpPCA_scores_v' version '.mat']);

if exist(filename,'file'),
    load(filename);
else,
    pfindex = 1;

    if ~exist('Trials','var'),
        MjgER2016_load_data();
    end

    if ~exist('pfd','var'),
        [pfd,tags,eigVec,eigVar,eigScore,validDims,unitSubsets,unitIntersection,zrmMean,zrmStd] = ...
            req20180123_ver5(Trials,sessionList,version,false,false);
    end

    % LOAD shuffled drz restricted behavior fields 
    pfdShuffled =  cf(@(t,g) MTAApfs(t,'tag',g), Trials, repmat({['HBPITCHxBPITCH_v',version,'-shuffled']},size(Trials)));

    % GET bhv rate maps
    rmaps = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'), pfd(:,pfindex),units');
    clu =  cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:), pfd(:,pfindex),units');
    tlu =  cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
    rmaps = cat(2, rmaps{:});
    clu = cat(2, clu{:});
    tlu = cat(2, tlu{:});
    clu = [tlu',clu'];
    [~,rind] = sortrows(clu);
    rmaps = rmaps(:,rind);
    rmaps = rmaps(validDims{pfindex},unitSubsets{pfindex});
    rmaps(isnan(rmaps)) = 0;

    rmapsShuffledMean = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'), pfdShuffled',units');
    clu =  cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:), pfdShuffled',units');    
    tlu =  cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
    rmapsShuffledMean = cat(2, rmapsShuffledMean{:});   
    clu = cat(2, clu{:});    
    tlu = cat(2, tlu{:});    
    clu = [tlu',clu'];
    [~,rind] = sortrows(clu);
    rmapsShuffledMean = rmapsShuffledMean(:,rind);
    rmapsShuffledMean = rmapsShuffledMean(validDims{pfindex},unitSubsets{pfindex});
    rmapsShuffledMean(isnan(rmapsShuffledMean)) = 0;

    rmapsShuffled = cf(@(p,u) p.data.rateMap(:,ismember(p.data.clu,u),:), pfdShuffled',units');
    clu =  cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:), pfdShuffled',units');    
    tlu =  cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
    rmapsShuffled = cat(2, rmapsShuffled{:});   
    clu = cat(2, clu{:});    
    tlu = cat(2, tlu{:});    
    clu = [tlu',clu'];
    [~,rind] = sortrows(clu);
    rmapsShuffled = rmapsShuffled(:,rind,:);
    rmapsShuffled = rmapsShuffled(validDims{pfindex},unitSubsets{pfindex},:);
    rmapsShuffled(isnan(rmapsShuffled)) = 0;

    D = cov(rmapsShuffledMean');
    LR = eigVec{pfindex};
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
    for i = 1:pfd{1}.parameters.numIter
        FSrS(:,:,end+1) =  ((rmapsShuffled(:,:,i)'-rsMean)./rsStd) * FSCFr;
    end
    fsrsMean = mean(FSrS,3);
    fsrsStd = std(FSrS,[],3);

    fsrcz = (FSrC-fsrsMean)./fsrsStd;
     
    save(filename,...
         'pfdShuffled',...
         'rmapsShuffledMean',...
         'rmapsShuffled',...
         'FSCFr',...
         'FSrM',...
         'FSrS',...
         'fsrsMean',...
         'fsrsStd',...
         'fsrsMean',...
         'fsrsStd',...
         'rmaps',...
         'FSrC',...
         'fsrcz',...
         '-v7.3'...
         );
end