
%clear('-global','AP');
global AP


% LOAD Trial data
% Vars
MjgER2016_load_data();


sessionListName = 'MjgER2016';
tag             = 'hbpptbpFS1v3';

% COMPUTE behavior rate maps over theta phases
AP.req20181106.Trials          = Trials;
AP.req20181106.units           = units;
AP.req20181106.sessionListName = sessionListName;
AP.req20181106.tag             = tag;

% COMPUTE factors of behavior rate over theta phases
AP.req20181119.sessionListName = sessionListName;
AP.req20181119.tag             = tag;



% LOAD behavior fields
[pfd,tags,eigVec,eigVar,eigScore,validDims,unitSubsets,unitIntersection,zrmMean,zrmStd] = ...
    req20180123_ver5(Trials,sessionList,'13',false,false);

MjgER2016_load_bhv_erpPCA_scores();


% COMPUTE behavior rate maps over theta phases
% COMPUTE factors of behavior rate over theta phases
pfs = req20181106(); 
[eigVecs, eigScrs, eigVars, eSpi, FSrBhvThp, metaData] = req20181119(); 

uind = {ismember(cluSessionMap(:,1),[1:2]);...
        ismember(cluSessionMap(:,1),[3:5]);...
        ismember(cluSessionMap(:,1),[6:7]);...
        ismember(cluSessionMap(:,1),[8:12]);...
        ismember(cluSessionMap(:,1),[13:16]);...
        ismember(cluSessionMap(:,1),[17:23])};

sid = {'er01';'ER06';'Ed10';'jg04';'jg04';'jg05'};

csmUnitSubset = cluSessionMap(unitSubsets{1},:);
csmEigScrs = eigScrs(unitSubsets{1},:,:);
csmEigVecs = eigVecs(unitSubsets{1},:,:,:);
csmEigVars = eigVars(unitSubsets{1},:,:);



clear('i');

[~,csmPhzPrefMax] = max(sq(csmEigScrs(:,:,:)),[],3);

csmPhzPrefMean = angle(sum(reshape(bsxfun(@times,reshape(bsxfun(@plus,csmEigScrs,abs(min(csmEigScrs,[],3))),[],size(csmEigScrs,3)),exp(-i.*pfs{1}.adata.bins{2})'),size(csmEigScrs)),3)./sum(bsxfun(@plus,csmEigScrs,abs(min(csmEigScrs,[],3))),3));

figure,plot(pfs{1}.adata.bins{2}(abs(csmPhzPrefMax(:,1)))+randn(size(csmPhzPrefMax(:,1)))/10,...
            csmPhzPrefMean(:,1),'.');
% $$$ figure,plot(sq(csmEigScrs(:,1,:))')

figure,plot(log10(csmEigVars(:,1)),log10(csmEigVars(:,2)),'.')

figure();
hold('on');
plot(circ_dist(csmPhzPrefMean(:,1),csmPhzPrefMean(:,2)),csmEigVars(:,1),'.b')
plot(circ_dist(csmPhzPrefMean(:,1),csmPhzPrefMean(:,2)),csmEigVars(:,2),'.r')


figure,
plot(circ_mean([csmPhzPrefMean(:,1),csmPhzPrefMean(:,2)]')+pi,csmEigVars(:,1),'.b')
