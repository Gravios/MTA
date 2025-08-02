
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

pftp = cf(@(Trial) MTAApfs(Trial,'tag','tp'), Trials);

rmaps = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'), pftp,units);
clu =  cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:), pfd(:,1),units');
tlu =  cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
rmaps = cat(2, rmaps{:});
clu = cat(2, clu{:});
tlu = cat(2, tlu{:});
clu = [tlu',clu'];
[~,rind] = sortrows(clu);
rmaps = rmaps(:,rind);
rmaps(isnan(rmaps)) = 0;
rmaps = rmaps';
 



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
csmESpi = eSpi(unitSubsets{1},:,:);
csmThpRate = permute(rmaps(unitSubsets{1},:),[1,3,2]);



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

ind = csmEigVars(:,1)>40&csmEigVars(:,2)>30;
ind = csmEigVars(:,1)>65;
ind = nniz(csmPhzPrefMean(:,1:2));
alpha = linspace(-pi,pi,100);
cppmsVar = zeros([numel(alpha),2]);
for s = 1:numel(alpha),
    cppms = [circ_dist(csmPhzPrefMean(ind,1),alpha(s)),circ_dist(csmPhzPrefMean(ind,2),alpha(s))]; 
    swpInd = bsxfun(@plus,[1,2],bsxfun(@times,double(gt(abs(cppms(:,1)),abs(cppms(:,2)))),[1,-1]));
    for a = 1:size(cppms,1),
        cppms(a,:) = cppms(a,swpInd(a,:));
    end
    cppmsVar(s,:) = circ_var(cppms);
end





csmBsiCpx = sum(reshape(bsxfun(@times,reshape(bsxfun(@plus,csmESpi,abs(min(csmESpi,[],3))),[],size(csmESpi,3)),exp(-i.*pfs{1}.adata.bins{2})'),size(csmESpi)),3)./sum(bsxfun(@plus,csmESpi,abs(min(csmESpi,[],3))),3);


csmBsiPhz = angle(csmBsiCpx);
csmBsiMag = abs(csmBsiCpx);



csmThpCpx = sum(reshape(bsxfun(@times,reshape(bsxfun(@plus,csmThpRate,abs(min(csmThpRate,[],3))),[],size(csmThpRate,3)),exp(-i.*pfs{1}.adata.bins{2})'),size(csmThpRate)),3)./sum(bsxfun(@plus,csmThpRate,abs(min(csmThpRate,[],3))),3);

csmThpPhz = angle(csmThpCpx);
csmThpMag = abs(csmThpCpx);

csmESpi

figure,
uid = 6;
subplot(221);plot(csmBsiPhz(uind{uid}(unitSubsets{1})),csmBsiMag(uind{uid}(unitSubsets{1})),'.')
subplot(222);plot(csmThpPhz(uind{uid}(unitSubsets{1})),csmThpMag(uind{uid}(unitSubsets{1})),'.')
subplot(223);plot(csmBsiPhz(uind{uid}(unitSubsets{1})),csmThpPhz(uind{uid}(unitSubsets{1})),'.')
subplot(224);plot(csmBsiMag(uind{uid}(unitSubsets{1})),csmThpMag(uind{uid}(unitSubsets{1})),'.')

figure
for u = 1:numel(uind)
    subplot(numel(uind),3,u*3-2);
    ind = uind{u}(unitSubsets{1})&(csmEigVars(:,1)>40&csmEigVars(:,2)>20)&FSrC(:,2)<0;
    %bar(pfs.adata.bins{2},histc(pfs.adata.bins{2}(abs(csmPhzPrefMax(ind,1))),pfs.adata.bins{2}),'histc');
    bar(pfs{1}.adata.bins{2},histc(csmPhzPrefMean(ind,1),pfs{1}.adata.bins{2}),'histc');    
    axis('tight');    
    title(sid{u});
    subplot(numel(uind),3,u*3-1);
    %bar(pfs{1}.adata.bins{2},histc(pfs{1}.adata.bins{2}(abs(csmPhzPrefMax(ind,2))),pfs{1}.adata.bins{2}),'histc');    
    bar(pfs{1}.adata.bins{2},histc(csmPhzPrefMean(ind,2),pfs{1}.adata.bins{2}),'histc');
    axis('tight');
    title(sid{u});
    subplot(numel(uind),3,u*3);    
    ind = uind{u}(unitSubsets{1})&(csmEigVars(:,2)<20)&FSrC(:,2)<0;    
    bar(pfs{1}.adata.bins{2},histc(csmPhzPrefMean(ind,1),pfs{1}.adata.bins{2}),'histc');    
    axis('tight');    
    title(sid{u});
% $$$     
% $$$     hist2([csmPhzPrefMean(ind,1),...
% $$$            csmPhzPrefMean(ind,2)], ...
% $$$           pfs{1}.adata.bins{2},pfs{1}.adata.bins{2});
% $$$     hist2([pfs{1}.adata.bins{2}(abs(csmPhzPrefMax(ind,1))),...
% $$$            pfs{1}.adata.bins{2}(abs(csmPhzPrefMax(ind,2)))], ...
% $$$           pfs{1}.adata.bins{2},pfs{1}.adata.bins{2});
    
end

% FIGURE 6 Bhv theta decomposition Summary


