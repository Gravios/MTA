

MjgER2016_load_data();
MjgER2016_supfig_bfrmCrossValidation();


% COMPUTE complex bhv rate map (over theta phase)
csmEigScrComplex = sum(reshape(bsxfun(@times,reshape(bsxfun(@plus,csmEigScrs,abs(min(csmEigScrs,[],3))),[],size(csmEigScrs,3)),exp(-i.*pfs.adata.bins{2})'),size(csmEigScrs)),3)./sum(bsxfun(@plus,csmEigScrs,abs(min(csmEigScrs,[],3))),3);


csmEigScrs = eigScrsNnmf;
csmEigVecs = eigVecsNnmf;
csmEigVars = eigVarsNnmf;
csmEigVars = bsxfun(@rdivide,csmEigVars,sum(csmEigVars,2,'omitnan'));
for j = 1:size(csmEigVars,1),
    [csmEigVars(j,:),sind] = sort(csmEigVars(j,:),2,'descend');
    csmEigScrs(j,:,:) = csmEigScrs(j,:,sind);
    csmEigVecs(j,:,:,:) = csmEigVecs(j,sind,:,:);
end
% COMPUTE complex nnmf score 
csmEigScrComplex = sq(sum(csmEigScrs.*repmat(exp(-i.*pfs.adata.bins{2})',[size(csmEigScrs,1),1, ...
                    size(csmEigScrs,3)]),2)./sum(csmEigScrs,2));



ind = diag(rhoB)>0.5&diag(rhoS)>0.5;
figure,plot3(fsrcz(ind,1),fsrcz(ind,2),fsrcz(ind,3),'.')


figure,
for u = unitSubset,
for j = 1:3,    
    subplot(2,3,j);  imagescnan({pfs.adata.bins{[1,3]},sq(eigVecsFact(u,j,:,:))'},[],[],true);axis('xy');
    title(num2str(eigVarsFact(u,j)));
    subplot(2,3,j+3);imagescnan({pfs.adata.bins{[1,3]},sq(eigVecsNnmf(u,j,:,:))'},[],[],true);axis('xy');
    title(num2str(eigVarsNnmf(u,j)));    
end
waitforbuttonpress();
end
    




csthresh = 2.1;
csoffset = 1.1;
cluSessionMapSubset = cluSessionMap(unitSubset,:);
cluSessionMapSubset_R = cluSessionMapSubset(  fsrcz(:,1)>csthresh                               ...
                                            & fsrcz(:,2)<(csthresh-csoffset)                    ...
                                            & fsrcz(:,3)<(csthresh-csoffset),:);
cluSessionMapSubset_L = cluSessionMapSubset(  fsrcz(:,3)>csthresh                               ...
                                            & fsrcz(:,1)<(csthresh-csoffset)                    ...
                                            & fsrcz(:,2)<(csthresh-csoffset),:);
cluSessionMapSubset_H = cluSessionMapSubset(  fsrcz(:,2)>csthresh                               ...
                                            & fsrcz(:,1)<(csthresh-csoffset)                    ...
                                            & fsrcz(:,3)<(csthresh-csoffset),:);
cluSessionMapSubset_N = cluSessionMapSubset(  fsrcz(:,2)<csthresh                               ...
                                            & fsrcz(:,1)<csthresh                               ...
                                            & fsrcz(:,3)<csthresh,:);
cluSessionMapSubset_C = cluSessionMapSubset(~ismember(cluSessionMapSubset,                      ...
                                                      [cluSessionMapSubset_L;                   ...
                                                       cluSessionMapSubset_H;                   ...
                                                       cluSessionMapSubset_R;                   ...
                                                       cluSessionMapSubset_N],                  ...
                                                      'rows'),                                  ...
                                            :);
cluSessionMapSubsets = {cluSessionMapSubset_R,cluSessionMapSubset_H,cluSessionMapSubset_L,cluSessionMapSubset_N};


sesIds = [3:5,8:12,17:23];
sesIds = [8:12,17:23];
sesInds = ismember(cluSessionMap(:,1),sesIds);
reqInds = ismember(cluSessionMap,cluSessionMap(unitSubset,:),'rows');
indR = ismember(cluSessionMap,cluSessionMapSubset_R,'rows');
indH = ismember(cluSessionMap,cluSessionMapSubset_H,'rows');
indL = ismember(cluSessionMap,cluSessionMapSubset_L,'rows');
indN = ismember(cluSessionMap,cluSessionMapSubset_N,'rows');
indC = ismember(cluSessionMap,cluSessionMapSubset_C,'rows');


ind = diag(rhoB)>0.5&diag(rhoS)>0.5;
guinds = false([size(cluSessionMap,1),1]);
guinds(unitSubset) = logical(ind)&csmEigVars(unitSubset,1)>0.4&csmEigVars(unitSubset,2)<0.4;
%guinds = true([size(cluSessionMap,1),1]);
guinds = guinds&reqInds&sesInds;
figure,
subplot(531);hist(angle(csmEigScrComplex(indR&guinds,1)),-pi:pi/6:pi);
subplot(532);hist(angle(csmEigScrComplex(indR&guinds,2)),-pi:pi/6:pi);
subplot(533);hist(angle(csmEigScrComplex(indR&guinds,3)),-pi:pi/6:pi);
subplot(534);hist(angle(csmEigScrComplex(indH&guinds,1)),-pi:pi/6:pi);
subplot(535);hist(angle(csmEigScrComplex(indH&guinds,2)),-pi:pi/6:pi);
subplot(536);hist(angle(csmEigScrComplex(indH&guinds,3)),-pi:pi/6:pi);
subplot(537);hist(angle(csmEigScrComplex(indL&guinds,1)),-pi:pi/6:pi);
subplot(538);hist(angle(csmEigScrComplex(indL&guinds,2)),-pi:pi/6:pi);
subplot(539);hist(angle(csmEigScrComplex(indL&guinds,3)),-pi:pi/6:pi);
subplot(5,3,10);hist(angle(csmEigScrComplex(indN&guinds,1)),-pi:pi/6:pi);
subplot(5,3,11);hist(angle(csmEigScrComplex(indN&guinds,2)),-pi:pi/6:pi);
subplot(5,3,12);hist(angle(csmEigScrComplex(indN&guinds,3)),-pi:pi/6:pi);
subplot(5,3,13);hist(angle(csmEigScrComplex(indC&guinds,1)),-pi:pi/6:pi);
subplot(5,3,14);hist(angle(csmEigScrComplex(indC&guinds,2)),-pi:pi/6:pi);
subplot(5,3,15);hist(angle(csmEigScrComplex(indC&guinds,3)),-pi:pi/6:pi);
af(@(a) ylim(a,[0,25]), findobj(gcf,'Type','Axes'));
af(@(a) xlim(a,[-pi,pi]), findobj(gcf,'Type','Axes'));

cf(@(t) t.load('nq'), Trials);
nq = cf(@(t) StructArray(t.nq,1), Trials)
nq = cf(@(n,u) n(u), nq,units);
nq = cat(1,nq{:});

% GET rhoB from MjgER2016_supfig_bfrmCrossValidation

% SELECT behaviorally and spatially stable units from the unitSubset;
ind = diag(rhoB)>0.5&diag(rhoS)>0.5;
figure,plot3(fsrcz(ind,1),fsrcz(ind,2),fsrcz(ind,3),'.')







