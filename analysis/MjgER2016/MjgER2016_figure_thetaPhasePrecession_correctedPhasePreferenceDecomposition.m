% MjgER2016_figure_thetaPhasePrecession_correctedPhasePreferenceDecomposition.m


sigma = 150;
pftHZTPD = cf(@(s) cf(@(T,u) MTAApfs(T,u,'tag',['ddtp-','s',num2str(sigma),'-',s]), Trials, units), stateLabels);
pftTPZDa = cf(@(pfs) ...
              cf(@(p,u) ...
                 p.data.rateMap(:,ismember(p.data.clu,u),:),  ...
                 pfs,units), ...
              pftHZTPD);

% $$$ u =3;figure;plot(sq(max(reshape(pftTPZDa{3}{20}(:,u,:),[40,16,101]))),'c');hold on, ...
% $$$    plot(sq(max(reshape(pftTPZDa{4}{20}(:,u,:),[40,16,101]))),'b');
% $$$ u =3;figure;plot(mean(RectFilter(sq(max(reshape(pftTPZDa{3}{20}(:,u,:),[40,16,101]))),3,1,'circular'),2),'c');
% $$$ hold('on'); plot(mean(RectFilter(sq(max(reshape(pftTPZDa{4}{20}(:,u,:),[40,16,101]))),3,1,'circular'),2),'b');

clu     = cf(@(pfs) cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:),                        pfs,units), pftHZTPD);
tlu     = cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
pftTPZDa = cf(@(pfs) cat(2,pfs{:}), pftTPZDa);
clu     = cf(@(c)   cat(2,c{:}),clu);
assert(all(cell2mat(cf(@(c) all(c), cf(@(c,o) c==o, clu,circshift(clu,1,2))))),'Clus do not match');
pftTPZDa = cat(4,pftTPZDa{:});
clu = clu{1};
tlu=cat(2,tlu{:});
clu=[tlu',clu'];
[clu,rind] = sortrows(clu);
clu     = clu(unitSubset,:);
pftTPZDa = pftTPZDa(:,rind,:,:);
pftTPZDa = pftTPZDa(:,unitSubset,:,:);
pftTPZDa = reshape(permute(pftTPZDa,[1,5,2,3,4]),[pftHZTPD{1}{1}.adata.binSizes',...
                    numel(unitSubset),pftHZTPD{1}{1}.parameters.numIter,numel(stateLabels)]);

% For each unit, for each state, for each phase, find the maximum rate and position along the gdz


pind =  all(sq(sum(any(isnan(pftTPZDa(:,:,:,1,2:4)),2)))<5,2);

[grate,gpos] = max(permute(sq(mean(RectFilter(permute(pftTPZDa(:,:,pind&sesInds(unitSubset),:,:),[2,1,3,4,5]),3,1,'circular'),4,'omitnan')),[2,1,3,4,5]));

% $$$ [grate,gpos] = max(permute(sq(mean(permute(pftTPZDa(:,:,pind&sesInds(unitSubset),:,:),[2,1,3,4,5]),4,'omitnan')),[2,1,3,4,5]));

grate = sq(grate(:,:,:,2:4));
gpos = sq(gpos);
[prate,pphz] = max(grate);
clusub = clu(sesInds(unitSubset),:);

npftTPZa = reshape(grate,16,[]);
npftTPZa(~nniz(npftTPZa(:))) = 0;
nind = sum(npftTPZa==0)==16;
npftTPZa(:,nind) =[];

[W,H,D] = nnmf(npftTPZa,3);
% $$$ figure,plot(W(phzOrder,:))
fscrCPPD = zeros([numel(nind),3]);
fscrCPPD(~nind,:) = H';
fscrCPPD = reshape(fscrCPPD,[],3,3);

msr = sq(max(sq(grate(phzOrder,:,:))));
[msr,msi] = sort(msr,2,'descend');
dsr = sq(mean(sq(grate(phzOrder,:,:))));
[dsr,dsi] = sort(dsr,2,'descend');

