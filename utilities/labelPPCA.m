function labelPPCA(Trial)

fet = fet_lgr(Trial);
fet.data = nunity(fet.data);
nind = nniz(fet);
[U,S,V] = svd(cov(fet(nind,:)));

hfig = figure(839);

cpnts = ClusterPP(hfig);

PC1 = [fet(nind,:)*V(:,1),fet(nind,:)*V(:,2)];
nind = nniz(PC1);

figure
hist2(PC1,[prctile(PC1(nind,1),[1,99]),100],linspace(-6,4.5,100)),caxis([0,300])



figure
per = Trial.stc{'r'};
hist2([fet(per,:)*V(:,1),fet(per,:)*V(:,2)],linspace(-8,5,100),linspace(-6,4.5,100)),caxis([0,300])

fet.data
nind = Trial.stc{'a-r'};
[U,S,V] = svd(cov(fet(nind,:)));
figure
hist2([fet(nind,:)*V(:,1),fet(nind,:)*V(:,2)],linspace(-8,5,100),linspace(-6,4.5,100)),caxis([0,300])


figure
per = Trial.stc{'n'};
hist2([fet(per,:)*V(:,1),fet(per,:)*V(:,2)],linspace(-8,5,100),linspace(-6,4.5,100)),caxis([0,300])

