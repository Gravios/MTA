function labelPPCA(Trial)
%Testing Vars
Trial = MTATrial('jg05-20120317');
Trial.stc

fet = fet_lgr(Trial);
fet.data = nunity(fet.data);
nind = nniz(fet);
[U,S,V] = svd(cov(fet(nind,:)));

% $$$ hfig = figure(839);
% $$$ 
% $$$ cpnts = ClusterPP(hfig);


PC1 = Trial.xyz.copy;
PC1.data = [fet.data*V(:,1),fet.data*V(:,2)];
nind = nniz(PC1);

figure
p1b = prctile(PC1(nind,1),[1,99]);
p2b = prctile(PC1(nind,2),[1,99]);
hist2(PC1,linspace(p1b(1)-3,p1b(2)+3,100),linspace(p2b(1)-3,p2b(2)+3,100)),caxis([0,300])

figure,
hist2(PC1(Trial.stc{'r'},:),linspace(p1b(1)-3,p1b(2)+3,100),linspace(p2b(1)-3,p2b(2)+3,100)),caxis([0,300])


%% Remove Rearing
rper = Trial.stc{'r'}.cast('TimeSeries',fet.sampleRate);
fet.data(rper.data(1:end-1)==1) = nan;
nind = nniz(fet);
[U,S,V] = svd(cov(fet(nind,:)));

PC1 = Trial.xyz.copy;
PC1.data = [fet.data*V(:,1),fet.data*V(:,2)];
nind = nniz(PC1);

figure
p1b = prctile(PC1(nind,1),[1,99]);
p2b = prctile(PC1(nind,2),[1,99]);
hist2(PC1,linspace(p1b(1)-3,p1b(2)+3,100),linspace(p2b(1)-3,p2b(2)+3,100)),caxis([0,300])

figure,
hist2(PC1(Trial.stc{'m'},:),linspace(p1b(1)-3,p1b(2)+3,100),linspace(p2b(1)-3,p2b(2)+3,100)),caxis([0,300])


%% Remove Grooming
rper = Trial.stc{'m'}.cast('TimeSeries',fet.sampleRate);
fet.data(rper.data(1:end-1)==1) = nan;
nind = nniz(fet);
[U,S,V] = svd(cov(fet(nind,:)));

PC1 = Trial.xyz.copy;
PC1.data = [fet.data*V(:,1),fet.data*V(:,2)];
nind = nniz(PC1);

figure
p1b = prctile(PC1(nind,1),[1,99]);
p2b = prctile(PC1(nind,2),[1,99]);
hist2(PC1,linspace(p1b(1)-3,p1b(2)+3,100),linspace(p2b(1)-3,p2b(2)+3,100)),caxis([0,300])

figure,
hist2(PC1(Trial.stc{'w'},:),linspace(p1b(1)-3,p1b(2)+3,100),linspace(p2b(1)-3,p2b(2)+3,100)),caxis([0,300])


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

