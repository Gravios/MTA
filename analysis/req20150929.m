
%% REQ 201509299
%
% REQ 1  Mahalanobis distance 
% REQ 2
% REQ 4


Trial = MTATrial('jg05-20120317');
nsr = 15;
[fts,Nmean,Nstd] = fet_tsne(Trial,nsr,true);
ind = Trial.stc{'w'};
cvs = fts(ind,:)'*fts(ind,:);
ms = mean(fts(ind,:));

md = zeros([fts.size(1),1]);
for i= 1:fts.size(1),
md(i) = sqrt((fts(i,:)-ms)/cvs*(fts(i,:)-ms)');
end



figure
hold on,plot(log10(md))
Lines(Trial.stc{'r',nsr}(:),[],'r');

mmd = MTADxyz('data',log10(md),'sampleRate',nsr);


eds = linspace(-2,-.2,100);
figure,hold on
ind = Trial.stc{'r'};
ha = bar(eds,histc(mmd(ind),eds),'histc');
ha.FaceColor = 'r';
ha.FaceAlpha = .5;
ha.EdgeColor = 'r';

ind = Trial.stc{'a-r'};
hs = bar(eds,histc(mmd(ind),eds),'histc');
hs.FaceColor = 'c';
hs.FaceAlpha = .5;
hs.EdgeColor = 'c';



SessionList.m