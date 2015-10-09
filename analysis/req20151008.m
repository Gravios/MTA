
Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev2');
fet = fet_tsne(Trial,Trial.xyz.sampleRate,false);

x = 7;
y = 15;
edx = linspace(-.2,1,100);
%edx = linspace(-.8,2,100);
%edy = linspace(-.2,1,100);
edy = linspace(-6,0,100);

xl = .5;
yl = .5;
%yl = -2.5,'r');Lines(.5,[],'r');

figure,
subplot(221)
ind = Trial.stc{'a-r-m-k'};
hist2([(fet(ind,x)),log10(fet(ind,y))],edx,edy),caxis([0,250])
title(ind.label);
%Lines([],yl,'r');Lines(xl,[],'r');
subplot(222)
ind = Trial.stc{'w'};
hist2([(fet(ind,x)),log10(fet(ind,y))],edx,edy),caxis([0,150])
title(ind.label);
%Lines([],yl,'r');%Lines(xl,[],'r');
subplot(223)
%ind = [Trial.stc{'r'}(:,2)-20,Trial.stc{'r'}(:,2)+20];
ind = Trial.stc{'p'};
hist2([(fet(ind,x)),log10(fet(ind,y))],edx,edy),caxis([0,150])
title(ind.label);
%Lines([],yl,'r');%Lines(xl,[],'r');
subplot(224)
ind = Trial.stc{'n'};
hist2([(fet(ind,x)),log10(fet(ind,y))],edx,edy),caxis([0,150])
title(ind.label);
Lines([],yl,'r');Lines(xl,[],'r');


edx = linspace(-.8,2,100);
ind = Trial.stc{'a-r-w-k-m'};
bar(edx,histc(fet(ind,x),edx),'histc')
ind = Trial.stc{'w'};
bar(edx,histc(fet(ind,x),edx),'histc')
