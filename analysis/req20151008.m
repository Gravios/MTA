
Trial = MTATrial('jg05-20120317');
fet = fet_tsne(Trial,Trial.xyz.sampleRate,false);

x = 4;
y = 15;
edx = linspace(-.8,2,100);
edy = linspace(-6,0,100);
ind = Trial.stc{'a-r-m-k'};
figure,
subplot(221)
hist2([log10(fet(ind,x)),log10(fet(ind,y))],edx,edy),caxis([0,250])
subplot(222)
ind = Trial.stc{'w'};
hist2([log10(fet(ind,x)),log10(fet(ind,y))],edx,edy),caxis([0,150])
subplot(224)
ind = Trial.stc{'n'};
hist2([log10(fet(ind,x)),log10(fet(ind,y))],edx,edy),caxis([0,150])

