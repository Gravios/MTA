function pRHM(Trial)

xyz = Trial.load('xyz');
ang = Trial.ang.copy;
ang.create(Trial,xyz);
rhm = fet_rhm(Trial,[],'wcsd');
ang.resample(rhm);
figure
subplot(121)
hist2([log10(rhm(Trial.stc{'h'},40)),ang(Trial.stc{'h'},5,7,2)],-10:.2:-1,-1.5:.05:1.5),caxis([0,40])
subplot(122)
hist2([log10(rhm(Trial.stc{'l'},40)),ang(Trial.stc{'l'},5,7,2)],-10:.2:-1,-1.5:.05:1.5),caxis([0,40])

