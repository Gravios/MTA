function pRHM(Trial)


[rhm,fs,ts] = fet_rhm(Trial,[],'mtchglong');

figure();
imagesc(ts,fs,log10(rhm.data)');
axis('xy');
colormap('jet');
caxis([-8,-4]);

% $$$ figure();
% $$$ subplot(121);
% $$$ hist2([log10(rhm(Trial.stc{'h'},40)),ang(Trial.stc{'h'},5,7,2)],-10:.2:-1,-1.5:.05:1.5),caxis([0,40]);
% $$$ subplot(122);
% $$$ hist2([log10(rhm(Trial.stc{'l'},40)),ang(Trial.stc{'l'},5,7,2)],-10:.2:-1,-1.5:.05:1.5),caxis([0,40]);

