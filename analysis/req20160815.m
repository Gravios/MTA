Trial = MTATrial('Ed03-20140624');

rhm = fet_rhm(Trial,[],'mta');
frhm = rhm.copy;
frhm.filter('ButFilter',3,[6,15],'bandpass');
prhm = rhm.phase([5,15]);

xyz = Trial.load('xyz');

fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.4,'low');
vxy = fxyz.vel([1],[1,2]);

ang = create(MTADang,Trial,xyz);

dang = xyz.copy;
dang.data = ang(:,1,4,3)-ang(:,2,4,3);
dang.filter('ButFilter',3,[3,8],'bandpass');
pdang = dang.phase([2,8]);

bang = xyz.copy;
bang.data = circ_dist(ang(:,1,4,1),ang(:,2,4,1));
bang.filter('ButFilter',3,[2,8],'bandpass');
pbang = bang.phase([2,8]);

ind = Trial.stc{'w'};
eds = linspace([-pi,pi,50]);
figure,
hist2([prhm(ind),pbang(ind)],eds,eds);

figure
hist2([prhm(ind),pdang(ind)],eds,eds);

figure
hist2([pdang(ind),pbang(ind)],eds,eds);


dmin = LocalMinimaN(-dang.data,-0.6,10);
figure,
subplot(1,2,1)
rose( prhm(dmin(:,1)))
subplot(1,2,2)
rose(pdang(dmin(:,1)))



dmin = LocalMinimaN(-bang.data,-0.01,10);
figure,
subplot(1,2,1)
rose( prhm(dmin(:,1)))
subplot(1,2,2)
rose(pdang(dmin(:,1)))


dmin = LocalMinimaN(-frhm.data,-.5,10);
subplot(1,2,1)
rose(pdang(dmin(:,1)))
subplot(1,2,2)
rose( prhm(dmin(:,1)))

figure,
rose(pdang(nniz(pdang)))
figure,
rose(prhm(nniz(prhm)))


figure,plot(prhm.data)
hold on,plot(pdang.data)






figure,



[ys,fs,ts] = fet_spec
