Trial = MTATrial('jg05-20120310');

rhm = fet_rhm(Trial,[],'mta');
frhm = rhm.copy;
frhm.filter('ButFilter',3,[6,12],'bandpass');
prhm = frhm.phase([6,12]);

xyz = Trial.load('xyz');

fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.4,'low');
vxy = fxyz.vel([1],[1,2]);

ang = create(MTADang,Trial,xyz);

frange = [1,5];

dang = xyz.copy;
dang.data = ang(:,1,4,3)-ang(:,2,4,3);
dang.filter('ButFilter',3,frange,'bandpass');
%dang.data = [0;diff(dang.data)];
pdang = dang.phase(frange);

bang = xyz.copy;
bang.data = circ_dist(ang(:,1,4,1),ang(:,2,4,1));
bang.filter('ButFilter',3,frange,'bandpass');
%bang.data = [0;diff(bang.data)];
pbang = bang.phase(frange);

ind = Trial.stc{'w'};
ind = vxy.data>3;
eds = linspace([-pi,pi,50]);


figure,
hist2([prhm(ind),pbang(ind)],eds,eds);

figure
hist2([prhm(ind),pdang(ind)],eds,eds);

figure
hist2([pdang(ind),pbang(ind)],eds,eds);

dmin = LocalMinimaN(-dang.data,-1,10);
%dmin = LocalMinimaN(-dang.data,-0.2,10);
figure,
rose(prhm(dmin(:,1)))
% $$$ figure,plot(prhm.data)
% $$$ hold on,plot(pdang.data)
% $$$ Lines(Trial.stc{'w'}(:),[],'m');



dmin = LocalMinimaN(-frhm.data,-.5,10);
figure
rose(pdang(dmin(:,1)))
figure
rose(pbang(dmin(:,1)))




dmin = LocalMinimaN(-bang.data,-0.005,10);
figure,
subplot(1,2,1)
rose( prhm(dmin(:,1)))
subplot(1,2,2)
rose(pdang(dmin(:,1)))


figure,
rose(pdang(nniz(xyz)))
figure,
rose(prhm(nniz(xyz)))








figure,



[ys,fs,ts] = fet_spec
