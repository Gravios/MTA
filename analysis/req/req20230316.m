
tid = 20;
rigidBodyMarkers = {'head_back','head_left','head_front','head_right'};

Trial = Trials{tid}

xyz = Trial.load('xyz','trb');
xyz.filter('ButFilter',4,20,'low');

[rigidBodyBasis, headCenterOfMass] = compute_rigidbody_basis_timeseries( xyz, rigidBodyMarkers);

scom = sq(headCenterOfMass);
dcom = circshift(scom,-3)-circshift(scom,3);
vcom = clip(multiprod(rigidBodyBasis,dcom,[2,3],[2]),-100,100);

ind = 200:400
figure,plot3(vcom(ind,1),vcom(ind,2),vcom(ind,3));

ind = 200:800
figure,plot(vcom(ind,2),vcom(ind,3));

fvcom = vcom;
fvcom(nniz(vcom),:) = ButFilter(vcom(nniz(vcom),:),4,2.5/(120*0.5),'low');

ind = 1200:1400
figure();
hold('on');
plot3( vcom(ind,1), vcom(ind,2), vcom(ind,3),'b');
plot3(fvcom(ind,1),fvcom(ind,2),fvcom(ind,3),'r');

pch = fet_hbp_hba(Trial);


sind = cast([Trial.stc{'x+p-n&a'}],'TimeSeries');
sind = logical(sind.data);

figure
hist2([pch(sind,2)-Trial.meta.correction.headYaw,fvcom(sind,2)],linspace(-1.25,1.25,30),linspace(-15,15,32),'xprob');
caxis([0,0.1])

figure,
plot(vcom(:,2)-fvcom(:,2));
hold('on');
plot(vcom(:,3)-fvcom(:,3));


[ys,fs,ts] = mtcsdglong([vcom(:,2)-fvcom(:,2),vcom(:,3)-fvcom(:,3)],2^7,xyz.sampleRate,2^6,2^5,3,[],[],[2,18]);


figure,
subplot(211);
imagesc(ts,fs,phase(ys(:,:,1,2))');
colormap(gca(),'hsv');
axis('xy');
subplot(212);
imagesc(ts,fs,log10(real(ys(:,:,1,1).*ys(:,:,2,2)))');
colormap(gca(),'jet');
axis('xy');


ts = ts+(2^5/2)/xyz.sampleRate;
ssr = 1/diff(ts(1:2));
pad = round([ts(1),mod(xyz.size(1)-round(2^5/2),2^5)/xyz.sampleRate].*ssr)-[1,0];
szy = size(ys);
rhm = MTADlfp('data',cat(1,zeros([pad(1),szy(2:end)]),ys,zeros([pad(2),szy(2:end)])),'sampleRate',ssr);

ts = cat(2,([1:pad(1)]./ssr),ts',([pad(1)+size(ts,1)]+[1:pad(2)])./ssr)';

fpch = resample(copy(pch),rhm);


sind = resample(cast([Trial.stc{'x+p&a'}],'TimeSeries'),rhm);
sind = logical(sind.data);

figure
hist2([fpch(sind,2)-Trial.meta.correction.headYaw,phase(rhm(sind,5,1,2))],linspace(-1.25,1.25,32),linspace(-pi,pi,32),'xprob');

figure
hist2([fpch(sind,2)-Trial.meta.correction.headYaw,phase(rhm(sind,6,1,2))],linspace(-1.25,1.25,32),linspace(-pi,pi,32),'');


figure,plot(ts,phase(rhm(:,6,1,2)))


rrhm = resample(copy(rhm),pch);

hold('on');
plot([1:size(rrhm)]./rrhm.sampleRate,phase(rrhm(:,6,1,2)))


sind = cast([Trial.stc{'x&a'}],'TimeSeries');
sind = logical(sind.data)& sqrt(sum(xyz(:,5,[1,2]).^2,3))<300;

figure
hist2([pch(sind,2)-Trial.meta.correction.headYaw,circ_dist(phase(rrhm(sind,6,1,2)),pi)],linspace(-1.25,1.25,16),linspace(-pi,pi,16),'yprob');
colormap(gca,'jet')
caxis([0.01,0.1])

spk = Trial.load('spk',xyz.sampleRate,'x&a&t');

unit = 25;
[mr,mxp] = pft.maxRate(unit)
res = spk(unit);
hpf = sq(xyz(:,7,[1,2])-xyz(:,5,[1,2]));
hpf = bsxfun(@rdivide,hpf,sqrt(sum(hpf.^2,2)));
hpf = cat(3,hpf,multiprod(permute(repmat([0,1;-1,0],[1,1,size(xyz,1)]),[3,1,2]), hpf,[2,3],[2]));
pvec = multiprod(hpf,bsxfun(@minus,mxp,sq(xyz(:,5,[1,2]))),[2,3],[2]);


figure,
scatter(pch(res,2)-Trial.meta.correction.headYaw,circ_dist(phase(rrhm(res,7,1,2)),pi),20,pvec(res,2),'filled');;
colormap(gca,'jet')
caxis([-200,200]);



figure,
scatter(pch(res,2)-Trial.meta.correction.headYaw,vcom(res,2),20,pvec(res,2),'filled');;
colormap(gca,'jet')
caxis([-250,250]);



% need the simple feature allo centric point in the head's frame of reference



figure
hist2([real(log10(real(rrhm(sind,6,1,1).*rrhm(sind,6,2,2)))),circ_dist(phase(rrhm(sind,6,1,2)),pi)],linspace(-4.25,-0.5,16),linspace(-pi,pi,16),'yprob');


vang = atan2(vcom(:,3),vcom(:,2));
vang = circ_dist(circshift(vang,-1),circshift(vang,1));


sind = cast([Trial.stc{'x&a'}],'TimeSeries');
sind = logical(sind.data)& sqrt(sum(xyz(:,5,[1,2]).^2,3))<300;
figure,
plot(vang(sind),pch(sind,2),'.');




unit = [20];
[mr,mxp] = pft.maxRate(unit)
res = spk(unit);
hpf = sq(xyz(:,7,[1,2])-xyz(:,5,[1,2]));
hpf = bsxfun(@rdivide,hpf,sqrt(sum(hpf.^2,2)));
hpf = cat(3,hpf,multiprod(permute(repmat([0,1;-1,0],[1,1,size(xyz,1)]),[3,1,2]), hpf,[2,3],[2]));
pvec = multiprod(hpf,bsxfun(@minus,mxp,sq(xyz(:,5,[1,2]))),[2,3],[2]);


vang = circshift(atan2(vcom(:,2),vcom(:,1)),0);
figure,
subplot(211);
scatter(pch(res,2)-Trial.meta.correction.headYaw,vang(res),20,pvec(res,2),'filled');;
colormap(gca,'jet')
caxis([-250,250]);
subplot(212);
scatter(pch(res,2)-Trial.meta.correction.headYaw,vang(res),20,pvec(res,1),'filled');;
colormap(gca,'jet')
caxis([-250,250]);


xyz = preproc_xyz(Trial,'trb');

figure,
subplot(211);
sind = [Trial.stc{'hloc'}];
histogram(xyz(sind,'hcom',3),linspace(0,300,30));
subplot(212);
sind = [Trial.stc{'lloc'}];
histogram(xyz(sind,'hcom',3),linspace(0,300,30));

figure
sind = [Trial.stc{'loc'}];
histogram(xyz(sind,'hcom',3),linspace(0,300,30));


figure
subplot(211);
sind = [Trial.stc{'hpause'}];
histogram(xyz(sind,'hcom',3),linspace(0,300,30));
subplot(212);
sind = [Trial.stc{'lpause'}];
histogram(xyz(sind,'hcom',3),linspace(0,300,30));
