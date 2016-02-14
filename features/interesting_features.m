Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev2');
xyz = Trial.load('xyz');
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,.5,'low');
ang = create(MTADang,Trial,fxyz);



name = 'spine and head wag'; label = 'shwag'; key = 'w';
zv = MTADfet.encapsulate(Trial,...
                         [circ_dist(ang(:,1,3,1),ang(:,2,4,1)),...
                          circ_dist(ang(:,2,4,1),ang(:,3,7,1))],...
                         xyz.sampleRate,...
                         name,label,key);


dspec = struct('nFFT',2^8,'Fs',zv.sampleRate,...
               'WinLength',2^7,'nOverlap',2^7*.875,...
               'FreqRange',[1,30]);
[ys,fs,ts] = fet_spec(Trial,zv,'mtchglong',true,[],dspec,true);

figure,imagesc(ts,fs,ys(:,:,1,2)'),axis xy,colormap jet
figure,imagesc(ts,fs,log10(ys(:,:,2,2))'),axis xy,colormap jet

eds = linspace(-9,1,100);
figure,hold on
ind = Trial.stc{'a-p'};
hs = bar(eds,histc(log10(ys(ind,2,2,2)),eds),'histc');hs.FaceColor = 'c';hs.FaceAlpha = .5;
ind = Trial.stc{'p'};
hs = bar(eds,histc(log10(ys(ind,2,2,2)),eds),'histc');hs.FaceColor = 'r';hs.FaceAlpha = .5;

name = 'zv and wag pw'; label = 'zvwp'; key = 'w';
zfet = MTADfet.encapsulate(Trial,...
                         log10(nanmedian(ys(:,1:fs<7,1,1),2)),...
                         ys.sampleRate,...
                         name,label,key);
zfet.data(~nniz(zfet),:)=-20;
zfet.resample(bfet);



%exp in normalizing distances
mar = [1,4];
NBINS = 100;
VEL_HISTOGRAM_BOUNDARIES = linspace(-3,2,NBINS);
ANG_HISTOGRAM_BOUNDARIES = linspace(-pi/2,pi/2,NBINS);


Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev2');
xyz = Trial.load('xyz');
xyz.addMarker('bcom',[.7,0,.7],{},...
              xyz.com(xyz.model.rb({'spine_lower','pelvis_root','spine_middle'})));
xyz.addMarker('hcom',[.7,0,.7],{},...
              xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));
ang = create(MTADang,Trial,xyz);
vxy = xyz.copy;
vxy.filter('ButFilter',3,2.4,'low');
vxy = vxy.vel({'bcom','hcom'},[1,2]);
vxy.data(vxy.data<1e-3) = 1e-3;
vxy.data = log10(vxy.data);
ind = Trial.stc('a-m-r').cast('TimeSeries');
[~,ind_v_b] = histc(vxy(:,'bcom'),VEL_HISTOGRAM_BOUNDARIES);
[~,ind_v_h] = histc(vxy(:,'hcom'),VEL_HISTOGRAM_BOUNDARIES);
[~,ind_p_s] = histc(ang(:,3,4,2),ANG_HISTOGRAM_BOUNDARIES);
mind = nniz([ind_v_b,ind_v_h,ind_p_s])&ind.data;
manifoldIndex = [ind_v_b(mind),ind_v_h(mind),ind_p_s(mind)];
acvar = ang(mind,mar(1),mar(2),3);
%acvar = xyz(mind,mar,3);
mz = accumarray(manifoldIndex,acvar,[NBINS,NBINS,NBINS],@nanmean);
vz = accumarray(manifoldIndex,acvar,[NBINS,NBINS,NBINS],@nanstd);


% $$$ figure,subplot(1,2,1),imagesc(mz'),subplot(1,2,2),imagesc(vz')
% $$$ figure,hist(mz(nniz(mz(:))&vz(:)<10),100)


Trial = MTATrial('Ed03-20140624');
Trial.load('stc','hand_labeled_rev2_alt');
%Trial = MTATrial('Ed01-20140707');
%Trial.load('stc','hand_labeled_rev1');


xyz = Trial.load('xyz');
xyz.addMarker('bcom',[.7,0,.7],{},...
              xyz.com(xyz.model.rb({'spine_lower','pelvis_root','spine_middle'})));
xyz.addMarker('hcom',[.7,0,.7],{},...
              xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));
ang = create(MTADang,Trial,xyz);
vxy = xyz.copy;
vxy.filter('ButFilter',3,2.4,'low');
vxy = vxy.vel({'bcom','hcom'},[1,2]);
vxy.data(vxy.data<1e-3) = 1e-3;
vxy.data = log10(vxy.data);
ind = Trial.stc{'a-m-r'}.cast('TimeSeries',xyz.sampleRate);
[~,ind_v_b] = histc(vxy(:,'bcom'),VEL_HISTOGRAM_BOUNDARIES);
[~,ind_v_h] = histc(vxy(:,'hcom'),VEL_HISTOGRAM_BOUNDARIES);
mind = nniz([ind_v_b,ind_v_h])&ind.data;
manifoldIndex = [ind_v_b(mind),ind_v_h(mind),ind_p_s(mind)];
acvar = ang(mind,mar(1),mar(2),3);
%acvar = xyz(mind,mar,3);
mzo = accumarray(manifoldIndex,acvar,[NBINS,NBINS,NBINS],@nanmean);
vzo = accumarray(manifoldIndex,acvar,[NBINS,NBINS,NBINS],@nanstd);



%figure,subplot(1,2,1),imagesc(mz'),axis xy,subplot(1,2,2), ...
%    imagesc(mzo'),axis xy

ind = vz(:)<10&vzo(:)<.2&nniz(mzo(:))&nniz(mz(:));
figure,hist(mz(ind)-mzo(ind),100)
Lines(nanmedian(mz(ind)-mzo(ind)),[],'r');


figure,hist(ang(mind,mar(1),mar(2),3),100);



Trial = MTATrial('jg05-20120317');
Trial = MTATrial('Ed03-20140624');
Trial = MTATrial('Ed01-20140707');
xyz = Trial.load('xyz');
vxy = xyz.copy;
vxy.filter('ButFilter',3,2.4,'low');
vxy = vxy.vel({'spine_lower'},[1,2]);
vxy.data(vxy.data<1e-3) = 1e-3;
vxy.data = log10(vxy.data);
figure,hist(xyz(vxy(:,1)<-.5&nniz(xyz),1,3),100)









eds = linspace(-pi/2,pi/2,100);
figure,hold on 

Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev2');

xyz = Trial.load('xyz');
xyz.addMarker('bcom',[.7,0,.7],{},...
              xyz.com(xyz.model.rb({'spine_lower','pelvis_root','spine_middle'})));
xyz.addMarker('hcom',[.7,0,.7],{},...
              xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));
ang = create(MTADang,Trial,xyz);

ind = Trial.stc{'w'};
hs = bar(eds,histc(ang(ind,4,7,2),eds),'histc'); hs.FaceColor = 'c'; hs.FaceAlpha = .5;
Trial = MTATrial('Ed03-20140624');
Trial.load('stc','hand_labeled_rev2_alt');

xyz = Trial.load('xyz');
xyz.addMarker('bcom',[.7,0,.7],{},...
              xyz.com(xyz.model.rb({'spine_lower','pelvis_root','spine_middle'})));
xyz.addMarker('hcom',[.7,0,.7],{},...
              xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));
ang = create(MTADang,Trial,xyz);

ind = Trial.stc{'w'};
hs = bar(eds,histc(ang(ind,4,7,2),eds),'histc'); hs.FaceColor = 'r'; hs.FaceAlpha = .5;















figure,plot()
figure,plot(ang(:,'bcom','hcom',2))


fxyz = xyz.copy;
fxyz.filter('ButFilter',3,1,'low');
fang = create(MTADang,Trial,fxyz);

name = 'head and body com'; label = 'hbcom'; key = 'c';
zv = MTADfet.encapsulate(Trial,...
                         [ang(:,'bcom','hcom',2)-fang(:,'bcom','hcom',2)],...
                         xyz.sampleRate,...
                         name,label,key);


dspec = struct('nFFT',2^9,'Fs',zv.sampleRate,...
               'WinLength',2^8,'nOverlap',2^8*.875,...
               'FreqRange',[1,30]);
[ys,fs,ts] = fet_spec(Trial,zv,'mtchglong',true,[],dspec,true);


figure,imagesc(ts,fs,ys(:,:,1,2)'),axis xy,colormap jet

figure,imagesc(ts,fs,log10(ys(:,:,1,1))'),axis xy,colormap jet
figure,imagesc(ts,fs,log10(ys(:,:,2,2))'),axis xy,colormap jet

eds = linspace(-.2,1.5,100);
figure,hold on
ind = Trial.stc{'s'};
hs = bar(eds,histc(ang(ind,1,3,2),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='c';
ind = Trial.stc{'p'};
hs = bar(eds,histc(ang(ind,1,3,2),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='y';
ind = Trial.stc{'w'};
hs = bar(eds,histc(ang(ind,1,3,2),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='b';
ind = Trial.stc{'m'};
hs = bar(eds,histc(ang(ind,1,3,2),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='m';
ind = Trial.stc{'n'};
hs = bar(eds,histc(ang(ind,1,3,2),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='g';
ind = Trial.stc{'r'};
hs = bar(eds,histc(ang(ind,1,3,2),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='r';




tfet = MTADxyz('data',bsxfun(@rdivide,ang_b_vel,10.^fang(:,3,4,2)),'sampleRate',fang.sampleRate);
tfet.filter('ButFilter',3,2,'low');
tfet.data = mean(tfet.data,2);
tfet.data(tfet.data<1e-6) = 1e-6;
tfet.data = log10(tfet.data);

% $$$ figure,plot(tfet);Lines(Trial.stc{'w'}(:),[],'b');Lines(Trial.stc{'n'}(:)+1,[],'g');



eds = linspace(-6,1,100);
figure,hold on
ind = Trial.stc{'s'};
hs = bar(eds,histc(tfet(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='c';
ind = Trial.stc{'p'};
hs = bar(eds,histc(tfet(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='y';
ind = Trial.stc{'w'};
hs = bar(eds,histc(tfet(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='b';
ind = Trial.stc{'m'};
hs = bar(eds,histc(tfet(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='m';
ind = Trial.stc{'n'};
hs = bar(eds,histc(tfet(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='g';
ind = Trial.stc{'r'};
hs = bar(eds,histc(tfet(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='r';



Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev2');
xyz = Trial.load('xyz');
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.4,'low');
figure,plot(diff(fxyz(:,1,3)));Lines(Trial.stc{'w'}(:),[],'b');Lines(Trial.stc{'n'}(:)+1,[],'g');




fxyz = xyz.copy;
fxyz.filter('ButFilter',3,1.5,'low');

bfet = Trial.xyz.copy;
tsh = 20;
afet = Trial.xyz.copy;
afet.data = circshift(fxyz.data,-tsh)-circshift(fxyz.data,tsh);
afet.data = reshape(afet.data,[],3);
afet.data = permute(bsxfun(@dot,permute(reshape(repmat(fxyz(:,4,:)-fxyz(:,1,:),[1,fxyz.size(2),1]),[],3),[2,1]),permute(afet.data,[2,1])),[2,1]);
afet.data = reshape(afet.data,[],fxyz.size(2));
bfet.data = mean(afet(:,1:4),2)./1000.*log10(var(afet(:,1:4),[],2));
bfet.filter('ButFilter',3,2.4,'low');

lbfet = bfet.copy;

lbfet.data = abs(log10(bfet.data(:,1)+10))

figure,plot(lbfet.data)
Lines(Trial.stc{'w'}(:),[],'b');Lines(Trial.stc{'n'}(:)+1,[],'g');


eds = linspace(0,3,100);
figure,hold on
ind = Trial.stc{'s'};
hs = bar(eds,histc(lbfet(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='c';
ind = Trial.stc{'p'};
hs = bar(eds,histc(lbfet(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='y';
ind = Trial.stc{'w'};
hs = bar(eds,histc(lbfet(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='b';
ind = Trial.stc{'m'};
hs = bar(eds,histc(lbfet(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='m';
ind = Trial.stc{'n'};
hs = bar(eds,histc(lbfet(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='g';
ind = Trial.stc{'r'};
hs = bar(eds,histc(lbfet(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='r';

% $$$ mwp = [];
% $$$ for per = Trial.stc{'w'}.data'
% $$$     mwp(end+ =maxlbfet(per(1):per(
    



name = 'hRL_shake'; label = 'hlrs'; key = 'h';
zv = MTADfet.encapsulate(Trial,...
                         circ_dist(circshift(fang(:,8,6,1),-1),fang(:,8,6,1)),...
                         xyz.sampleRate,...
                         name,label,key);


dspec = struct('nFFT',2^8,'Fs',zv.sampleRate,...
               'WinLength',2^7,'nOverlap',2^7*.875,...
               'FreqRange',[1,30]);
[ys,fs,ts] = fet_spec(Trial,zv,'mtchglong',true,[],dspec,true);

figure,imagesc(ts,fs,log10(ys(:,:,1,1))'),axis xy,colormap jet

figure,imagesc(ts,fs,ys(:,:,1,2)'),axis xy,colormap jet

,Lines(StcHL{'m',xyz.sampleRate}(:),[],'m');





f = 9;
figure;hold on;
eds =linspace(20,160,100),
ind = stcjg{'w'};
hs = bar(eds,histc(ofet(ind,f),eds),'histc');,hs.FaceAlpha=.5;hs.FaceColor='c';
ind = stced{'w'};
hs = bar(eds,histc(mfet(ind,f),eds),'histc');,hs.FaceAlpha=.5;hs.FaceColor='r';


f = 7;
figure;hold on;
eds =linspace(-pi/2,pi/2,100),
ind = stcjg{'a'};
hs = bar(eds,histc(ofet(ind,f),eds),'histc');,hs.FaceAlpha=.5;hs.FaceColor='c';
ind = stced{'a'};
hs = bar(eds,histc(mfet(ind,f),eds),'histc');,hs.FaceAlpha=.5;hs.FaceColor='r';





%% head stuff

Trial = MTATrial('Ed01-20140707');
Trial.load('stc','hand_labeled_rev1_Ed');
xyz = Trial.load('xyz');

xyz.addMarker('bcom',[.7,0,.7],{},...
    xyz.com(xyz.model.rb({'spine_lower','pelvis_root','spine_middle'})));
xyz.addMarker('scom',[.7,0,.7],{},...
    xyz.com(xyz.model.rb({'pelvis_root','spine_middle','spine_upper'})));
xyz.addMarker('hcom',[.7,0,.7],{},...
    xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));

fxyz = xyz.copy;
fxyz.filter('ButFilter',3,1);

xyz.addMarker('fsl',[.7,0,.7],{},fxyz(:,'spine_lower',:));
xyz.addMarker('fbcom',[.7,0,.7],{},fxyz(:,'bcom',:));
xyz.addMarker('fscom',[.7,0,.7],{},fxyz(:,'scom',:));
xyz.addMarker('fhcom',[.7,0,.7],{},fxyz(:,'hcom',:));


ang = create(MTADang,Trial,xyz);
figure,plot(circ_dist(ang(:,'fbcom','fscom',1),ang(:,'head_back','head_front',1)))

figure,plot(circ_dist(ang(:,'fbcom','fscom',1),ang(:,'fscom','hcom',1)))

Lines(Trial.stc{'w'}(:),[],'c');
Lines(Trial.stc{'r'}(:),[],'r');



proxy_ss = ss.copy;
ss = proxy_ss.copy;

ss.filter('ButFilter',3,10);

dss = sqrt(sum(diff(ss.data).^2,3));

figure,imagesc(log10(dss')),
caxis([-1,2]),colormap copper
Lines(Trial.stc{'w'}(:),[],'c');
Lines(Trial.stc{'r'}(:),[],'r');
Lines(Trial.stc{'m'}(:),[],'m');
Lines(Trial.stc{'n'}(:),[],'g');

figure,imagesc(diff(dss)'),axis xy
caxis([-.2,.2]),colormap copper
Lines(Trial.stc{'w'}(:)+.3,[],'c');
Lines(Trial.stc{'r'}(:)-.2,[],'r');
Lines(Trial.stc{'m'}(:)+.7,[],'m');
Lines(Trial.stc{'n'}(:)+.1,[],'g');


figure,plot(sum(dss,2)),
Lines(Trial.stc{'w'}(:),[],'c');
Lines(Trial.stc{'r'}(:),[],'r');
Lines(Trial.stc{'m'}(:),[],'m');
Lines(Trial.stc{'n'}(:),[],'g');


ss = proxy_ss.copy;
ss.filter('ButFilter',3,20);

dss = sqrt(sum(diff(ss.data).^2,3));

name = 'speed 3dss'; label = 'dss'; key = 'd';
zv = MTADfet.encapsulate(Trial,...
                         sum(dss,2),...
                         xyz.sampleRate,...
                         name,label,key);

dspec = struct('nFFT',2^8,'Fs',zv.sampleRate,...
               'WinLength',2^7,'nOverlap',2^7*.875,...
               'FreqRange',[1,10]);
[ys,fs,ts] = fet_spec(Trial,zv,'mtchglong',false,[],dspec,false);

figure,imagesc(ts,fs,log10(ys.data)'),axis xy, colormap jet

eds = linspace(-5,4,100);
figure,hold on;
ind = Trial.stc{'s'};
hs = bar(eds,histc(mean(log10(ys(ind,round(fs)==4)),2),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='c';
ind = Trial.stc{'p'};
hs = bar(eds,histc(mean(log10(ys(ind,round(fs)==4)),2),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='r';

s = 20;
figure,plot(abs(circ_dist(circshift(ang(:,1,4,1),-s), ...
                                 circshift(ang(:,1,4,1),s)))),Lines(Trial.stc{'n'}(:),[],'g');


fxyz = Trial.load('xyz');
fxyz.filter('ButFilter',3,0.5);
ang = create(MTADang,Trial,fxyz)
s = 15;
name = 'angspeed'; label = 'as'; key = 'a';
zv = MTADfet.encapsulate(Trial,...
                         abs(circ_dist(circshift(ang(:,1,5,1),-s),circshift(ang(:,1,5,1),s))),...
                         xyz.sampleRate,...
                         name,label,key);
zv.filter('ButFilter',3,1);
zv.data(zv.data<1e-5) = 1e-5;
zv.data = log10(zv.data);

eds = linspace(-5,4,100);
figure,hold on;
ind = Trial.stc{'p'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='c';
ind = Trial.stc{'n'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='g';
ind = Trial.stc{'p'};
hs = bar(eds,histc(zv(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='m';


try
    man = Trial.load('fet','lsppc');
catch err
    gen_fet_lsppc(Trial);    
    man = Trial.load('fet','lsppc');
end
man.filter('ButFilter',3,2,'low');

edx = linspace(-0.2,1,100);
edy = linspace(-5,.5,100);
figure,
subplot(2,2,1)
ind = Trial.stc{'a'};
hist2([man(ind),zv(ind)],edx,edy)
subplot(2,2,2)
ind = Trial.stc{'n'};
hist2([man(ind),zv(ind)],edx,edy)
subplot(2,2,3)
ind = Trial.stc{'p'};
hist2([man(ind),zv(ind)],edx,edy)
subplot(2,2,4)
ind = Trial.stc{'w'};
hist2([man(ind),zv(ind)],edx,edy)
ForAllSubplots('grid on;caxis([0,200])')



eds = linspace(-.3,pi/2,100);
figure,hold on
ind = Trial.stc{'a-r'};
hs = bar(eds,histc(circ_dist(ang(ind,1,2,2),ang(ind,2,3,2)),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'c';
ind = Trial.stc{'r'};
hs = bar(eds,histc(circ_dist(ang(ind,1,2,2),ang(ind,2,3,2)),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'r';



Trial = MTATrial('jg05-20120317');
Trial = MTATrial('jg05-20120317');
figure,plot(diff(fang(:,3,4,2))),Lines(Trial.stc{'r'}(:),[],'r');


zrf =Trial.xyz.copy;
zrf.data = diff(fang(:,3,4,2));


zrs = zrf.segs(Trial.stc{'r'}(:,1)-120,240);
ss = [];
for s = 1:size(zrs,2),
    [~,ss(end+1)] = max(convn(zrs(:,s),zrs(90:150,2),'valid'));
end
ss = ss-90;
zre = zrf.segs(Trial.stc{'r'}(:,2)-120,240);
se = [];
for s = 1:size(zrs,2),
    [~,se(end+1)] = max(convn(zre(:,s),zre(90:150,2),'valid'));    
end
se = se-90;


zrsc = zrf.segs(Trial.stc{'r'}(:,1)-ss'-120,240);
zrse = zrf.segs(Trial.stc{'r'}(:,2)-se'-120,240);

zct = nanmedian(zrsc(:,:),2);
zet = nanmedian(zrse(:,:),2);

rsc = convn(zrf.data,zct,'same');
rse = convn(zrf.data,zet,'same');

rtc = LocalMinimaN(-rsc,-0.0025,20);
rte = LocalMinimaN(-rse,-0.0025,20);

rtp = sort([rtc(:,1)-60;rte(:,1)+60]);

rtper = [rtp(1:end-1),rtp(2:end)];

figure,hold on
mh = [];
for p = rtper',
    mh(end+1) = median(xyz(p',7,3));
    plot(p(1):p(2),repmat(mh(end),1,1+p(2)-p(1)))
end
Lines(Trial.stc{'r'}(:),[],'r');




rper = Trial.stc{'r'};

cc = bsxfun(@minus,rtp(:,1)',rper(:,1));




%% Turning Feature
Trial = MTATrial('jg05-20120317');
xyz = Trial.load('xyz');
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,1.5,'low');
fang = create(MTADang,Trial,fxyz);

fxyz = xyz.copy;

fxyz.filter('ButFilter',3,2.4);
fvxy = fxyz.vel(1,[1,2]);


figure,
ind = Trial.stc{'a-m-r'};
hist2([ang(ind,1,2,2),abs(circ_dist(ang(ind,1,3,1),ang(ind,3,7,1)))],linspace(0.35,pi/2,100),0:0.01:pi/2);
caxis([0,200]);

dxyz = xyz.copy;
dxyz.data(:,1:4,:) = ss(:,[5,35,65,95],:);
dxyz.filter('ButFilter',3,20,'low');
dang = create(MTADang,Trial,dxyz);
%figure,plot(circ_dist(dang(:,1,2,2))
figure,plot(dang(:,1,2,2))
Lines(Trial.stc{'w'}(:),[],'b');
Lines(Trial.stc{'n'}(:),[],'g');
Lines(Trial.stc{'r'}(:),[],'r');
Lines(Trial.stc{'m'}(:),[],'m');



name = 'lower spine pitch'; label = 'SLpitch'; key = 'p';
zv = MTADfet.encapsulate(Trial,...
                         dang(:,1,2,2),...
                         xyz.sampleRate,...
                         name,label,key);

dspec = struct('nFFT',2^9,'Fs',zv.sampleRate,...
               'WinLength',2^8,'nOverlap',2^8*.875,...
               'FreqRange',[1,30]);
[ys,fs,ts] = fet_spec(Trial,zv,'mtchglong',true,[],dspec,true);


figure,imagesc(ts,fs,log10(ys.data)');
axis xy,caxis([-9,-6]),colormap jet
Lines(Trial.stc{'w',1}(:),[],'b');
Lines(Trial.stc{'n',1}(:),[],'g');
Lines(Trial.stc{'r',1}(:),[],'r');
Lines(Trial.stc{'m',1}(:),[],'m');




sh = 1;
sang = [circ_dist(circshift(fang(:,1,2,1),-sh),circshift(fang(:,1,2,1),sh)),...
        circ_dist(circshift(fang(:,1,3,1),-sh),circshift(fang(:,1,3,1),sh)),...
        circ_dist(circshift(fang(:,1,4,1),-sh),circshift(fang(:,1,4,1),sh)),...
        circ_dist(circshift(fang(:,1,5,1),-sh),circshift(fang(:,1,5,1),sh)),...        
        circ_dist(circshift(fang(:,1,7,1),-sh),circshift(fang(:,1,7,1),sh))...                
       ];
zav = Trial.xyz.copy;
zav.data =  log10(mean(abs(sang),2)./(var(sang,[],2)+1));

figure,plot(mean(abs(sang),2)./(var(sang,[],2)+1))
Lines(Trial.stc{'w'}(:),[],'b');
Lines(Trial.stc{'n'}(:),[],'g');

eds = linspace(-6,0,100);
figure,hold on
ind = Trial.stc{'a-n-w-r'};
hs = bar(eds,histc(zav(ind),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'c';
ind = Trial.stc{'n'};
hs = bar(eds,histc(zav(ind),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'r';
Lines([],0.0025,'k');

figure,plot(fang(:,1,7,3))
Lines(Trial.stc{'w'}(:),[],'b');
Lines(Trial.stc{'n'}(:),[],'g');
Lines(Trial.stc{'r'}(:),[],'r');



edx = linspace(-3,2,100);
edy = linspace(-6,0,100);
figure,
ind = Trial.stc{'p'};
hist2([log10(fvxy(ind)),zav(ind)],edx,edy)
caxis([0,300]);




ss = Trial.load('fet','3dss');
fs = ss.copy;
fs.filter('ButFilter',3,.4);
sd = sqrt(sum((fs.data-circshift(fs.data,-1,2)).^2,3));


sv = Trial.xyz.copy;
sv.data = sum(sd(:,2:end-1),2)./sd(:,end);

figure,plot(sv.data)
Lines(Trial.stc{'w'}(:),[],'b');
Lines(Trial.stc{'n'}(:),[],'g');
Lines(Trial.stc{'r'}(:),[],'r');
Lines(Trial.stc{'m'}(:),[],'m');


key = 'r';
eds = linspace(-0.1,.3,100);
figure,hold on
ind = Trial.stc{['a-' key]};
hs = bar(eds,histc(log10(sv(ind)),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'c';
ind = Trial.stc{key};
hs = bar(eds,histc(log10(sv(ind)),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'r';


edx = linspace(-.1,.3,100);
edy = linspace(1.4,2.5,100);
figure,
ind = Trial.stc{'r'};
hist2([log10(sv(ind)),log10(fang(ind,1,4,3).*cos(fang(ind,1,4,2)))],edx,edy)
caxis([0,300]);


Trial = MTATrial('jg05-20120317');
xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,30);

wfet = [0;ButFilter(diff(xyz(:,1,3)),3,[2,10]/(0.5*xyz.sampleRate),'bandpass')];

figure,
plot(diff(wfet))
Lines(Trial.stc{'w'}(:),[],'r');


name = 'lower spine Z acc'; label = 'lsza'; key = 'a';
zv = MTADfet.encapsulate(Trial,...
                         diff(wfet),...
                         xyz.sampleRate,...
                         name,label,key);

dspec = struct('nFFT',2^7,'Fs',zv.sampleRate,...
               'WinLength',2^6,'nOverlap',2^6*.875,...
               'FreqRange',[1,30]);
[ys,fs,ts] = fet_spec(Trial,zv,'mtchglong',false,[],dspec,true);

figure,imagesc(ts,fs,log10(ys.data)');
axis xy,caxis([-9,-1]),colormap jet


key = 'w';
eds = linspace(-9,0,100);
figure,hold on
ind = Trial.stc{['a-n-r-' key]};
hs = bar(eds,histc(log10(ys(ind,5)),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'c';
ind = Trial.stc{key};
hs = bar(eds,histc(log10(ys(ind,5)),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'r';


key = 'w';
eds = linspace(-3,2,100);
figure,hold on
ind = Trial.stc{['a-n-r-' key]};
hs = bar(eds,histc(log10(fvxy(ind)),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'c';
ind = Trial.stc{key};
hs = bar(eds,histc(log10(fvxy(ind)),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'r';


vp = [];
for s = Trial.stc{'w',fvxy.sampleRate}.data',
    vp(end+1) = log10(median(fvxy(s')));
end

yp = [];
for s = Trial.stc{'w',ys.sampleRate}.data',
    yp(end+1) = log10(median(ys(s',5)));
end

dp = [];
for s = Trial.stc{'w'}.data',
    dp(end+1) = log10(diff(s'));
end

figure,
hist2([dp,yp],linspace(0.2,2,20),linspace(-5,-2,20))


xyz = Trial.load('xyz');
fxyz = xyz.copy;
fxyz.resample(ys);
fxyz.filter('ButFilter',3,2.4);
fvxy = fxyz.vel(1,[1,2]);

edx = linspace(-3,2,100);
edy = linspace(-8,0,100);
figure,
ind = Trial.stc{'w'};
hist2([log10(fvxy(ind)),log10(ys(ind,5))],edx,edy)
caxis([0,30]);



name = 'zv and wag pw'; label = 'zvwp'; key = 'w';
zfet = MTADfet.encapsulate(Trial,...
                         log10(nanmedian(ys(:,1:fs<7,1,1),2)),...
                         ys.sampleRate,...
                         name,label,key);
zfet.data(~nniz(zfet),:)=-20;
zfet.resample(bfet);



sl = SessionList('hand_labeled_jg');
model = ['MTAC_BATCH-fet_tsne_rev13_SR_12_NORM_1'...
         '_REF_Ed03-20140625.cof.all_STC_hand_labeled_rev1_Ed'...
         '_NN_100_NI_200_NN_multiPN_RAND_WSBNT-wrnpms'];

vp = [];
dp = [];
wd = [];
for t = sl,
    
    Trial = MTATrial.validate(t);
    %Trial.load('stc',model);
    xyz = Trial.load('xyz');
    fxyz = xyz.copy;
    fxyz.filter('ButFilter',3,2.4);
    fvxy = fxyz.vel(1,[1,2]);
    
    for s = Trial.stc{'w'}.data',
        vp(end+1) = log10(median(fvxy(s')));
        wd(end+1) = log10(sqrt(sum((fxyz(s(2),1,[1,2])-fxyz(s(1),1,[1,2])).^2,3)));
        dp(end+1) = log10(diff(s'));
    end
    
end


edx = linspace(0,3,30);
edy = linspace(0,2,30);
figure(1),
subplot(211)
plot(dp,vp,'.')
title({'Session Set: hand\_labeled\_jg',...
       'Labeling Method: Hand Labeled',...
       ...%'Labeling Method: Neural Network WSBNT',...       
       'Walk: mean speed vs period duration'})
ylabel('mean(log10(lowerSpine Speed))');xlabel('log10(samples)')
xlim([0,3]);ylim([0,2])
subplot(212)
hist2([dp',vp'],edx,edy);
ylabel('mean(log10(lowerSpine Speed))');xlabel('log10(samples)')
xlim([0,3]);ylim([0,2])



edx = linspace(0,3,30);
edy = linspace(0,3,30);
figure(1),
subplot(211); plot(dp,wd,'.')
title({'Session Set: hand\_labeled\_jg',...
       'Labeling Method: Hand Labeled',...
       ...%'Labeling Method: Neural Network WSBNT',...
       'Walk: total distance vs period duration'})
xlim([0,3]);ylim([0,3])
ylabel('total log10(dist)');xlabel('log10(samples)')
subplot(212); hist2([dp',wd'],edx,edy);
xlim([0,3]);ylim([0,3])
ylabel('total log10(dist)');xlabel('log10(samples)')





zv = afet.copy;
%zv.data = mean(afet(:,1:3),2).*(1+var(afet(:,1),[],2));
zv.data = afet(:,1);
zv.data = log10(abs(zv(:,1))).*sign(zv(:,1));

key = 'w';
eds = linspace(-4,5,100);
figure,hold on
ind = Trial.stc{['a-n-r-' key]};
hs = bar(eds,histc(zv(ind),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'c';
ind = Trial.stc{key};
hs = bar(eds,histc(zv(ind),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'r';

