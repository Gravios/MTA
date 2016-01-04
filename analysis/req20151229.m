Trial = MTATrial('Ed03-20140624');
Trial.load('stc','hand_labeled_rev1');
Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev2');
xyz = Trial.load('xyz');



% COM Body Center of Mass
rbb = xyz.model.rb({'spine_lower','pelvis_root','spine_middle'});
xyz.addMarker('bcom',[.7,0,.7],{{'head_back','head_front',[0,0,255]}},xyz.com(rbb));


% FXYZ Filtered Marker Positions {Low Pass 2.5 Hz}
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,.5,'low');

xyz.addMarker('fbl',[.7,0,.7],{{'fbl','pelvis_root',[0,0,255]}},fxyz(:,1,:));

fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.5,'low');


% ANG InterMarker Spherical Coordinates
ang = create(MTADang,Trial,xyz);

%figure,plot(circ_dist(circshift(ang(:,10,12,1),-2),circshift(ang(:,10,12,1),2)))

% FANG Filtered Intermarker angles 
fang = create(MTADang,Trial,fxyz);


                 
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,.1,'high');                 
figure,plot(diff(fxyz(:,[1,2,3,4],3)))
Lines(Trial.stc{'n'}(:),[],'g');Lines(Trial.stc{'w'}(:),[],'c');

fvfz = xyz.copy;
fvfz.data = diff(fxyz(:,[1,2,3,4],3));
fvfz.filter('ButFilter',3,2,'low');                 
figure,plot(fvfz.data)
Lines(Trial.stc{'n'}(:),[],'g');Lines(Trial.stc{'w'}(:),[],'c');

figure,plot(sqrt(sum([circshift(xyz(:,1,[1,2]),-30)-circshift(xyz(:,1,[1,2]),30)].^2,3)))
Lines(Trial.stc{'n'}(:),[],'g');Lines(Trial.stc{'w'}(:),[],'c');

xyz = Trial.load('xyz');

fws = 1:.2:5;
sind = Trial.stc{'w'};
nind = Trial.stc{'p'};
vdps = [];
nnfor f=fws;
    vfx = xyz.copy;
    vfx.data = vfx(:,1,[1,2]);
    vfx.filter('ButFilter',3,f,'low');
    vfx = vfx.vel(1,[1,2]);
    vdps(end+1) = (mean(vfx(sind))-mean(vfx(nind)))/sqrt(.5*(var(vfx(sind))+var(vfx(nind))));
end


%shift = 30;
shifts = 5:2:80;
sind = Trial.stc{'w'};
nind = Trial.stc{'p'};
dps = [];
for shift = shifts,
    vcs = xyz.copy;
    vcs.data = sqrt(sum([circshift(xyz(:,1,[1,2]),-shift)-circshift(xyz(:,1,[1,2]),shift)].^2,3));
    dps(end+1) = (mean(vcs(sind))-mean(vcs(nind)))/sqrt(.5*(var(vcs(sind))+var(vcs(nind))));
end

figure,plot(shifts,dps)

shift = 45;
vcs = xyz.copy;
vcs.data = sqrt(sum([circshift(xyz(:,:,[1,2]),-shift)-circshift(xyz(:,:,[1,2]),shift)].^2,3));


%figure,plot(log10(vcs.data))
%figure,plot(log10(vfx.data))

hx = bfet; hxi = 1; edx = linspace(-3,2,50);
hy = man; hyi = 1; edy = linspace(-.2,1,50);

hx = vcs; hxi = 1; edx = linspace(-3,3,100);
hy = vcs; hyi = 7; edy = linspace(-3,3,100);

hx = vfx; hxi = 1; edx = linspace(-1,2,100);
hy = vcs; hyi = 1; edy = linspace(-3,3,100);

clev = [0,100];
figure,
subplot(221)
ind = Trial.stc{'a-w-n'};
hist2([log10(hx(ind,hxi)),hy(ind,hyi)],edx,edy)
caxis(clev)
subplot(222)
ind = Trial.stc{'w'};
hist2([log10(hx(ind,hxi)),hy(ind,hyi)],edx,edy)
caxis(clev)
subplot(223)
ind = Trial.stc{'a'};
hist2([log10(hx(ind,hxi)),hy(ind,hyi)],edx,edy)
caxis(clev)
subplot(224)
ind = Trial.stc{'n'};
hist2([log10(hx(ind,hxi)),hy(ind,hyi)],edx,edy)
caxis(clev)

hx = bfet; hxi = 1; edx = linspace(-3,2,50);
hy = man; hyi = 1; edy = linspace(-.2,1,50);

hx = bfet; hxi = 1; edx = linspace(-3,2,50);
hy = fvelxy; hyi = 1; edy = linspace(-3,2,50);

hx = vcs; hxi = 1; edx = linspace(-3,3,100);
hy = vcs; hyi = 7; edy = linspace(-3,3,100);

hx = vfx; hxi = 1; edx = linspace(-1,2,100);
hy = vcs; hyi = 1; edy = linspace(-3,3,100);

clev = [0,100];
figure,
subplot(221)
ind = Trial.stc{'a-w-n'};
hist2(log10([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
caxis(clev)
subplot(222)
ind = Trial.stc{'w'};
hist2(log10([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
caxis(clev)
subplot(223)
ind = Trial.stc{'a'};
hist2(log10([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
caxis(clev)
subplot(224)
ind = Trial.stc{'n'};
hist2(log10([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
caxis(clev)

f = 2.4;
vfx = xyz.copy;
vfx.data = vfx(:,1,[1,2]);
vfx.filter('ButFilter',3,f,'low');
vfx = vfx.vel(1,[1,2]);

shift = 45;
vcs = xyz.copy;
vcs.data = sqrt(sum([circshift(xyz(:,1,[1,2]),-shift)-circshift(xyz(:,1,[1,2]),shift)].^2,3));
vcs.resample(newSampleRate);

fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.4,'low');
fvelxy = fxyz.vel([],[1,2]);
fvelxy.resample(newSampleRate);
%fvelxy.filter('ButFilter',3,2.4,'low');

fvelxy = xyz.vel([],[1,2]);
fvelxy.resample(newSampleRate);
fvelxy.filter('ButFilter',3,2.4,'low');



hvar = vcs;
eds = linspace(-3,3,100);

hvar = vfx;
hvar = fvelxy;
eds = linspace(-3,2,100);
eds = linspace(-.5,2,100);

eds = linspace(-2,2.4,100);
eds = linspace(-2,3,100);
hind = 1;

figure,hold on
ind = Trial.stc{'p'};
hs = bar(eds,histc(log10(hvar(ind,hind)),eds),'histc');
hs.FaceColor = 'c';
hs.FaceAlpha = .4;

ind = Trial.stc{'s'};
hs = bar(eds,histc(log10(hvar(ind,hind)),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;

ind = Trial.stc{'w'};
hs = bar(eds,histc(log10(hvar(ind,hind)),eds),'histc');
hs.FaceColor = 'g';
hs.FaceAlpha = .4;









%%%% Turning Feature
Trial = MTATrial('jg05-20120317');
xyz = Trial.load('xyz');
ang = create(MTADang,Trial,xyz);
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,0.1,'low');
fang = create(MTADang,Trial,fxyz);
figure,
hold on,plot(circ_dist(ang(:,1,4,1),fang(:,1,4,1)))
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,0.3,'low');
fang = create(MTADang,Trial,fxyz);
hold on,plot(circ_dist(ang(:,1,4,1),fang(:,1,4,1)))
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,0.5,'low');
fang = create(MTADang,Trial,fxyz);
hold on,plot(circ_dist(ang(:,1,4,1),fang(:,1,4,1)))
Lines(Trial.stc{'n'}(:),[],'g');Lines(Trial.stc{'w'}(:),[],'c');

figure,plot(mean([circ_dist(ang(:,1,3,1),fang(:,1,4,1)),circ_dist(ang(:,1,4,1),fang(:,1,4,1))],2))
Lines(Trial.stc{'n'}(:),[],'g');Lines(Trial.stc{'w'}(:),[],'c');



fxyz = xyz.copy;
fxyz.filter('ButFilter',3,0.2,'low');
fang = create(MTADang,Trial,fxyz);
figure,plot(abs(diff(circ_dist(ang(:,1,4,1),fang(:,1,4,1)))))
Lines(Trial.stc{'n'}(:),[],'g');Lines(Trial.stc{'w'}(:),[],'c');

%%%% Walking Feature
Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev2');
xyz = Trial.load('xyz');

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

figure,plot(bfet.data)
Lines(Trial.stc{'w'}(:)+1,[],'c');
Lines(Trial.stc{'n'}(:),[],'g');

edx = linspace(-50,200,200);
figure,hold on
ind = Trial.stc{'w'}; hs = bar(edx,histc(bfet(ind),edx),'histc'); hs.FaceColor = 'c'; hs.FaceAlpha = .4;
ind = Trial.stc{'n'}; hs = bar(edx,histc(bfet(ind),edx),'histc'); hs.FaceColor = 'g'; hs.FaceAlpha = .4;
ind = Trial.stc{'p'}; hs = bar(edx,histc(bfet(ind),edx),'histc'); hs.FaceColor = 'r'; hs.FaceAlpha = .4;




aft = mat2cell(afet.data,size(afet,1),[1,1]);
[afet.data,bfet.data] = cart2pol(aft{:});
% $$$ afet.data = reshape(afet.data,[],xyz.size(2));
bfet.data = reshape(bfet.data,[],xyz.size(2));
% $$$ m = MTADxyz('data',circ_dist(afet(:,1),ang(:,1,4,1)),'sampleRate',Trial.xyz.sampleRate);
% $$$ m.data = circ_dist(circshift(m.data,-5),circshift(m.data,5));
% $$$ m.data = [diff(m.data);0];
%bfet.resample(newSampleRate);

%bfet.filter('ButFilter',3,5,'low')



name = 'body_traj_magnitude'; label = 'bfet'; key = 'b';
fet = MTADfet.encapsulate(Trial,...
                          bfet.data,...
                          bfet.sampleRate,...
                          name,label,key);



figure,hold on
ind = Trial.stc{'p'};
hs = bar(eds,histc(log10(hvar(ind,hind)),eds),'histc');
hs.FaceColor = 'c';
hs.FaceAlpha = .4;

ind = Trial.stc{'s'};
hs = bar(eds,histc(log10(hvar(ind,hind)),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;

ind = Trial.stc{'w'};
hs = bar(eds,histc(log10(hvar(ind,hind)),eds),'histc');
hs.FaceColor = 'g';
hs.FaceAlpha = .4;




%% Grooming stuff

name = 'center of mass of head markers'; label = 'hcom'; key = 'h';
hcom = MTADfet.encapsulate(Trial,...
                           xyz.com(xyz.model.rb(xyz.model.find('head_[^t]'))),...
                           xyz.sampleRate,...
                           name,label,key);

fhcom = hcom.copy;
fhcom.filter('ButFilter',3,.2,'low');

figure,plot(fhcom(:,1,3))
hold on,plot(hcom(:,1,3))

figure,plot(hcom(:,1,3)-fhcom(:,1,3));
Lines(Trial.stc{'m'}(:),[],'m');

figure,plot(sq(hcom.data-fhcom.data));
Lines(Trial.stc{'m'}(:),[],'m');

vfxy = xyz.copy;
vfxy.filter('ButFilter',3,0.1,'low');
vfxy = vfxy.vel([],[1,2]);

fvxy = xyz.copy;
fvxy.filter('ButFilter',3,2.5,'low');
fvxy = fvxy.vel([],[1,2]);

%fvxy = xyz.copy;
%fvxy = fvxy.vel([],[1,2]);
%fvxy.filter('ButFilter',3,.1,'low');

hx = vfxy; hxi = 7; edx = linspace(-3,2,100);
hy = fvxy; hyi = 7; edy = linspace(-3,2,100);

clev = [0,400];
figure,
subplot(221)
ind = Trial.stc{'a-p-m'};
hist2(log10([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
caxis(clev)
subplot(222)
ind = Trial.stc{'p'};
hist2(log10([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
caxis(clev)
subplot(223)
ind = Trial.stc{'a'};
hist2(log10([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
caxis(clev)
subplot(224)
ind = Trial.stc{'m'};
hist2(log10([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
caxis(clev)



figure,hold on
ind = Trial.stc{'p'};
hs = bar(edx,histc(log10(hx(ind,hxi)),edx),'histc');
hs.FaceColor = 'c';
hs.FaceAlpha = .4;

ind = Trial.stc{'s'};
hs = bar(edx,histc(log10(hx(ind,hxi)),edx),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;

ind = Trial.stc{'m'};
hs = bar(edx,histc(log10(hx(ind,hxi)),edx),'histc');
hs.FaceColor = 'g';
hs.FaceAlpha = .4;



figure,
hxi = 7;
hyi = 7;
ind = Trial.stc{'m'};
hist2(log10([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
caxis(clev)


cfxyz = xyz.copy;
cfxyz.data = [hcom.data,fhcom.data];
cang = create(MTADang,Trial,cfxyz);

figure,plot(circ_dist(circshift(cang(:,1,2,2),-1),circshift(cang(:,1,2,2),1)))
Lines(Trial.stc{'m'}(:),[],'m');
% Stop what you are doing this isn't what you are suppose to be
% working on.

name = 'body bend xy', label='scurv';key='c';
gfet = MTADfet.encapsulate(Trial,...
                           10.^abs(circ_dist(ang(:,1,3,1),ang(:,1,7,1))),...
                           ang.sampleRate,...
                           name,label,key);

name = 'BLBU dist', label='hunch';key='h';
hfet = MTADfet.encapsulate(Trial,...
                           10.^abs(cos(ang(:,1,4,2)).*ang(:,1,4,3)),...
                           ang.sampleRate,...
                           name,label,key);




hx = hfet; hxi = 1; edx = linspace(50,160,100);
hy = gfet; hyi = 1; edy = linspace(0,2,100);


clev = [0,100];
figure,
subplot(231)
ind = Trial.stc{'a-p-m'};
hist2(log10([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
caxis(clev)
subplot(232)
ind = Trial.stc{'p'};
hist2(log10([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
caxis(clev)
subplot(233)
ind = Trial.stc{'w'};
hist2(log10([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
caxis(clev)
subplot(234)
ind = Trial.stc{'a'};
hist2(log10([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
caxis(clev)
subplot(235)
ind = Trial.stc{'m'};
hist2(log10([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
caxis(clev)
subplot(236)
ind = Trial.stc{'r'};
hist2(log10([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
caxis(clev)



figure,hold on
ind = Trial.stc{'p'};
hs = bar(edx,histc(hx(ind,hxi),edx),'histc');
hs.FaceColor = 'c';
hs.FaceAlpha = .4;

ind = Trial.stc{'s'};
hs = bar(edx,histc(hx(ind,hxi),edx),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;

ind = Trial.stc{'m'};
hs = bar(edx,histc(hx(ind,hxi),edx),'histc');
hs.FaceColor = 'g';
hs.FaceAlpha = .4;



name = 'mean_av'; label = 'angvel'; key = 'v';
mbang = MTADfet.encapsulate(Trial,...
                            circ_dist(circ_mean([ang(:,11,2,1),...
                                                 ang(:,11,3,1),...
                                                 ang(:,11,4,1)],...
                                                [],2),...
                                      ang(:,1,10,1)),...
                            ang.sampleRate,...
                            name,label,key);
%mbang.filter('ButFilter',3,[.1,5],'bandpass');
mbang.filter('ButFilter',3,10,'low');

mbang.data = diff(mbang.data);
mbang.filter('ButFilter',3,0.2,'high');
figure,plot(mbang.data)
Lines(Trial.stc{'m'}(:),[],'m');
Lines(Trial.stc{'w'}(:),[],'c');
Lines(Trial.stc{'n'}(:)+1,[],'g');

hold on,
figure,plot(circshift(var(mbang.segs(1:mbang.size(1),90,0))'.*100,45))

vbang = mbang.copy;
vbang.data = circshift(var(mbang.segs(1:mbang.size(1),80,0))'.*100,40);
vbang.data = log10(vbang.data);

cbang = mbang.copy;
cbang.data = ang(:,1,11,3);
cbang.filter('ButFilter',3,2.4,'low');
cbang.data(cbang.data<1e-4) = 1e-4;
cbang.data = log10(cbang.data);

hx = vbang; hxi= 1; edx = linspace(-7,2,100);
figure,hold on
ind = Trial.stc{'w'};
hs = bar(edx,histc(hx(ind,hxi),edx),'histc');
hs.FaceColor = 'c';
hs.FaceAlpha = .4;

ind = Trial.stc{'s'};
hs = bar(edx,histc(hx(ind,hxi),edx),'histc');
hs.FaceColor = 'm';
hs.FaceAlpha = .4;

ind = Trial.stc{'p'};
hs = bar(edx,histc(hx(ind,hxi),edx),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;

ind = Trial.stc{'n'};
hs = bar(edx,histc(hx(ind,hxi),edx),'histc');
hs.FaceColor = 'g';
hs.FaceAlpha = .4;




hx = ys; hxi = 2; edx = linspace(-10,-0,100);
hy = fvxy; hyi = 1; edy = linspace(-3,2,100);

figure,
subplot(231)
ind = Trial.stc{'a-w-n'};
hist2(log10([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
subplot(232)
ind = Trial.stc{'p'};
hist2(log10([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
subplot(233)
ind = Trial.stc{'w'};
hist2(log10([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
subplot(234)
ind = Trial.stc{'a'};
hist2(log10([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
subplot(235)
ind = Trial.stc{'n'};
hist2(log10([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
subplot(236)
ind = Trial.stc{'r'};
hist2(log10([hx(ind,hxi),hy(ind,hyi)]),edx,edy)

ForAllSubplots(['caxis([0,120]),grid on'])


[ys,fs,ts] = fet_spec(Trial,mbang,'mtchglong',true);

ys.data = ys(:,1:20);
ys.resample(vbang);

figure,imagesc(ts,fs,log10(ys.data)'),axis xy;colormap jet


%% Walk z axis vertical spectrum

%Z dynamics
fxyz = xyz.copy;

name = 'mean_av'; label = 'angvel'; key = 'v';
mbang = MTADfet.encapsulate(Trial,...
                            circ_dist(circ_mean([ang(:,11,2,1),...
                                                 ang(:,11,3,1),...
                                                 ang(:,11,4,1)],...
                                                [],2),...
                                      ang(:,1,10,1)),...
                            ang.sampleRate,...
                            name,label,key);
%mbang.filter('ButFilter',3,[.1,5],'bandpass');
mbang.filter('ButFilter',3,10,'low');
mbang.data = diff(mbang.data);

vbang = mbang.copy;
vbang.data = circshift(var(mbang.segs(1:mbang.size(1),60,0))'.*100,30);

name = 'lower spine Z speed'; label = 'lszs'; key = 'z';
zv = MTADfet.encapsulate(Trial,...
                         [diff(xyz(:,1,3)),vbang.data],...
                         xyz.sampleRate,...
                         name,label,key);

%figure,plot(zv.data)

dspec = struct('nFFT',2^7,'Fs',zv.sampleRate,...
               'WinLength',2^6,'nOverlap',2^6*.875,...
               'FreqRange',[1,30]);
[ys,fs,ts] = fet_spec(Trial,zv,'mtchglong',true,[],dspec);


figure,sp = [];
sp(1)=subplot(211);imagesc(ts,fs,log10(ys(:,:,1,1))'),axis xy;colormap jet
sp(2)=subplot(212);imagesc(ts,fs,log10(ys(:,:,2,2))'),axis xy;colormap jet
linkaxes(sp,'xy');
ForAllSubplots(['grid on'])



name = 'zv and wag pw'; label = 'zvwp'; key = 'w';
zvw = MTADfet.encapsulate(Trial,...
                         log10([nanmedian(ys(:,1:fs<7,1,1),2),nanmedian(ys(:,1:fs<7,2,2),2)]),...
                         ys.sampleRate,...
                         name,label,key);
zvw.data(~nniz(zvw),:)=-20;
zvw.resample(xyz);

lfvxy = copy(fvxy);
lfvxy.data = log10(lfvxy.data);
hx = lfvxy; hxi = 1; edx = linspace(-3,2,100);
hy = pzvw; hyi = 1; edy = linspace(6,20,100);

figure,
subplot(231)
ind = Trial.stc{'a-w-n'};
hist2(([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
subplot(232)
ind = Trial.stc{'p'};
hist2(([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
subplot(233)
ind = Trial.stc{'w'};
hist2(([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
subplot(234)
ind = Trial.stc{'a'};
hist2(([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
subplot(235)
ind = Trial.stc{'m'};
hist2(([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
subplot(236)
ind = Trial.stc{'n'};
hist2(([hx(ind,hxi),hy(ind,hyi)]),edx,edy)

ForAllSubplots(['caxis([0,120]),grid on'])


[U,S,V] = svd(cov(zvw(Trial.stc{'a'},:)));
pzvw = zvw.copy;
pzvw.data =zvw.data* V(:,1);

figure,plot(zvw.data)
Lines(Trial.stc{'w'}(:),[],'c');


hx = pzvw; hxi= 1; edx = linspace(6,20,100);
figure,hold on
ind = Trial.stc{'w'};
hs = bar(edx,histc(hx(ind,hxi),edx),'histc');
hs.FaceColor = 'c';
hs.FaceAlpha = .4;

ind = Trial.stc{'s'};
hs = bar(edx,histc(hx(ind,hxi),edx),'histc');
hs.FaceColor = 'm';
hs.FaceAlpha = .4;

ind = Trial.stc{'p'};
hs = bar(edx,histc(hx(ind,hxi),edx),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;

ind = Trial.stc{'n'};
hs = bar(edx,histc(hx(ind,hxi),edx),'histc');
hs.FaceColor = 'g';
hs.FaceAlpha = .4;



mpw = [];
for ind = Trial.stc{'w'}.data',
    mpw(end+1) = mean(pzvw(ind'));
end


mpp = [];
for ind = Trial.stc{'p'}.data',
    mpp(end+1) = mean(pzvw(ind'));
end

figure,hold on
hs = bar(edx,histc(mpw,edx),'histc'); hs.FaceColor = 'c'; hs.FaceAlpha = .4;
hs = bar(edx,histc(mpp,edx),'histc'); hs.FaceColor = 'r'; hs.FaceAlpha = .4;


tsh = 40;
txyz = xyz.copy;
txyz.data = circshift(xyz.data,-tsh)-circshift(xyz.data,tsh);
x = mat2cell(sq(txyz(:,1,:)),txyz.size(1),[1,1,1]);
tang = cart2sph(x{:});
tt = xyz.copy;
tt.data = unwrap(circ_dist(ang(:,1,4,1),tang));
tt.filter('ButFilter',3,[.01,30],'bandpass');

figure,plot(diff(tt.data))


figure,plot(diff())

