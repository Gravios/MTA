
Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev2');
xyz = Trial.load('xyz');





                 
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
fxyz.filter('ButFilter',3,1,'low');

ang = create(MTADang,Trial,xyz);
tsh = 1;
afet = Trial.xyz.copy;
bfet = Trial.xyz.copy;
afet.data = circshift(xyz(:,:,[1,2]),-tsh)-circshift(fxyz(:,:,[1,2]),tsh);
afet.data = reshape(afet.data,[],2);
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
vfxy.filter('ButFilter',3,1,'low');
vfxy = vfxy.vel([],[1,2]);

fvxy = xyz.copy;
fvxy = fvxy.vel([],[1,2]);
fvxy.filter('ButFilter',3,1,'low');

hx = vfxy; hxi = 5; edx = linspace(-3,2,100);
hy = fvxy; hyi = 5; edy = linspace(-1,2,100);

clev = [0,400];
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
ind = Trial.stc{'m'};
hist2(log10([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
caxis(clev)

figure,
hxi = 7;
hyi = 7;
ind = Trial.stc{'m'};
hist2(log10([hx(ind,hxi),hy(ind,hyi)]),edx,edy)
caxis(clev)

