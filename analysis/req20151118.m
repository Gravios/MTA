% req20151118.m
%
% #DESC
%
% Examination of rats' casting motions during exploration
%
%  Goals:
%    
%    1. Identify types of head casting
%    2. Characterize phiysical characteristics of each casting type
%    3. Determine relationship between casting and sniffing
%


Trial = MTATrial('Ed01-20140707');
Trial = MTATrial('jg05-20120317');


xyz  = Trial.load('xyz');
% create a ridgid body model
rb = Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'});
% find the center of mass of the model
hcom = xyz.com(rb);
% add coordinates of the model's center of mass to the xyz object
xyz.addMarker('hcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},hcom);
xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},ButFilter(hcom,3,[2]./(Trial.xyz.sampleRate/2),'low'));
% remove high frequency noise 
%xyz.filter('ButFilter',3,55);
xyz.data = xyz(:,:,[1,2,3]);




ang = create(MTADang,Trial,xyz)
%ang.resample(30);

ncp = Trial.load('lfp',2)
ncp.resample(xyz);

vxy = xyz.vel([],[1,2]);
vxy.filter('ButFilter',3,2.4,'low');
%vxy.resample(ang);

ms=[1,3];
mo=[2,4];
shift = 0;

cang = MTADang('data',circ_dist(circshift(ang(:,ms(1),ms(2),1),-shift),circshift(ang(:,mo(1),mo(2),1),shift)),'sampleRate',ang.sampleRate);
cang.filter('ButFilter',3,2,'high');
cang.data = sq(sum(sqrt(cang.segs(1:size(cang,1),90,0).^2)))';


[rhm,fs,ts]  = fet_rhm(Trial,[],'mtchglong',false);
rhm.resample(ang);

ind = Trial.stc{'w'};
figure, hist2(log10([cang(ind),rhm(ind,40)]),linspace(-1.5,1.4,75),linspace(-7,-1,75))

ind = Trial.stc{'w'};
figure, hist2(log10([cang(ind),vxy(ind,1)]),linspace(-1.5,1.4,75),linspace(-1,2,75))


ind = Trial.stc{'w'};
figure, hist2([log10([vxy(ind,1)]),ang(ind,5,7,2)],linspace(-1,2,75),linspace(-1,1.4,75))


figure,plot(circ_dist(circshift(ang(:,4,10,1),-1),circshift(ang(:,4,10,1),1)),clip(ang(:,6,10,3),-40,40),'.')

ind = Trial.stc{'p'};
figure, hist2(log10([cang(ind),rhm(ind,40,2,2)]),linspace(-1.5,2,100),linspace(-8,-1,100))




cang = MTADang('data',diff(ang(:,8,10,3).*ang(:,8,10,3)),'sampleRate',ang.sampleRate);
cang.filter('ButFilter',3,5,'low');
cang.data = sq(sum(sqrt(cang.segs(1:size(cang,1),90,0).^2)))';
cang.resample(30);
ind = Trial.stc{'a'};

figure, 
eds = {linspace(0,4,75),linspace(-1,2,75)};
ind = Trial.stc{'a'};
subplot(121)
hist2(log10([cang(ind),vxy(ind,1)]),eds{1},eds{2})
ind = Trial.stc{'w'};
subplot(122)
hist2(log10([cang(ind),vxy(ind,1)]),eds{1},eds{2})



ang.filter('ButFilter',3,[2,20],'bandpass');
figure,plot(circ_dist(ang(1:end-1,11,10,1),ang(2:end,11,10,1)))
figure,plot(diff(circ_dist(ang(1:end-1,5,7,2),ang(2:end,5,7,2)).*10))
%figure,plot(diff(circ_dist(ang(1:end-1,5,10,2),ang(2:end,5,11,2)).*10))


cang = MTADang('data',ang(:,5,10,3),'sampleRate',ang.sampleRate);
cang.filter('ButFilter',3,10,'low');
bang = MTADang('data',ang(:,5,7,2),'sampleRate',ang.sampleRate);
bang.filter('ButFilter',3,5,'high');

figure,plot(diff(cang.data-30).*10)
hold on,plot(diff(bang.data).*10)
hold on,plot(diff(ncp.data)/5000)

figure, plot(ang(:,1,3,2)-ang(:,1,2,2))


figure,plot(ang(:,1,3,2))

wins = 80;
bang = MTADang('data',ang(:,1,4,2),'sampleRate',ang.sampleRate);
bang.filter('ButFilter',3,[4,10],'bandpass');
bang.data = sq(sum(sqrt(bang.segs(1:size(bang,1),wins,0).^2)))';
bang.data = circshift(bang.data,-round(wins/2));

figure,plot(bang.data)
Lines(Trial.stc{'w'}(:)-60,[],'b');



figure,
eds = {linspace(-3,1,75),linspace(-1,2,75)};
subplot(121);
ind = Trial.stc{'a'};
hist2(log10([bang(ind),vxy(ind,1)]),eds{1},eds{2})
subplot(122);
ind = Trial.stc{'w'};
hist2(log10([bang(ind),vxy(ind,1)]),eds{1},eds{2})



cang = MTADang('data',diff(ang(:,1,1,3).*ang(:,8,10,3)),'sampleRate',ang.sampleRate);
cang.filter('ButFilter',3,5,'low');
cang.data = sq(sum(sqrt(cang.segs(1:size(cang,1),90,0).^2)))';
cang.resample(30);
ind = Trial.stc{'a'};


%% DRV of pitch
d = 2;
cang = ang.copy;
cang.data = [ang(:,1,2,d),...
             ang(:,1,3,d),...
             ang(:,1,5,d),...
             ang(:,2,3,d),...
             ang(:,2,4,d),...
             ang(:,3,4,d),...
             ang(:,3,5,d),...
             ang(:,4,5,d),...
             ang(:,4,7,d),...
             ang(:,5,7,d)];
cang.filter('ButFilter',3,5,'low');
figure,plot(diff(cang.data).*cang.sampleRate);
figure,imagesc(diff(cang.data)'.*cang.sampleRate);caxis([-10,10])


%% DRV of Yaw
d = 1;
cang = ang.copy;
cang.data = [ang(:,1,2,d),...
             ang(:,1,3,d),...
             ang(:,1,5,d),...
             ang(:,2,3,d),...
             ang(:,2,4,d),...
             ang(:,3,4,d),...
             ang(:,3,5,d),...
             ang(:,4,5,d),...
             ang(:,4,7,d),...
             ang(:,5,7,d)];
cang.filter('ButFilter',3,5,'low');

cang.data = circ_dist(cang.data,circshift(cang.data,1)).*cang.sampleRate;
cang.data(cang.data<1e-9) = 1e-9;
cang.data = log10(cang.data);
figure,plot(cang.data);
figure,imagesc(circ_dist(cang.data,circshift(cang.data,1))'.*cang.sampleRate);caxis([-10,10])
