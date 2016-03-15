
Session = MTASession('OP01-20160201','sof',true,'','vicon','nn');
QuickTrialSetup(Session);
Trial = MTATrial('OP01-20160201','all','sof');

xyz = Trial.load('xyz');
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,20,'low');
ang = create(MTADang,Trial,xyz);
fang = create(MTADang,Trial,fxyz);

%name='dpitchh';label='dph';key='p';
%fet = MTADfet.encapsulate(Trial,circ_dist(fang(:,2,6,1),circshift(fang(:,2,6,1),-10)),fang.sampleRate,name,label,key);

vxy = xyz.copy;
vxy.data = ButFilter(vxy.data,3,2.4/(xyz.sampleRate/2),'low');
vxy = vxy.vel([1,6],[1,2]);
% 
% figure,
% t = [1:fang.size(1)]/xyz.sampleRate;
% plot(t(1:end-1),diff(ButFilter(circ_dist(fang(:,1,3,1),fang(:,1,2,1)),3,4/(0.5*ang.sampleRate),'low'))),
% hold on,plot([t],nunity(vxy(:,:),[],[],[],[],1)/50)
% xlim([  20082.6729409171 , 24422.8100556091 ]/xyz.sampleRate)
% ylim([-.03,.08])

%tag = 'Head_HeightVSlower_spine_pitch';
v1 = vxy.copy;
v1.data = log10(vxy(:,1));
edx = linspace(-3,2,70);
edy = linspace(-3,2,70);
v2 = vxy.copy;
v2.data = log10(vxy(:,2));

ind = nniz(vxy);
figure,hist2([v1(ind),v2(ind)],edx,edy)


%% RAT 
Trial = MTATrial('jg05-20120310');
xyz = Trial.load('xyz');
%fxyz.filter('ButFilter',3,20,'low');
%ang = create(MTADang,Trial,xyz);
%fang = create(MTADang,Trial,fxyz);

%name='dpitchh';label='dph';key='p';
%fet = MTADfet.encapsulate(Trial,circ_dist(fang(:,2,6,1),circshift(fang(:,2,6,1),-10)),fang.sampleRate,name,label,key);

vxy = xyz.copy;
vxy.data = ButFilter(vxy.data,3,2.4/(xyz.sampleRate/2),'low');
vxy = vxy.vel([1,6],[1,2]);

%tag = 'Head_HeightVSlower_spine_pitch';
v1 = vxy.copy;
v1.data = log10(vxy(:,1));
edx = linspace(-3,2,100);
edy = linspace(-3,2,100);
v2 = vxy.copy;
v2.data = log10(vxy(:,2));

ind = nniz(vxy);
figure,hist2([v1(ind),v2(ind)],edx,edy)



