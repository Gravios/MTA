
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





Session = MTASession('nk00-20160315','sof',true,'','vicon','nn');
QuickTrialSetup(Session);
Trial = MTATrial('nk00-20160315','all','sof');


xyz = Trial.load('xyz');
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,20,'low');
%ang = create(MTADang,Trial,xyz);
fang = create(MTADang,Trial,fxyz);
t = [1:fang.size(1)]/xyz.sampleRate;
vxy = xyz.copy;
vxy.data = ButFilter(vxy.data,3,2.4/(xyz.sampleRate/2),'low');
vxy = vxy.vel([1,9],[1,2]);
%name='dpitchh';label='dph';key='p';
%fet = MTADfet.encapsulate(Trial,circ_dist(fang(:,2,6,1),circshift(fang(:,2,6,1),-10)),fang.sampleRate,name,label,key);
name='spine_wag';label='spw';key='w';
fet = MTADfet.encapsulate(Trial,circ_dist(fang(:,1,3,1),circshift(fang(:,1,4,1),-10)),fang.sampleRate,name,label,key);



 [nvxy,vMean,vStd] = nunity(vxy(:,:),[],[],[],[],1);
 figure,
 subplot(211),plot(t(1:end-1),diff(ButFilter(circ_dist(fang(:,1,2,1),fang(:,1,4,1)),3,4/(0.5*fang.sampleRate),'low'))),
 hold on,plot([t],nvxy/50)
 xlim([40,70])
 ylim([-0.03,0.05])
 title('Scaled Features Derived from the XY Plane: Mouse')
 xlabel('Time (s)')
 legend('Spine Curvature','Body Speed', 'Head Speed');
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


def_spec = struct('nFFT',2^9,'Fs',fet.sampleRate,  'WinLength',2^8,'nOverlap',2^8*.875,'FreqRange',[1,20]);


% scanning of head
name='head scanning';label='scn';key='n';
fet = MTADfet.encapsulate(Trial,circ_dist(ang(:,1,4,1),ang(:,7,9,1)),ang.sampleRate,name,label,key);
fet.resample(120);
[ys,fs,ts] = fet_spec(Trial,fet,'mtchglong',true,'defspec',def_spec,'overwrite',true);
figure,
subplot(2,1,1), imagesc(ts,fs,log10(ys.data')),axis xy, colormap jet,caxis([-6,-4])


% Sniff/whisk pitch of head
name='sniff_pitch';label='sniffp';key='p';
fet = MTADfet.encapsulate(Trial,fang(:,7,9,2),fang.sampleRate,name,label,key);
[ys,fs,ts] = fet_spec(Trial,fet,'mtchglong',true,'defspec',def_spec,'overwrite',true);
figure, subplot(2,1,1),cla, imagesc(ts,fs,log10(ys.data')),axis xy, colormap jet,caxis([-6,-4])
xlim([40,70])

[rhm,fs,ts] = fet_rhm(Trial,[],'mtchglong',true,'defspec',def_spec,'overwrite',true);
figure,imagesc(ts,fs,log10(rhm.data')),axis xy, colormap jet,

%% RAT 
Trial = MTATrial('Ed03-20140625');
xyz = Trial.load('xyz');
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,20,'low');
%ang = create(MTADang,Trial,xyz);
fang = create(MTADang,Trial,fxyz);
t = [1:xyz.size(1)]/xyz.sampleRate;
vxy = xyz.copy;
vxy.data = ButFilter(vxy.data,3,2.4/(xyz.sampleRate/2),'low');
vxy = vxy.vel([1,7],[1,2]);

%name='dpitchh';label='dph';key='p';
%fet = MTADfet.encapsulate(Trial,circ_dist(fang(:,1,6,1),circshift(fang(:,2,6,1),-10)),fang.sampleRate,name,label,key);

 
 subplot(212),cla,plot(t(1:end-1),diff(ButFilter(circ_dist(fang(:,1,2,1),fang(:,1,4,1)),3,4/(0.5*fang.sampleRate),'low'))),
 hold on,plot([t],nunity(vxy(:,:),[],vMean,vStd,[],1)/50)
 xlim([0,30])
 ylim([-0.03,0.05])
 title('Scaled Features Derived from the XY Plane: Mouse')
 xlabel('Time (s)')
 legend('Spine Curvature','Body Speed', 'Head Speed');
 

%tag = 'Head_HeightVSlower_spine_pitch';
v1 = vxy.copy;
v1.data = log10(vxy(:,1));
edx = linspace(-3,2,100);
edy = linspace(-3,2,100);
v2 = vxy.copy;
v2.data = log10(vxy(:,2));

ind = nniz(vxy);
figure,hist2([v1(ind),v2(ind)],edx,edy)

name='head scanning';label='scn';key='n';
fet = MTADfet.encapsulate(Trial,circ_dist(ang(:,1,4,1),ang(:,5,7,1)),ang.sampleRate,name,label,key);
[ys,fs,ts] = fet_spec(Trial,fet,'mtchglong',true,'defspec',def_spec,'overwrite',false);
subplot(2,1,2), imagesc(ts,fs,log10(ys.data')),axis xy, colormap jet,caxis([-6,-4])


name='sniff_pitch';label='sniffp';key='p';
fet = MTADfet.encapsulate(Trial,fang(:,7,9,2),fang.sampleRate,name,label,key);
[ys,fs,ts] = fet_spec(Trial,fet,'mtchglong',true,'defspec',def_spec,'overwrite',false);,
subplot(2,1,2),cla, imagesc(ts,fs,log10(ys.data')),axis xy, colormap jet,caxis([-6,-4])
xlim([0,30])

