
%% figure 1
Trial = MTATrial('Ed10-20140812');
Trial = MTASession('jg05-20120317');
Trial = MTASession('er06-20130612');
Trial.stc.updateMode('auto_wbhr');
Trial.stc.load;

xyz = Trial.load('xyz');
xyz.filter(gtwin(.3,xyz.sampleRate));
ang = Trial.ang.copy;
ang.create(Trial,xyz);


%% Fig:1:E
% marker error
hfig = figure(838884);
edges = 43:.05:45;
%rTrial.stc{'w'}
rper = Trial.stc{'r'}.cast('TimeSeries');
wper = Trial.stc{'w'}.cast('TimeSeries');
vel = xyz.vel(7,[1,2]);
%ind = ~(rper.data|wper.data);
ind = vel<2;
N = histc(ang(ind,5,7,3),edges);
bar(edges,N,'histc')
title('Distance between the back and front head markers')
xlabel('Inter Marker Distance (mm)');
ylabel('Count');

vel = xyz.vel([1:9]);
vel.data(vel<0.01) = 0.01;
vel.data = log10(vel.data);

m = {'head_back','head_right'};mb = [];
vbins = -.5:.05:2;
[~,mb(:,1)] = histc(vel(:,m{1}),vbins);
[~,mb(:,2)] = histc(vel(:,m{2}),vbins);

mb = MTADxyz('data',mb,'sampleRate',xyz.sampleRate);

ind = resample(Trial.stc{'r'}.cast('TimeSeries'),xyz)&nniz(mb);



A = accumarray(mb(ind,:),ang(ind,m{1},m{2},3),repmat(numel(vbins),[1,2]),@mean);
B = accumarray(mb(ind,:),ang(ind,m{1},m{2},3),repmat(numel(vbins),[1,2]),@std);
S = accumarray(mb(ind,:),ind(ind),repmat(numel(vbins),[1,size(mb,2)]),@sum);
A(S<100)=nan;
B(S<100)=nan;

% mean marker distance vs speed
figure,
subplot(121)
imagescnan({vbins,vbins,A'},prctile(A(nniz(A(:))),[5,95]),[],true,[0,0,0]),axis xy,
subplot(122)
imagescnan({vbins,vbins,B'},[0.4,5],[],true,[0,0,0]),axis xy

figure,
hist(ang(Trial.stc{'a'},'head_back','head_front',3),20:.02:37)
figure,
hist(ang(Trial.stc{'a'},'head_right','head_front',3),20:.02:37)
figure,
hist(ang(Trial.stc{'a'},'head_right','head_left',3),20:.02:37)
figure,
hist(ang(Trial.stc{'a'},'head_back','head_right',3),10:.02:37)


%% Fig:1:F Breathing (maybe)

Trial = MTASession('Ed10-20140815');

xyz = Trial.load('xyz').filter(gtwin(.25,Trial.xyz.sampleRate));;
ang = Trial.ang.copy;
ang.create(Trial,xyz);


ind = 866330:867110;
mar = {'pelvis_root','spine_upper'};

hfig = figure(8482838);
sp=[];
sp(1) = subplot(211);
plot(round((ind-ind(1))/xyz.sampleRate,2),ang(ind,mar{1},mar{2},3));
%ylim([131.9,132.6])
title('Distance Between the Pelvis and the Upper Spine')
ylabel('Distance (mm)')
sp(2) = subplot(212);
plot(round((ind-ind(1))/xyz.sampleRate,2),lfp.data(ind));
title('Nasal Cavity Pressure Sensor')
ylabel('NCP (mV)')
ylim([-4000,2800])
xlabel('Time (s)')
linkaxes(sp,'x')
xlim([.25,6]);

saveas(hfig,'/gpfs01/sirota/homes/gravio/Documents/Manuscripts/Vicon_Methods_2015/Fig1F-alt1.png','png');
saveas(hfig,'/gpfs01/sirota/homes/gravio/Documents/Manuscripts/Vicon_Methods_2015/Fig1F-alt1.eps','eps');


% $$$ xyz = Trial.load('xyz');
% $$$ % lfp = Trial.lfp.copy;
% $$$ % lfp.create(Trial,66);
% $$$ % lfp.resample(xyz);
% $$$ % 
% $$$ % 
% $$$ % 
% $$$ % specParms = struct('nFFT',2^11,...
% $$$ %                    'Fs',lfp.sampleRate,...
% $$$ %                    'WinLength',2^10,...
% $$$ %                    'nOverlap',2^10*.875,...
% $$$ %                    'FreqRange',[1,15]);
% $$$ % 
% $$$ % 
% $$$ % lfp.data = cat(2,lfp.data,xyz(:,1,3));
% $$$ % lfp.data = cat(2,lfp.data,xyz(:,3,3));
% $$$ 
% $$$ 
% $$$ xyz.filter(gtwin(.5,xyz.sampleRate));
% $$$ ang = Trial.ang.copy;
% $$$ ang.create(Trial,xyz);
% $$$ 
% $$$ figure,
% $$$ plot(sq(ButFilter(xyz(:,[1:3],3)));
% $$$ 
% $$$ figure,plot(ang(:,1,3,3));
% $$$ 
% $$$ 
% $$$ rhm = fet_rhm(Trial);
% $$$ lfp.data = cat(2,lfp.data,rhm);
% $$$ 
% $$$ [ys,fs,ts,phi,fstat] = fet_spec(Trial,lfp,'mtchglong','overwrite',true);
% $$$ 
% $$$ figure,imagesc(ts,fs,ys(:,:,1,2)'),axis xy



%% Figure 2 Trajectories and behavioral labeling
Trial = MTATrial('jg05-20120317');
labelMode = 'auto_wbhr';
labelMode = 'manual_mknsrw';
labelMode = 'hand_labeled';

xyz = Trial.load('xyz');
Trial.stc.load(Trial,labelMode); 


%% Fig:2:A

% repair origins of epochs ( use to be in sampling rate of the
% object, now in the absolute timeline in seconds
% $$$ Trial.stc.updateMode(labelMode);
% $$$ Trial.stc.updatePath(Trial.spath);
% $$$ Trial.stc.load;
% $$$ Trial.stc.updatePath(Trial.spath)
% $$$ Trial.stc.updateMode('hand_labeled');
% $$$ Trial.stc.save(1);
% $$$ for i = 1:numel(Trial.stc.states),
% $$$     %Trial.stc.states{i}.origin = Trial.stc.states{i}.origin/Trial.stc.states{i}.sampleRate;
% $$$     Trial.resync(Trial.stc.states{i});
% $$$ end
% $$$ Trial.stc.save(1);


% jg05-20120310
% 31600 - 33700 turn -> walk -> rear

% jg05-20120317
% 26800 - 27800 turn -> walk -> rear
%# create coordinates


ftit = 'Turning';
ind = 26764;
perind = (ind-140):(ind);
figure,hold on
for i= [1:4,5,7];
    plot3(xyz(perind,i,1),xyz(perind,i,2),xyz(perind,i,3))
end
plotSkeleton(xyz,perind(end),'surface');
zlim([0,300]);
title(ftit)
set(gca,'YTickLabelMode','manual');
set(gca,'YTickLabel',{});
set(gca,'XTickLabelMode','manual');
set(gca,'XTickLabel',{});
set(gca,'ZTickLabelMode','manual');
set(gca,'ZTickLabel',{});

% features

figure,imagesc((1:size(fet,1))/xyz.sampleRate,1:size(fet,2),fet'),caxis([-2,3])
Lines(Trial.stc{'r'}(:)/xyz.sampleRate,[],'r',[],3);
Lines(Trial.stc{'w'}(:)/xyz.sampleRate,[],'k',[],3);
colormap jet
xlim([31600,33700]/xyz.sampleRate)
title('Behavioral Segmentation Features');
set(gca,'YTickLabelMode','manual');
set(gca,'YTickLabel',{});
ylabel('Features');
xlabel('Time (s)')





%% Figure 3 JPDFs


Trial = MTATrial('jg05-20120310');

% Vars of interest 
xyz = Trial.load('xyz');
ang = Trial.load('ang');
vxy = vel(Trial.load('xyz').filter(gtwin(.5,xyz.sampleRate)),{'spine_lower','head_front'},[1,2])
vxy.data(vxy.data<.001) = 0.001;
rol = fet_roll(Trial,[],'default');


%% Fig:3:A - Rearing from everything else

tag = 'Rearing-Everything';
v1 = Trial.ang.copy;
v1.data = ang(:,2,3,3);
v2 = Trial.ang.copy;
v2.data = ang(:,3,4,2);
bhv_JPDF(Trial,v1,v2,70,70,...
         'Distance(Pelvis,Spine Middle) (mm)',...
         'Upper Spine Pitch (rad)',{'a-r','r'},tag)


tag = 'pelvis2spineM_vs_headroll';
v1 = Trial.ang.copy;
v1.data = ang(:,2,3,3);pp
v2 = rol.copy;
bhv_JPDF(Trial,v1,v2,70,70,...
         'Distance(Pelvis,Spine Middle) (mm)',...
         'Head roll (rad)','rwhl',tag)


%% Fig:3:B - walk from everything else
tag = 'speed_SpineL_vs_HeadF';
v1 = vxy.copy;
v1.data = log10(vxy(:,1));
v2 = vxy.copy;
v2.data = log10(vxy(:,2));
bhv_JPDF(Trial,v1,v2,70,70,...
         'Speed (Lower Spine) log10(mm)',...
         'Speed (Head Front) log10(mm)',{'a-r','w','a-w'},tag)



%% Fig:3:C - walk to high and low walk

tag = 'Walk_High-low';
v1 = Trial.ang.copy;
v1.data = ang(:,5,7,3);
v2 = Trial.ang.copy;
v1.data = ang(:,4,5,3);
bhv_JPDF(Trial,v1,v2,70,70,...
         'Distance(SpineU,HeadB) (mm)',...
         'Head roll (rad)',{'w','h','l'},tag)


tag = 'HeightSL_vs_headPitch';
v1 = Trial.xyz.copy;
v1.data = xyz(:,1,3);
v2 = Trial.ang.copy;
v2.data = ang(:,5,7,2);
bhv_JPDF(Trial,v1,v2,70,70,...
         'Hight(SpineL) (mm)',...
         'Head Pitch (rad)','whl',tag)



%% Figure 4 Fine movement characterization
Trial = MTATrial('jg05-20120310');
xyz = Trial.load('xyz');
xyz.filter(gtwin(.2,xyz.sampleRate));
ang = Trial.ang.copy;
ang.create(Trial,xyz);

figure,
mp = [1,2;2,3;3,4;4,5;5,7]';
vang = [];
for i = mp,
vang(:,end+1) = circ_dist(ang(:,i(1),i(2),1),circshift(ang(:,i(1),i(2),1),-20));
end
plot(diff(vang))
%xlim([31600,33700]);
Lines(Trial.stc{'r'}(:),[],'r',[],3);
Lines(Trial.stc{'w'}(:),[],'k',[],3);



%% Figure 5 RHM (rythmic head motion) feature versus NCP (nasal cavity pressure)
Trial = MTATrial('Ed10-20140815');

%generate features
[rhm,fs,ts] = fet_rhm(Trial,[],'wcsd');
ncp = fet_ncp(Trial,[],'wcsd',66);
%plot features with linked axes

figure,
sp(1) = subplot(211);
imagesc(ts,fs,log10(rhm.data)'),axis xy,caxis([-5,-3.1])
title('Rhythmic Head Motion (RHM)')
ylabel('frequency (Hz)')
xlabel('Time (s)')
sp(2) = subplot(212);
imagesc(ts,fs,log10(ncp.data)'),axis xy,caxis([3,4.2])
title('Nasal Cavity Pressure (NCP)')
ylabel('frequency (Hz)')
xlabel('Time (s)')
linkaxes(sp,'xy');
xlim([700,860])


Trial = MTATrial('Ed10-20140812');

lfp = Trial.lfp.copy;
lfp.load(Trial,[8,35]);

specParms = struct('nFFT',2^9,...
                    'Fs',lfp.sampleRate,...
                    'WinLength',2^8,...
                    'nOverlap',2^8*.875,...
                    'FreqRange',[20,150]);

[ys,fs,ts] = fet_spec(Trial,lfp,'mtchglong',true,lfp.sampleRate,specParms);

