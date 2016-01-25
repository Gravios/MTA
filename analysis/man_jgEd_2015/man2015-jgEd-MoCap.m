
%% Behavioral Segmentation Paper 
%
% Authors: Justin Graboski, Eduardo Hernandez, Anton Sirota
%


%% FIG1:  Recording Set up and processing pipeline  |
% Eduardo has this                                  |
%___________________________________________________|



%% FIG2: Features and Segmentation of behaviors------------------|
%  A: Trajectories of behaving rat                               |
%  B: Feature matrix demonstrating the features we used          |
%  C: Labels corresponding to the hand labeled data              |
%  D:Labels corresponding to the neural network labeled data     |
%  E: t-sne dimensionality reduction method                      |
%  F: Tabel of labeling stats between and animals and between    |
%     labelers                                                   |
%________________________________________________________________|



%% FIG3: Examples of feature dynamics and head body independence ---|
% A  SPEED head and body                                            |
% B  DIRECTION head and body                                        |
% C  PITCH SLPR and SMSU                                            |
% D  Intermarker Distance of head and upper spine                   |
% E  Hand Labels                                                    |
% F  RHM Rhythmic Head Motion                                       |
% G  Mutual information between marker speeds                       |
% H  Time lag of maximum mutual information (ms)                    |
% __________________________________________________________________|



%% FIG4 - Heirarichal Segmentation 

Trial = MTATrial('jg05-20120317');
Stc = Trial.load('stc','hand_labeled_rev1'); 
figPath = '/storage/gravio/manuscripts/man2015-jgEd-MoCap/Figures/Figure_4';

xyz = Trial.load('xyz');

figure
plotSkeleton(xyz,Trial.stc{'r'}(1,1)+40,'t_per',[-60,0])
plotSkeleton(xyz,Trial.stc{'w'}(1,1)+40,'t_per',[-60,0])

clf,plotSkeleton(xyz,Trial.stc{'m'}(2,1)+120,'t_per',[-160,0])

sts = Trial.stc{'n'}.data;
i = 9;
while i~=-1,
    clf,plotSkeleton(xyz,sts(i,1)+40,'t_per',[-80,0]);
    i = figure_controls(gcf,i,'flags','-v');
end


%Segmentation JPDF rear
ang = create(MTADang,Trial,xyz);
vxy = Trial.xyz.copy;
vxy.data = [ang(:,1,2,2),xyz(:,7,3)];
vxy.data = ButFilter(vxy.data,3,4/(xyz.sampleRate/2),'low');
tag = 'Head_HeightVSlower_spine_pitch';
v1 = vxy.copy;
v1.data = vxy(:,1);
v2 = vxy.copy;
v2.data = log10(vxy(:,2));
bhv_JPDF(Trial,v1,v2,70,70,[.25,1.6],[1.4,2.5],...
         'Pitch (SLPR) radians',...
         'Height (Head Front) log10(mm)',{'a-r','r'},...
         tag)


sts = 'rw';
stc = 'rw';
hedgs    = {linspace(.25,1.6,75)};
hedgs(2) = {linspace(1.4,2.5,75)};
edgs    = {linspace(.25,1.6,75)};
edgs(2) = {linspace(1.4,2.5,75)};
[edgs{:}] = get_histBinCenters(edgs);
[X,Y] = meshgrid(edgs{:});

hold on,
for i = 1:numel(sts),
    b = [vxy(Trial.stc{sts(i)},1),log10(vxy(Trial.stc{sts(i)},2))];
    o = hist2(b,hedgs{1},hedgs{2});
    F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
    o = conv2(o,F,'same');
    contour(X,Y,o',[20,20],'linewidth',2.5,'Color',stc(i))
end


%Segmentation JPDF rear alt1
ang = create(MTADang,Trial,xyz);
vxy = Trial.xyz.copy;
vxy.data = [ang(:,3,4,2),xyz(:,7,3)];
vxy.data = ButFilter(vxy.data,3,4/(xyz.sampleRate/2),'low');
tag = 'Head_HeightVSlower_spine_pitch';
v1 = vxy.copy;
v1.data = vxy(:,1);
v2 = vxy.copy;
v2.data = log10(vxy(:,2));
bhv_JPDF(Trial,v1,v2,70,70,[-1,1.7],[1.4,2.6],...
         'Pitch (SMSU) radians',...
         'Height (Head Front) log10(mm)',{'a-r','r'},...
         tag)


sts = 'rw';
stc = 'rw';
hedgs    = {linspace(-1,1.7,75)};
hedgs(2) = {linspace(1.4,2.6,75)};
edgs    = {linspace(-1,1.7,75)};
edgs(2) = {linspace(1.4,2.6,75)};
[edgs{:}] = get_histBinCenters(edgs);
[X,Y] = meshgrid(edgs{:});

hold on,
for i = 1:numel(sts),
    b = [vxy(Trial.stc{sts(i)},1),log10(vxy(Trial.stc{sts(i)},2))];
    o = hist2(b,hedgs{1},hedgs{2});
    F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
    o = conv2(o,F,'same');
    contour(X,Y,o',[20,20],'linewidth',2.5,'Color',stc(i))
end




%Segmentation JPDF walk
vxy = xyz.vel([1,7],[1,2]);
vxy.data = ButFilter(vxy.data,3,4/(xyz.sampleRate/2),'low');
tag = 'speed_SpineL_vs_HeadF';
v1 = vxy.copy; v1.data = log10(vxy(:,1));
v2 = vxy.copy; v2.data = log10(vxy(:,2));
bhv_JPDF(Trial,v1,v2,70,70,[-.7,2],[-.7,2],...
         'Speed (Lower Spine) log10(mm)',...
         'Speed (Head Front) log10(mm)',{'a-r','w','a-w-r'},...
         tag)


sts = 'w';
stc = 'w';
hedgs    = {linspace(-.75,2,75)};
hedgs(2) = {linspace(-.75,2,75)};
edgs    = {linspace(-.75,2,75)};
edgs(2) = {linspace(-.75,2,75)};
[edgs{:}] = get_histBinCenters(edgs);
[X,Y] = meshgrid(edgs{:});

hold on,
for i = 1:numel(sts),
    b = log10(vxy(Trial.stc{sts(i)},:));
    o = hist2(b,hedgs{1},hedgs{2});
    F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
    o = conv2(o,F,'same');
    contour(X,Y,o',[20,20],'linewidth',2.5,'Color',stc(i))
end


% WALK High/Low
ang = create(MTADang,Trial,xyz);
vxy = Trial.xyz.copy;
vxy.data = [ang(:,5,7,2),log10(abs(xyz(:,1,3)))];
%vxy.data = ButFilter(vxy.data,3,4/(xyz.sampleRate/2),'low');
tag = 'Head_pitchVSlower_spine_pitch';
v1 = vxy.copy;
v1.data = vxy(:,1);
v2 = vxy.copy;
v2.data = vxy(:,2);
bhv_JPDF(Trial,v1,v2,70,70,[-1.75,1.75],[.5,2],...
         'Pitch (Head) radians',...
         'Height (Lower Spine) log10(mm)',{'a-r','w'},...
         tag)

sts = 'w';
stc = 'r';
hedgs    = {linspace(.25,1.6,75)};
hedgs(2) = {linspace(-1.8,1.8,75)};
edgs    = {linspace(.25,1.6,75)};
edgs(2) = {linspace(-1.8,1.8,75)};
[edgs{:}] = get_histBinCenters(edgs);
[X,Y] = meshgrid(edgs{:});

hold on,
for i = 1:numel(sts),
p    b = [vxy(Trial.stc{sts(i)},1),vxy(Trial.stc{sts(i)},2)];
    o = hist2(b,hedgs{1},hedgs{2});
    F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
    o = conv2(o,F,'same');
    contour(X,Y,o',[20,20],'linewidth',2.5,'Color',stc(i))
end

% EXP groom seg
nind = nniz(ang(:,1,4,3));
figure,hist2([ang(nind,1,4,3),ang(nind,2,3,3)],linspace(110,170,100),linspace(30,70,100)),caxis([0,100]);


%% Figure 5 Fine movement characterization
figPath = '/gpfs01/sirota/homes/gravio/Documents/Manuscripts/Vicon_Methods_2015/Figures/Figure_4';


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fig:5:A - Fast scans of head
%                            
% Description: {}



figure,plot(circ_dist(ang(:,5,7,1),ang(:,5,10,1)))
Lines(Trial.stc{'w'}(:),[],'k');
Lines(Trial.stc{'n'}(:),[],'g');
hold on,plot(nunity(diff(ButFilter(ang(:,5,10,3),3,[1,30]/(ang.sampleRate/2),'bandpass'))),'m')

dang = Trial.ang.copy;
dang.data = [circ_dist(ang(:,5,7,1),ang(:,5,10,1)),[0;diff(ButFilter(ang(:,5,10,3),3,[1,30]/(ang.sampleRate/2),'bandpass'))]];
sparm = struct('nFFT',2^9,...
               'Fs',dang.sampleRate,...
               'WinLength',2^7,...
               'nOverlap',2^7*.875,...
               'FreqRange',[1,50]);
[ys,fs,ts,ps] = fet_spec(Trial,dang,'mtchglong',true,'defspec',sparm);

figure,imagesc(ts,fs,(log10(ys(:,:,1,1)))');axis xy,
Lines(Trial.stc{'n',ys.sampleRate}(:),[],'g')
Lines(Trial.stc{'w',ys.sampleRate}(:),[],'k')

sang = ang.copy;
sang.resample(ys);

ind = Trial.stc{'a'};
figure,hist2([sang(ind,5,7,2),log10(ys(ind,40))],linspace(-1,1.6,100),linspace(-10,-3,100))




phs.data = circ_dist(ang(:,5,7,1),ang(:,5,10,1));
p = phs.phase([1,4]);
phs.data = [0;diff(ButFilter(ang(:,5,10,3),3,[1,30]/(ang.sampleRate/2),'bandpass'))];
s = phs.phase([6,13]);

ind = Trial.stc{'w'};
figure,hist2([p(ind,1),s(ind,:)],linspace(-pi,pi,100),linspace(-pi,pi,100))

mind = LocalMinima(-ButFilter(circ_dist(ang(:,5,7,1),ang(:,5,10,1)),3,[1,4]/(ang.sampleRate/2),'bandpass'),10,-.1);
aind = Trial.stc{'a'}.cast('TimeSeries');
pind = false([aind.size]);
pind(mind) = true;
nind = pind&logical(aind.data);

[Co,f] = Comodugram(dang(Trial.stc{'a'},:),2^9,dang.sampleRate,[1,20],2^8,[],'linear');


%older stuff
hfig = figure(401);

plot([circ_dist(ang(:,5,7,1),ang(:,2,3,1)),...
      circ_dist(ang(:,1,2,1),ang(:,3,4,1))]);

figure,plot(diff([circ_dist(ang(:,5,7,1),ang(:,2,3,1))]))

dang = Trial.ang.copy;
%dang.data = diff([0;[circ_dist(ang(:,5,7,1),ang(:,2,3,1))-...
%             circ_dist(ang(:,3,4,1),ang(:,1,2,1))]]);
dang.data = diff(circ_dist(ang(:,5,7,1),ang(:,2,4,1)));

sparm = struct('nFFT',2^9,...
               'Fs',dang.sampleRate,...
               'WinLength',2^7,...
               'nOverlap',2^7*.875,...
               'FreqRange',[1,50]);


[ys,fs,ts,ps] = fet_spec(Trial,dang,'mtchglong',true,'defspec',sparm);

figure,imagesc(ts,fs,(log10(ys(:,:,1,1)))');axis xy,
Lines(Trial.stc{'n',ys.sampleRate}(:),[],'g')


figure,imagesc(ts,fs,ys(:,:,1,2)');axis xy,
Lines(Trial.stc{'w',ys.sampleRate}(:),[],'k',[],3);

figure,sp = [];
sp(1) = subplot(211),imagesc(ts,fs,log10(ys(:,:,1,1))');axis xy,
sp(2) = subplot(212),imagesc(ts,fs,log10(ys(:,:,2,2))');axis xy,
linkaxes(sp,'xy');
Lines(Trial.stc{'w',ys.sampleRate}(:),[],'k',[],3);


%% FIG6 RHM NCP
figPath = '/storage/gravio/manuscripts/man2015-jgEd-MoCap/Figures/Figure_5';


% SUBPLOT 6.A
% NAME Time series RHM NCP scaled signals
% DESCRIPTION 
Trial = MTATrial('Ed05-20140529','all','ont');
rhm = fet_rhm(Trial);
ncp = fet_ncp(Trial);

hfig = figure(38447);
plot((1:rhm.size(1))./rhm.sampleRate,bsxfun(@times,nunity([rhm.data,ncp.data]),[3,1]))
xlabel('Time (s)');
ylabel('AU');
title({'Time Series of',['Rythmic Head Motion (RHM) and Nasal Cavity ' ...
                    'Pressure (NCP)']})
legend('RHM','NCP');

saveas(hfig,fullfile(figPath,...
                     ['Fig3-Features-sample_' num2str(s) '_' Trial.filebase '.png']),...
       'png');
saveas(hfig,fullfile(figPath,...
                     ['Fig3-Features-sample_' num2str(s) '_' Trial.filebase '.eps']),...
       'epsc');



% PSD Time Series
Trial = MTATrial('Ed05-20140529','all','ont');
[ncp,nfs,nts] = fet_ncp(Trial,[],'mtchglong',2);
[rhm,rfs,rts] = fet_rhm(Trial,[],'mtchglong');

sp = [];
hfig = figure(87372);
set(hfig,'position',[100,100,640,330])
set(hfig,'paperposition',[0,0,640,330])
% RHM Spectrogram
sp(1) = subplot(211);
imagesc(rts,rfs,log10(rhm.data)');axis xy, 
title('Rythmic Head Motion PSD')
ylabel('Frequency (Hz)')
%xlim([35,60])
caxis([-6,-3.5])
colorbar

% NCP Spectrogram
sp(2) = subplot(212);
imagesc(nts,nfs,log10(ncp.data)');axis xy, 
title('Nasal Cavity Pressure PSD')
ylabel('Frequency (Hz)')
xlabel('Time (s)')
%xlim([35,60])
caxis([2,6])
colorbar
linkaxes(sp,'xy')

%saveas(hfig,fullfile(figPath,['Fig5-RHM_NCP_PSD-' Trial.filebase '.png']),'png');
%saveas(hfig,fullfile(figPath,['Fig5-RHM_NCP_PSD-' Trial.filebase '.eps']),'epsc');


% SUBPLOT 6.F
% NAME RHM NCP phase diff
% DESCRIPTION > The trough of each breathing cycle was
% detected and created a JPDF The circular difference 
% between the RHM and NCP phases and the trough pressure 
Trial = MTATrial('Ed05-20140529','all','ont');
[ncp] = fet_ncp(Trial,'chans',2);
[rhm] = fet_rhm(Trial);
nphs = ncp.phase([6,14]);
rphs = rhm.phase([6,14]);
mind = LocalMinima(ncp.data,8,0);

mind = LocalMinima(ncp.data,8,0);
out = PPC([nphs(mind),rphs(mind)]);
figure,
hist2([circ_dist(nphs(mind),rphs(mind)),ncp(mind)],30,linspace(-6000,-1000,30))
xlabel('Phase Difference (radians)')
ylabel('Peak Negative Presure (A.U.)')
title(['NCP RHM Phase Difference at Inhalation [PPC: ' ...
       num2str(round(out,4)) ']'])


%% Group stats for RHM NCP
%
% Ed01
% Ed03
% Ed05
% Ed10
%


% RESP Breathing
ti = 1:ang.size(1);
ts = ti/ang.sampleRate;
xl = round([243267,244603]/ang.sampleRate);

figure,hold on,
plot(ts,nunity(ButFilter(ang(ti,2,4,3),3,[.5,30]/(ang.sampleRate/2),'bandpass')).*10)
plot(ts,nunity(ButFilter(ncp(ti)      ,3,[.5,30]/(ang.sampleRate/2),'bandpass')),'g')
xlim(xl)
title({'Respiration During Immobility','r(SLSU) and NCP bandpass [0.5,30]'})
ylabel('A.U.')
xlabel('Time (s)')
legend('r(SLSU)','NCP')


