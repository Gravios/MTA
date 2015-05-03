%% new Figure2 
Trial = MTATrial('jg05-20120317');
figPath = '/storage/gravio/manuscripts/man2015-jgEd-MoCap/Figures/Figure_2';

%Features 

% SPD Low Pass Filtered 4Hz
xyz = Trial.load('xyz');
vl = vel(xyz,1:8,[1,2]);
vl.data = ButFilter(vl.data,3,4/(vl.sampleRate/2),'low');

% RHM Spectrum 1-30Hz
dsp.nFFT= 512;
dsp.Fs = 119.881035;
dsp.WinLength = 128;
dsp.nOverlap = 112;
dsp.FreqRange = [1 30];
[rhm,fs,ts] = fet_rhm(Trial,[],'mtchglong','defspec',dsp);

% XYZ Low Pass Filtered 4Hz
xyz.data = ButFilter(xyz.data,3,[4]/(xyz.sampleRate/2),'low');

% ANG From Low Pass Filetered (4Hz) XYZ data
ang = create(Trial.ang.copy,Trial,xyz);


sfet = Trial.ang.copy;
% $$$ sfet.data = sum(abs([circ_dist(ang(:,1,3,1),ang(:,1,2,1)),...
% $$$                      circ_dist(ang(:,2,4,1),ang(:,2,3,1)),...
% $$$                      circ_dist(ang(:,3,5,1),ang(:,3,4,1)),...
% $$$                      circ_dist(ang(:,4,7,1),ang(:,4,5,1))]),2);
sfet.data = -nunity(circ_dist(ang(:,5,7,1),ang(:,1,4,1)));

swg = Trial.ang.copy;
swg.data = nunity(circ_dist(ang(:,1,4,1),ang(:,2,4,1)));


%% mutinfo crap
v = log10(vl.data(:,[1:8]));

edges = linspace(-.5,2,64);
sbound = -130:130;
ixy = zeros([numel(sbound),size(v,2),size(v,2)]);

padding = [0,0];
vind = logical(subsref(cast(resample(Trial.stc{'a'}+padding,xyz),'TimeSeries'),substruct('.',{'data'})));
nind = numel(vind);

s = 1;
for m = 1:size(v,2)
for o = 1:size(v,2)
for shift = sbound
[out,xb,yb,p]=hist2([v(vind,m),circshift(v(vind,o),shift)],edges,edges);
pxy = out./nind;
px = histc(v(vind,m),xb);
px = px(1:end-1)/nind;
py = histc(circshift(v(vind,o),shift),yb);
py = py(1:end-1)/nind;
ixy(s,m,o) = nansum(nansum(pxy.*log2(pxy./(px*py'))));
s = s+1;
end
s = 1;
end
end

[mixy,sixy] = max(ixy);
mixy = sq(mixy);
sixy = sq(sixy)-ceil(numel(sbound)/2);
sixy = sixy([1:4,7],[1:4,7]);



%% Fig Parameters
sts = 'rwnms';
stc = 'rcymg';
edgs    = {linspace(-.5,2,75)};
edgs(2) = {linspace(-.5,2,75)};
[edgs{:}] = get_histBinCenters(edgs);
[X,Y] = meshgrid(edgs{:});
ns = 6;

xpers = bsxfun(@plus,Trial.stc{'w',1}(1:3:end,1),[-15,15]);
xpers(xpers(:,1)<1,:) = [];
xpers((xpers(:,2)-Trial.sync(end))>0,:) = [];
s = 20;




hfig = figure(2);


% JPDF - Head/Body speed
subplot2(ns,4,[3,4],4);cla;
b = log10([median(vl(nniz(vl),[1:2]),2),median(vl(nniz(vl),[5:8]),2)]);
hist2(b,edgs{1},edgs{2});
xlabel('log10 body speed (cm/s)');
ylabel('log10 head speed (cm/s)');
title('JPDF of log10 head and body speeds');

hold on,
for i = 1:numel(sts),
    b = log10([median(vl(Trial.stc{sts(i)},1:2),2),median(vl(Trial.stc{sts(i)},5:8),2)]);
    o = hist2(b,linspace(-.5,2,75),linspace(-.5,2,75));
    F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
    o = conv2(o,F,'same');
    contour(X,Y,o',[20,20],'linewidth',1.5,'Color',stc(i))
end


ndag = double(~eye(8));
subplot2(ns,4,[1,2],4);
imagesc(mixy(:,:).*ndag);
colorbar
title('mutual information between marker speeds');
set(gca,'YtickMode','manual');
set(gca,'Ytick',1:8);
set(gca,'YtickLabelMode','manual');
set(gca,'YtickLabel',vl.model.ml('short'));

subplot2(ns,4,[5,6],4);
imagesc(sixy(:,:)/vl.sampleRate*1000);
colorbar
title('time lag of maximum mutual information (ms)')

set (hfig,'position',[0,0,1000,700])
set (hfig,'paperposition',[0,0,1000/100,700/100])
set (hfig,'PaperType','a3');
ns = 6;
periodInds = 1:size(xpers,1);
periodInds = 32;
periodInds = [14,24,29,35,36,48,54,55,84,89,92];
for s = periodInds
    ind = round(xpers(s,:).*xyz.sampleRate);
    ind = ind(1):ind(2);
    i = 1;
    
    %figure(201) % SPEED head and body
    subplot2(ns,4,i,[1:3]);i=i+1;cla
    plot(ind/vl.sampleRate,[median(vl(ind,1:2),2),median(vl(ind,5:8),2)]),axis tight
    title('xy speed of head and body')
    ylabel('Speed (cm/s)');
    ylim([0,80])
    set(gca,'XTickLabelMode','manual');
    set(gca,'XTickLabel',{});
    set(gca,'TickDir','out');
    
    %figure(202) % DIRECTION head and body
    subplot2(ns,4,i,[1:3]);i=i+1;cla
    plot(ind/vl.sampleRate,[sfet(ind),swg(ind)])
    axis tight
    title('diff ang(head,body) and spine waggle')
    ylabel('Circular difference normalized (AU)');
    set(gca,'XTickLabelMode','manual');
    set(gca,'XTickLabel',{});
    set(gca,'TickDir','out');

    
    %figure(203)%  PITCH SLPR and SMSU
    subplot2(ns,4,i,[1:3]);i=i+1;cla
    plot(ind/vl.sampleRate,[ang(ind,1,3,2),ang(ind,5,7,2)]),axis tight
    title('Pitch of Body and Head')
    ylabel('Pitch (radians)');
    ylim([-pi/2,pi/2])
    set(gca,'XTickLabel',{});
    set(gca,'TickDir','out');    
    set(gca,'TickDir','out');

    
    %figure(204)
    subplot2(ns,4,i,[1:3]);i=i+1;cla
    plot(ind/vl.sampleRate,[ang(ind,4,5,3)]),axis tight%,ang(ind,4,7,3)])
    title('Intermarker Distance of head and upper spine')
    ylabel('Distance (mm)');
    set(gca,'XTickLabelMode','manual');
    set(gca,'XTickLabel',{});
    set(gca,'TickDir','out');
    
    %figure(205)
    subplot2(ns,4,i,[1:3]);i=i+1;cla
    plotSTC(Trial.stc,1,'patch',{'rear','walk','turn','groom'},'rbgm');
    title('hand labeling')
    xlim(xpers(s,:));
    ylim([0,1])
    set(gca,'TickDir','out');
    set(gca,'XTickLabelMode','manual');
    set(gca,'XTickLabel',{});
   
    
    %figure(206) RHM Rhythmic Head Motion
    subplot2(ns,4,i,[1:3]);i=i+1;cla
    sind = round(xpers(s,:).*rhm.sampleRate);
    sind = sind(1):sind(2);
    title('RHM Rhythmic Head Motion')
    imagesc(ts(sind),fs,log10(rhm(sind,:))'),axis xy ,caxis([-6,-3.5])
    ylabel('Frequency (Hz)')
    xlabel('Time (s)')
    set(gca,'TickDir','out');
    %ca = colorbar;
    %set(ca,'positon',[0.13,0.11,0.568882978723404,0.102587412587413])

    
    
    saveas(hfig,fullfile(figPath,['Fig2-Features-sampleN_' num2str(s) '_' Trial.filebase '.png']),'png');
    saveas(hfig,fullfile(figPath,['Fig2-Features-sampleN_' num2str(s) '_' Trial.filebase '.eps']),'epsc');


end



%Alt MUT info

figure,
plot(round(sbound/xyz.sampleRate*1000,3),[ixy(:,7,1),ixy(:,7,2),ixy(:,7,3),ixy(:,7,4)])
hold on,plot(sbound(sixy(5,1)),mixy(7,2),
xlim([round(sbound(1)/xyz.sampleRate*1000,3),round(sbound(end)/xyz.sampleRate*1000,3)]);
legend('HF->SL','HF->PR','HF->SM','HF->SU')
title({'Mutual Information Between Marker Speeds','with Varying Time Lags'})
ylabel('Mutual Information (bits)')
xlabel('Time Lag (ms)')


% FIG2 Sup
pbins = linspace(-7,-3,100);
pfd = histc(log10(rhm(:,:)),pbins,1);
figure,imagesc(pbins,fs,pfd'),axis xy
title('RHM Power Distribution')
xlabel('RHM Power (a.u.)')
ylabel('Frequency (Hz)')


%% FIG3 Segmentation and Features

Trial = MTATrial('jg05-20120317');
Stc = Trial.load('stc','hand_labeled_rev1'); 
figPath = '/storage/gravio/manuscripts/man2015-jgEd-MoCap/Figures/Figure_3';
pPad = [0,0];

%exPer = [26664,27100];
%exPer = [26000,26480];
%exPer = [25200,25580];
%exPer = [55200,60000];
exPer = [51549, 55145];
%exPer = [129470,139470];




xyz = Trial.load('xyz').filter(gtwin(.05,Trial.xyz.sampleRate));
ang = create(MTADang,Trial,xyz);
%stateColors = 'brcgym';


hfig = figure(38239385);clf
set(hfig,'position',[836   110   775   702]);
set(hfig,'paperposition',[0,0,775/100,702/100])
%% Fig:3:A Skeleton examples


axes('Position', [0.1300,0.5869,0.7750,0.4500]);hold on;

% $$$ set(gca,'CameraPositionMode', 'manual'                    ,...
% $$$ 	'YLim', [-200 200],...
% $$$         'XLim', [-300 400],...
% $$$ 	'ZLim', [0 300],...
% $$$ 	'CameraPosition',     [-1909.49 3535.01 1621.08],...
% $$$ 	'CameraTargetMode',   'manual'                    ,...
% $$$ 	'CameraTarget',       [50 0 150]                ,...
% $$$  	'CameraUpVectorMode', 'manual'                    ,...
% $$$  	'CameraUpVector',     [0 0 1]                   ,...
% $$$  	'CameraViewAngleMode','manual'                  ,...
% $$$  	'CameraViewAngle',    [6.88708]);
% $$$ daspect([1,1,1])
set(gca,'CameraPositionMode', 'manual'                    ,...
	'XLim', [-300 400],...
	'YLim', [-300 400],...
	'ZLim', [0 300],...
        'CameraPosition', [2050.7 4543.64 1748.25],...
	'CameraPositionMode','manual',...
	'CameraTarget',[50 50 150],...
	'CameraTargetMode','manual',...
	'CameraUpVector',[0 0 1],...
	'CameraUpVectorMode','manual',...
	'CameraViewAngle',[6.31812],...
	'CameraViewAngleMode','manual')
daspect([1,1,1]);
        
        
pMode = 'line';  %'surface';
xyz = Trial.load('xyz').filter(gtwin(.05,Trial.xyz.sampleRate));

plotSkeleton(xyz,exPer(1)+1500,pMode,[],[0,500]);% rear
plotSkeleton(xyz,exPer(1)+2000,pMode,[],[0,300]);% rear
plotSkeleton(xyz,exPer(1)+2300,pMode,[],[0,0]);% rear


zlim([0,300]);



% Fig:3:B - feature matrix
fet = fet_lgr(Trial);
axes('Position',[ 0.1300,0.3397,0.7750,0.2000])
ts = (1:fet.size(1))./fet.sampleRate;
per = round(exPer./xyz.sampleRate)+pPad;
ind = ts>per(1)&ts<per(2);
ufet = nunity(fet);
imc = imagesc(ts(ind),1:fet.size(2),ufet(ind,:)');
caxis([-2,2]);
xlim(round(exPer./xyz.sampleRate)+pPad)
Lines(round([exPer(1)+1500,exPer(1)+2000,exPer(1)+2300]./xyz.sampleRate),[],'k');
%xlabel('Time (s)')
flabels = {'speed SL'   ,...
           'speed SM'   ,...
           'speed HF'   ,...
           'height SL'  ,...
           'Z-diff SL_HF',...
           'pitch SL_PR',...
           'pitch SM_SU',...
           'dist SL_PR' ,...
           'dist PR_SM' ,...
           'dist SU_HB' ,...
           'dist SL_HB' ,...
           'av_SLSM_SMHF'};
set(gca,'YTickLabelMode','manual',...
        'YTickMode','manual',...
        'YTick',1:numel(flabels),...
        'YTickLabel',flabels);
set(gca,'XTickLabelMode','manual',...
        'XTickMode','manual',...
        'XTick',[],...
        'XTickLabel',{});
set(gca,'TickDir','out');


% Fig:3:C - Expert Labels
stateLabels = {'walk','rear','turn','groom','sit'};
stateColors = 'brgmc';
Stc = Trial.load('stc','hand_labeled_rev1');
axes('Position',[ 0.1300,0.2859,0.7750,0.0400])
plotSTC(Stc,1,'patch',stateLabels,stateColors);
xlim(round(exPer./xyz.sampleRate)+pPad)
set(gca,'YTickLabelMode','manual',...
        'YTickMode','manual',...
        'YTick',.5,...
        'YTickLabel',{'EXP'});
set(gca,'XTickLabelMode','manual',...
        'XTickMode','manual',...
        'XTick',[],...
        'XTickLabel',{});
set(gca,'TickDir','out');



% Fig:2:D - Empirical Labels
stateLabels = {'walk','rear'};
stateColors = 'br';
Stc = Trial.load('stc','auto_wbhr');
axes('Position', [0.1300,0.2426,0.7750,0.0400])
plotSTC(Stc,1,'patch',stateLabels,stateColors);
xlim(round(exPer./xyz.sampleRate)+pPad)
set(gca,'YTickLabelMode','manual',...
        'YTickMode','manual',...
        'YTick',.5,...
        'YTickLabel',{'EMP'});
set(gca,'XTickLabelMode','manual',...
        'XTickMode','manual',...
        'XTick',[],...
        'XTickLabel',{});
set(gca,'TickDir','out');




% Fig:3:E - LGR Model Labels
stateLabels = {'walk','rear','turn','groom','sit'};
stateColors = 'brgmc';
Stc = Trial.load('stc','LGR-hand_labeled_rev1-wrsnkm');
%Stc = Trial.load('stc','LGR_wrsnkm');
axes('Position',[ 0.1300,0.1983,0.7750,0.0400])
plotSTC(Stc,1,'patch',stateLabels,stateColors);
xlim(round(exPer./xyz.sampleRate)+pPad)
set(gca,'YTickLabelMode','manual',...
        'YTickMode','manual',...
        'YTick',.5,...
        'YTickLabel',{'LGR'});
set(gca,'XTickLabelMode','manual',...
        'XTickMode','manual',...
        'XTick',[],...
        'XTickLabel',{});
set(gca,'TickDir','out');


% Fig:3:F - LDA Model Labels
stateLabels = {'walk','rear','turn','groom','sit'};
stateColors = 'brgmc';
Stc = Trial.load('stc','LDA_hand_labeled_rev1-wrsnkm');
%Stc = Trial.load('stc','LDA_wrsnkm');
axes('Position',[ 0.1300,0.1540,0.7750,0.0400])
plotSTC(Stc,1,'patch',stateLabels,stateColors);
xlim(round(exPer./xyz.sampleRate)+pPad)
set(gca,'YTickLabelMode','manual',...
        'YTickMode','manual',...
        'YTick',.5,...
        'YTickLabel',{'LDA'});
set(gca,'XTickLabelMode','manual',...
        'XTickMode','manual',...
        'XTick',[],...
        'XTickLabel',{});
set(gca,'TickDir','out');


% Fig:3:G - LDA Model Classifier scores
axes('Position', [0.1300,0.0585,0.7750,0.0853])
[Stc,d_state] = bhv_lda(Trial,false,'display',false);
plot((1:size(d_state,1))/30,d_state)
xlim(round(exPer./xyz.sampleRate)+pPad)
ylim([-10,-2])
ylabel('A.U.');
xlabel('Time (s)');
set(gca,'TickDir','out');

saveas(hfig,fullfile(figPath,['Fig3A-F_complete_20150421.png']),'png');
saveas(hfig,fullfile(figPath,['Fig3A-F_complete_20150421.eps']),'epsc');





%% FIG4 Heirarichal Segmentation 

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
    b = [vxy(Trial.stc{sts(i)},1),vxy(Trial.stc{sts(i)},2)];
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


