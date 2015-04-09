%% new Figure2 
Trial = MTATrial('jg05-20120317');
figPath = '/storage/gravio/manuscripts/man2015-jgEd-MoCap/Figures/Figure_2';

%Features 
xyz = Trial.load('xyz').filter(gtwin(.1,Trial.xyz.sampleRate));

vl = vel(xyz,1:8,[1,2]);
vl.data = ButFilter(vl.data,3,4/(vl.sampleRate/2),'low');

[rhm,fs,ts] = fet_rhm(Trial,[],'wcsd');


ang = create(Trial.ang.copy,Trial,xyz);

sfet = Trial.ang.copy;
sfet.data = sum(abs([circ_dist(ang(:,1,3,1),ang(:,1,2,1)),...
                     circ_dist(ang(:,2,4,1),ang(:,2,3,1)),...
                     circ_dist(ang(:,3,5,1),ang(:,3,4,1)),...
                     circ_dist(ang(:,4,7,1),ang(:,4,5,1))]),2);

swg = Trial.ang.copy;
swg.data = circ_dist(ang(:,1,4,1),ang(:,2,4,1));


%% mutinfo crap
v = log10(vl.data(:,[1:4,7]));

edges = linspace(-2,2,64);
sbound = -30:30;
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




%% Fig Parameters
sts = 'rwnms';
stc = 'rcymg';
edgs    = {linspace(-2,2,75)};
edgs(2) = {linspace(-2,2,75)};
[edgs{:}] = get_histBinCenters(edgs);
[X,Y] = meshgrid(edgs{:});
ns = 6;

xpers = bsxfun(@plus,Trial.stc{'w',1}(1:3:end,1),[-15,15]);
xpers(xpers(:,1)<1,:) = [];
xpers((xpers(:,2)-Trial.sync(end))>0,:) = [];
s = 20;




hfig = figure(2);


% JPDF - Head/Body speed
subplot2(ns,4,[3,4],4);
b = log10([median(vl(nniz(vl),[1:2]),2),median(vl(nniz(vl),[5:8]),2)]);
hist2(b,edgs{1},edgs{2});
xlabel('log10 body speed (cm/s)');
ylabel('log10 head speed (cm/s)');
title('JPDF of log10 head and body speeds');

hold on,
for i = 1:numel(sts),
    b = log10([median(vl(Trial.stc{sts(i)},1:2),2),median(vl(Trial.stc{sts(i)},5:8),2)]);
    o = hist2(b,linspace(-1.5,2,75),linspace(-1.5,2,75));
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


ns = 6;
for s = 1:size(xpers,1),
    ind = round(xpers(s,:).*xyz.sampleRate);
    ind = ind(1):ind(2);
    i = 1;
    
    %figure(201) % SPEED head and body
    subplot2(ns,4,i,[1:3]);i=i+1;cla
    plot(ind/vl.sampleRate,[median(vl(ind,1:2),2),median(vl(ind,5:8),2)]),axis tight
    title('xy speed of head and body')
    %xlabel('Time (s)')
    ylabel('Speed (cm/s)');
    ylim([0,80])
    set(gca,'XTickLabelMode','manual');
    set(gca,'XTickLabel',{});
    
    %figure(202) % DIRECTION head and body
    subplot2(ns,4,i,[1:3]);i=i+1;cla
    plot(ind/vl.sampleRate,[circ_dist(ang(ind,1,4,1),ang(ind,2,4,1)),circ_dist(ang(ind,5,7,1),ang(ind,1,4,1))])
    axis tight
    %plot(ind/vl.sampleRate,[ang(ind,1,3,1),ang(ind,4,7,1)]),axis tight
    title('diff ang(head,body) and spine waggle')
    %xlabel('Time (s)')
    ylabel('Circular difference (radians)');
    ylim([-pi,pi])
    set(gca,'XTickLabelMode','manual');
    set(gca,'XTickLabel',{});
    
    %figure(203)%  PITCH SLPR and SMSU
    subplot2(ns,4,i,[1:3]);i=i+1;cla
    plot(ind/vl.sampleRate,[ang(ind,3,4,2),ang(ind,5,7,2)]),axis tight
    title('Pitch of Body and Head')
    %xlabel('Time (s)')
    ylabel('Pitch (radians)');
    ylim([-pi/2,pi/2])
    set(gca,'XTickLabelMode','manual');
    set(gca,'XTickLabel',{});
    
    %figure(204)
    subplot2(ns,4,i,[1:3]);i=i+1;cla
    plot(ind/vl.sampleRate,[ang(ind,4,5,3)]),axis tight%,ang(ind,4,7,3)])
    title('Intermarker Distance of head and body')
    %xlabel('Time (s)')
    ylabel('Distance (mm)');
    set(gca,'XTickLabelMode','manual');
    set(gca,'XTickLabel',{});

    
    %figure(205)
    subplot2(ns,4,i,[1:3]);i=i+1;cla
    plotSTC(Trial.stc,1,'patch',{'rear','walk','turn','groom'},'rbgm');
    title('hand labeling')
    xlim(xpers(s,:));
    ylim([0,1])
    set(gca,'XTickLabelMode','manual');
    set(gca,'XTickLabel',{});
   
    
    %figure(206) RHM Rhythmic Head Motion
    subplot2(ns,4,i,[1:3]);i=i+1;cla
    sind = round(xpers(s,:).*rhm.sampleRate);
    sind = sind(1):sind(2);
    title('RHM Rhythmic Head Motion')
    imagesc(ts(sind),fs,log10(rhm(sind,:))'),axis xy ,caxis([-5.2,-2.5])
    xlabel('Time (s)')
    
    
    saveas(hfig,fullfile(figPath,['Fig2-Features-sample_' num2str(s) '_' Trial.filebase '.png']),'png');
    saveas(hfig,fullfile(figPath,['Fig2-Features-sample_' num2str(s) '_' Trial.filebase '.eps']),'epsc');


end
mind = Trial.stc{'m'};
xe = linspace(0.4,pi/2,100);
ye = linspace(-pi/4,pi/2,100);
figure,
hist2([ang(mind,1,2,2),ang(mind,3,4,2)],xe,ye)
xlabel('SLPR pitch')
ylabel('SMSU pitch')
caxis([0,500])


figure,
plot(ButFilter(circ_dist(ang(:,1,4,1),ang(:,2,4,1)),3,1/(ang.sampleRate/2),'high')),
hold on,
plot(circ_dist(ang(:,1,4,1),ang(:,2,4,1)),'r'),
Lines(Trial.stc{'w'}(:),[],'k');
Lines(Trial.stc{'n'}(:),[],'g');

plot(circ_dist(ang(:,5,7,1),ang(:,1,4,1)),'c')



figure,plot(circ_dist(ang(:,5,7,1),ang(:,1,4,1)))
Lines(Trial.stc{'w'}(:),[],'k');
Lines(Trial.stc{'n'}(:),[],'g');

ned = linspace(50,300,100);

figure
sts = 'arwnm';
nsts = numel(sts);
for i = 1:nsts
subplot(nsts,1,i);bar(ned,histc(ang(Trial.stc{sts(i)},1,7,3),ned),'histc'),title(Trial.stc{sts(i)}.label)
end

ned = linspace(0,8,100);
figure
sts = 'arwnm';
nsts = numel(sts);
for i = 1:nsts
subplot(nsts,1,i);bar(ned,histc(sfet(Trial.stc{sts(i)}),ned),'histc'),title(Trial.stc{sts(i)}.label)
end

%% new figure 3
Trial = MTATrial('jg05-20120317');
Stc = Trial.load('stc','hand_labeled_rev1'); 
figPath = '/storage/gravio/manuscripts/man2015-jgEd-MoCap/Figures/Figure_3';

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
ang = create(Trial.ang.copy,Trial,xyz);
vxy.data = [ang(:,1,2,2),xyz(:,7,3)];
vxy.data = ButFilter(vxy.data,3,4/(xyz.sampleRate/2),'low');
tag = 'Head_HeightVSlower_spine_pitch';
v1 = vxy.copy;
v1.data = vxy(:,1);
v2 = vxy.copy;
v2.data = log10(vxy(:,2));
bhv_JPDF(Trial,v1,v2,70,70,[-.7,2],[-.7,2],...
         'Pitch (Lower Spine) radians',...
         'Height (Head Front) log10(mm)',{'a-r','r'},...
         tag)


sts = 'w';
stc = 'r';
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



%Segmentation JPDF walk
vxy = xyz.vel([1,7],[1,2]);
vxy.data = ButFilter(vxy.data,3,4/(xyz.sampleRate/2),'low');
tag = 'speed_SpineL_vs_HeadF';
v1 = vxy.copy;
v1.data = log10(vxy(:,1));
v2 = vxy.copy;
v2.data = log10(vxy(:,2));
bhv_JPDF(Trial,v1,v2,70,70,[-.7,2],[-.7,2],...
         'Speed (Lower Spine) log10(mm)',...
         'Speed (Head Front) log10(mm)',{'a-r','w','a-w-r'},...
         tag)


sts = 'w';
stc = 'r';
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

