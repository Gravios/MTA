%% FIG4: Examples of feature dynamics and head body independence ---|
% A  SPEED head and body                                            |
% B  DIRECTION head and body                                        |
% C  PITCH SLPR and SMSU                                            |
% D  Intermarker Distance of head and upper spine                   |
% E  Hand Labels                                                    |
% F  RHM Rhythmic Head Motion                                       |
% G  Mutual information between marker speeds                       |
% H  Time lag of maximum mutual information (ms)                    |
% __________________________________________________________________|


MTAstartup('man_jgEd_2015');
Trial = MTATrial('jg05-20120317');
figPath = fullfile(Trial.path.project,'figures');
anlPath = fullfile(Trial.path.project,'analysis');
figPath = '/storage/gravio/manuscripts/man2015-jgEd-MoCap/Figures/Figure_2';

xyz = Trial.load('xyz');

% SPD Low Pass Filtered 4Hz
vl = xyz.copy;
vl = vl.vel(,1:8,[1,2]);
vl.filter('ButFilter',3,4);

% RHM Spectrum 1-30Hz
dsp.nFFT= 2^8;
dsp.Fs = 119.881035;
dsp.WinLength = 2^7;
dsp.nOverlap = 2^7.*0.8125;
dsp.FreqRange = [1 30];
[rhm,fs,ts] = fet_rhm(Trial,[],'mtchglong','defspec',dsp);

% XYZ Low Pass Filtered 4Hz
xyz.data = ButFilter(xyz.data,3,[4]/(xyz.sampleRate/2),'low');

% ANG From Low Pass Filtered (4Hz) XYZ data
ang = create(Trial.ang.copy,Trial,xyz);


sfet = Trial.ang.copy;
% $$$ sfet.data = sum(abs([circ_dist(ang(:,1,3,1),ang(:,1,2,1)),...
% $$$                      circ_dist(ang(:,2,4,1),ang(:,2,3,1)),...
% $$$                      circ_dist(ang(:,3,5,1),ang(:,3,4,1)),...
% $$$                      circ_dist(ang(:,4,7,1),ang(:,4,5,1))]),2);
sfet.data = -nunity(circ_dist(ang(:,5,7,1),ang(:,1,4,1)));

swg = Trial.ang.copy;
swg.data = nunity(circ_dist(ang(:,1,4,1),ang(:,2,4,1)));


%% mutinfo stuff
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

    
    %figure(204) Intermarker Distance of head and upper spine
    subplot2(ns,4,i,[1:3]);i=i+1;cla
    plot(ind/vl.sampleRate,[ang(ind,4,5,3)]),axis tight%,ang(ind,4,7,3)])
    title('Intermarker Distance of head and upper spine')
    ylabel('Distance (mm)');
    set(gca,'XTickLabelMode','manual');
    set(gca,'XTickLabel',{});
    set(gca,'TickDir','out');
    

    %figure(205) Hand labeles
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


% FIG3 Sup
pbins = linspace(-7,-3,100);
pfd = histc(log10(rhm(:,:)),pbins,1);
figure,imagesc(pbins,fs,pfd'),axis xy
title('RHM Power Distribution')
xlabel('RHM Power (a.u.)')
ylabel('Frequency (Hz)')


% FIG3 Sup

Trial = MTATrial('jg03-20110501');
Trial = MTATrial('Ed05-20140529','all','ont');
xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,20);
ang = create(MTADang,Trial,xyz);
wang = Trial.xyz.copy;
wang.data = diff(circ_dist(ang(:,1,3,1),ang(:,2,4,1)).*5);
wang.filter('ButFilter',3,20);
figure,plot(diff(xyz(:,1,3))-diff(ang(:,1,3,2)).*50);

Lines(Trial.stc{'w'}(:),[],'r');

