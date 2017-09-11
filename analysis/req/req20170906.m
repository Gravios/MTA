Trial = MTATrial.validate('Ed01-20140707.cof.all');
stc = Trial.load('stc','hand_labeled_rev2_jg');
Trial = MTATrial.validate('Ed05-20140529.ont.all');
stc = Trial.load('stc','hand_labeled_rev1_Ed');
Trial = MTATrial.validate('Ed05-20140528.cof.all');
stc = Trial.load('stc','msnn_ppsvd');

Trial = MTATrial.validate('jg05-20120317.cof.all');
stc = Trial.load('stc','hand_labeled_rev3_jg');

OwnDir = '/storage/gravio/nextcloud/';
FigDir = 'Shared/Behavior Paper/Figures/Figure_4/parts';



xyz = preproc_xyz(Trial);
xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},ButFilter(xyz(:,'hcom',:),3,[5]./(xyz.sampleRate/2),'low'));
xyz.addMarker('flhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},ButFilter(xyz(:,'hcom',:),3,[4]./(xyz.sampleRate/2),'low'));
xyz.addMarker('fvlhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},ButFilter(xyz(:,'hcom',:),3,[0.2]./(xyz.sampleRate/2),'low'));
nz = -cross(xyz(:,'head_back',:)-xyz(:,'hcom',:),xyz(:,'head_left',:)-xyz(:,'hcom',:));
nz = bsxfun(@rdivide,nz,sqrt(sum((nz).^2,3)));
nm = nz.*20+xyz(:,'hcom',:);
xyz.addMarker('htx',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);

xyz.addMarker('flbcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},ButFilter(xyz(:,'bcom',:),3,[1.2]./(xyz.sampleRate/2),'low'));

ang = create(MTADang,Trial,xyz);



vxy = xyz.vel([],[1,2]);
vxy.data(vxy.data<1e-3) = 1e-3;
vxy.data = log10(vxy.data);
ang = create(MTADang,Trial,xyz);
features = fet_bref(Trial);
features.filter('ButFilter',3,[1.2,12],'bandpass');
%ncp = fet_ncp(Trial);
rhm = fet_rhm(Trial);
rhbm = rhm.copy();
rhbm.data = features(:,24)-mean(features(:,[22]),2);
rhbm.filter('ButFilter',3,[5,14],'bandpass');
rpm = rhm.copy();
rpm.data = ang(:,'head_back','head_front',2);
rpm.filter('ButFilter',3,[5,14],'bandpass');

rrm = rhm.copy();
rrm.data = nunity(circ_dist(ang(:,'flhcom','head_right',1),ang(:,'flhcom','head_front',1)));
rrm.filter('ButFilter',3,[4,12],'bandpass');
rtm = rhm.copy();
rtm.data = nunity(circ_dist(ang(:,'flhcom','htx',2),ang(:,'flhcom','head_front',2)));
rtm.filter('ButFilter',3,[4,12],'bandpass');

rbm = rhm.copy();
rbm.data = nunity(circ_dist(ang(:,'flbcom','spine_lower',2),ang(:,'flbcom','spine_middle',2)));
rbm.filter('ButFilter',3,[1.2,12],'bandpass');

rsm = rhm.copy();
rsm.data = nunity(circ_dist(ang(:,'flbcom','spine_upper',2),ang(:,'flbcom','spine_middle',2)));
rsm.filter('ButFilter',3,[1.2,12],'bandpass');

hfig = figure();
subplot2(3,1,1:2,1);
plot([rbm.data*3,rsm.data*3,nunity(features(:,[17,19]))+5,nunity(features(:,[18,20,22]))-5]);
subplot2(3,1,3,1);
plotSTC(stc);
linkaxes(findobj(hfig,'Type','axes'),'x');



freqRange = [4,12];
prhbm = rhbm.phase(freqRange);
prhm = rhm.phase(freqRange);
%pncp = ncp.phase(freqRange);
prpm = rpm.phase(freqRange);
prrm = rrm.phase(freqRange);
prtm = rtm.phase(freqRange);

sper = stc{'w+p+n'};

ncpMins = LocalMinima(ncp.data,5,-2000);
ncpMins(~WithinRanges(ncpMins,sper.data)) = [];


figure,plot([nunity(ncp.data),[0;diff(rrm.data)]*100]);
figure,plot([nunity(ncp.data),[0;diff(rtm.data)]*100,nunity(rhm.data)])

figure,
subplot(241);rose(circ_dist(prhm(ncpMins),pncp(ncpMins)),100)
subplot(242);rose(circ_dist(prhbm(ncpMins),pncp(ncpMins)),100)
subplot(243);rose(circ_dist(prpm(ncpMins),pncp(ncpMins)),100)
subplot(244);rose(circ_dist(prrm(ncpMins),pncp(ncpMins)),100)
subplot(245);rose(circ_dist(prtm(ncpMins),pncp(ncpMins)),100)


figure,
hist2([circ_dist(prtm(ncpMins),pncp(ncpMins)),circ_dist(prrm(ncpMins),pncp(ncpMins))],linspace(-pi,pi,30),linspace(-pi,pi,30))

figure,
hist2([circ_dist(prtm(ncpMins),pncp(ncpMins)),circ_dist(prrm(ncpMins),pncp(ncpMins))],linspace(-pi,pi,30),linspace(-pi,pi,30))

figure,
hist2([circ_dist(prtm(ncpMins),pncp(ncpMins)),circ_dist(prhm(ncpMins),pncp(ncpMins))],linspace(-pi,pi,30),linspace(-pi,pi,30))

figure,
hist2([circ_dist(prpm(ncpMins),pncp(ncpMins)),circ_dist(prtm(ncpMins),pncp(ncpMins))],linspace(-pi,pi,30),linspace(-pi,pi,30))



figure,
hist2([circ_dist(prrm(ncpMins),pncp(ncpMins)),circ_dist(prhm(ncpMins),pncp(ncpMins))],linspace(-pi,pi,30),linspace(-pi,pi,30))


figure,
subplot(221);hist2([circ_dist(prhbm(ncpMins),pncp(ncpMins)),ang(ncpMins,'spine_lower','spine_upper',2)],30,linspace(0.1,.5,30))

subplot(222);hist2([circ_dist(prhm(ncpMins),pncp(ncpMins)),ang(ncpMins,'spine_lower','spine_upper',2)],30,linspace(0.1,.5,30))

subplot(223);hist2([circ_dist(prhbm(ncpMins),pncp(ncpMins)),ang(ncpMins,'spine_upper','hcom',2)],30,linspace(-0.8,.8,30))

subplot(224);hist2([circ_dist(prhm(ncpMins),pncp(ncpMins)),ang(ncpMins,'spine_upper','hcom',2)],30,linspace(-0.8,.8,30))


figure,
subplot(221);hist2([circ_dist(prhbm(ncpMins),pncp(ncpMins)),vxy(ncpMins,1)],30,linspace(-0.5,2,30))


ncpThresh = -1000;
ncpThresh = -2000;
[ncpTroughs,ncpTroughsVals] = LocalMinima(ncp.data,5,ncpThresh);
ncpTroughs(~WithinRanges(ncpTroughs,sper.data)) = [];
[ncpPeaks,ncpPeaksVals] = LocalMinima(-ncp.data,5,ncpThresh);
ncpPeaks(~WithinRanges(ncpPeaks,sper.data)) = [];
ncpPeaksVals = -ncpPeaksVals;

[ncpTroughs,ncpTroughsVals] = LocalMinima(rhm.data,5,-.2);
ncpTroughs(~WithinRanges(ncpTroughs,sper.data)) = [];
[ncpPeaks,ncpPeaksVals] = LocalMinima(-rhm.data,5,-0.2);
ncpPeaks(~WithinRanges(ncpPeaks,sper.data)) = [];
ncpPeaksVals = -ncpPeaksVals;

nbins = 30;
dp = [0;diff(ang(:,'head_back','head_front',2))];
dp = [0;diff(ang(:,'head_back','head_front',2))];
dp = [rhm.data];
dp = [ncp.data];
dp = [vxy(:,1)];
dp = [ang(:,1,4,2)];                         dedge = linspace(0,0.5,nbins);
dp = [ang(:,'spine_upper','hcom',2)];        dedge = linspace(-1.2,0.6,nbins);
dp = [ang(:,'head_back','head_front',2)];    dedge = linspace(-1.2,0.6,nbins);


dedge = linspace(-.05,.05,30);
dedge = linspace(-1e4,1e4,nbins);
dedge = linspace(-1,1,30);
dedge = linspace(-1,2,30); % vxy
dedge = linspace(-1,2,30); % vxy
dedge = linspace(-pi/2,pi/2,nbins);


bpedge = linspace(-pi,pi,nbins);


figure,
subplot(221);hist2([circ_dist(prrm(ncpTroughs),pncp(ncpTroughs)),dp(ncpTroughs,1)],pedge,dedge)
subplot(222);hist2([circ_dist(prrm(ncpPeaks),  pncp(ncpPeaks)  ),dp(ncpPeaks,1)],pedge,dedge)
subplot(223);hist2([circ_dist(prtm(ncpTroughs),pncp(ncpTroughs)),dp(ncpTroughs,1)],pedge,dedge)
subplot(224);hist2([circ_dist(prtm(ncpPeaks),  pncp(ncpPeaks)  ),dp(ncpPeaks,1)],pedge,dedge)

subplot(221);hist2([circ_dist(prhm(ncpTroughs),pncp(ncpTroughs)),dp(ncpTroughs,1)],pedge,dedge)
subplot(222);hist2([circ_dist(prhm(ncpPeaks),  pncp(ncpPeaks)  ),dp(ncpPeaks,1)],pedge,dedge)
subplot(223);hist2([circ_dist(prtm(ncpTroughs),pncp(ncpTroughs)),dp(ncpTroughs,1)],pedge,dedge)
subplot(224);hist2([circ_dist(prtm(ncpPeaks),  pncp(ncpPeaks)  ),dp(ncpPeaks,1)],pedge,dedge)

figure,
subplot(221);hist2([circ_dist(prtm(ncpTroughs),prhm(ncpTroughs)),dp(ncpTroughs,1)],pedge,dedge)
subplot(222);hist2([circ_dist(prtm(ncpPeaks),  prhm(ncpPeaks)  ),dp(ncpPeaks,1)],pedge,dedge)
subplot(223);hist2([circ_dist(prpm(ncpTroughs),prhm(ncpTroughs)),dp(ncpTroughs,1)],pedge,dedge)
subplot(224);hist2([circ_dist(prpm(ncpPeaks),  prhm(ncpPeaks)  ),dp(ncpPeaks,1)],pedge,dedge)
ForAllSubplots('caxis([0,40])')



Trials = af(@(t)  MTATrial.validate(t),  get_session_list('hand_labeled'));

[out,xedges,yedges,labels] = cf(@(t)  bhv_rhm_rrm(t),  Trials);
hfig = figure();
for s = 1:numel(Trials);
for i = 1:4,
    subplot2(numel(Trials),4,s,i);
    imagesc(xedges{s}{i},yedges{s}{i},out{s}{i}');
    axis('xy');
    title(labels{s}{i});
end
end

FigName = ['bvh_rhm_rrm-hand_labeled-individual'];
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));



aout =cat(1,out{:});

rout = {};
for i = 1:4,
    sout    = aout(:,i);
    rout{i} = sum(cat(3,sout{:}),3);
end


hfig = figure();
for i = 1:4,
    subplot(2,2,i);
    imagesc(xedges{1}{i},yedges{1}{i},rout{i}');
    axis('xy');
    title(labels{1}{i});
    caxis([0,200]);
end

FigName = ['bvh_rhm_rrm-hand_labeled-all'];
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));
