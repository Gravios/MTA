Trial = MTATrial.validate('Ed05-20140529.ont.all');
Trial.load('stc','hand_labeled_rev1_Ed');
%xyz = Trial.load('xyz','trb');
xyz = Trial.load('xyz','seh');
rb = Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'});
hcom = xyz.com(rb);
xyz.addMarker('fhcom',[128,255,128],{{'head_back','head_front',[0,0,1]}},...
               ButFilter(hcom,3,[2]./(Trial.xyz.sampleRate/2),'low'));
xyz.addMarker('hcom',[128,255,128],{{'head_back','head_front',[0,0,1]}},...
               hcom);
xyz.addMarker('fsu',[128,255,128],{{'head_back','head_front',[0,0,1]}},...
               ButFilter(xyz(:,4,:),3,[2]./(Trial.xyz.sampleRate/2),'low'));

ang = create(MTADang,Trial,xyz);
ncp = fet_ncp(Trial,xyz);
rhm = fet_rhm(Trial,xyz);
pth = ang.copy;
pth.data = ang(:,5,7,2);
pth.filter('ButFilter',3,[2,50],'bandpass');
pth.label = 'headPitch';
shd = ang.copy;
shd.data = ang(:,'fsu','hcom',3);
shd.filter('ButFilter',3,[2,50],'bandpass');
shd.label = 'distSpHe';
bth = ang.copy;
bth.data = ang(:,'fsu','hcom',2);
bth.label = 'neckPitch';

state = 'a-n-m-k-s';
fhn = 238499;

fhn = fhn+1;
bhv_mean_coherence(Trial,ncp,rhm,state,'figHnum',fhn); 
print(gcf,'-depsc2',fullfile(getenv('PROJECT'),'manuscripts/man2015-jgEd-MoCap/Figures/Figure_5',...
                     ['fig_ncpX' rhm.label '_coherence.eps']))

fhn = fhn+1;
bhv_mean_coherence(Trial,ncp,pth,state,'figHnum',fhn); 
print(gcf,'-depsc2',fullfile(getenv('PROJECT'),'manuscripts/man2015-jgEd-MoCap/Figures/Figure_5',...
                     ['fig_ncpX' pth.label '_coherence.eps']))

fhn = fhn+1;
bhv_mean_coherence(Trial,ncp,shd,state,'figHnum',fhn); 
print(gcf,'-depsc2',fullfile(getenv('PROJECT'),'manuscripts/man2015-jgEd-MoCap/Figures/Figure_5',...
                     ['fig_ncpX' shd.label '_coherence.eps']))

fhn = fhn+1;
bhv_mean_coherence(Trial,ncp,bth,state,'figHnum',fhn);  
print(gcf,'-depsc2',fullfile(getenv('PROJECT'),'manuscripts/man2015-jgEd-MoCap/Figures/Figure_5',...
                     ['fig_ncpX' bth.label '_coherence.eps']))




pncp = ncp.copy;
ppth = pth.copy;
pshd = shd.copy;
prhm = rhm.copy;

pncp = pncp.phase([5,15]);
ppth = ppth.phase([5,15]);
pshd = pshd.phase([5,15]);
prhm = prhm.phase([5,15]);


iang = ang(:,5,7,2)<-0.5;
ind = Trial.stc{'w+p'};
ind.cast('TimeSeries');
ind = ind.data&iang;

set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)

hfig = figure(2016062301);
set(hfig,'PaperPositionMode','auto');
set(hfig,'units','centimeters');
set(hfig,'Position',[0,0,20,8])


i= 1
axes('Units',   'centimeters',...
     'Position',[(i-1)*5+1,2,3,3]);
hist2([ppth(ind),pshd(ind)],50,50)
xlabel('head pitch phase [5-15Hz]');
ylabel('spine-head distance phase [5-15Hz]');
title('Trans VS Rot Phase')

i = i+1;
axes('Units',   'centimeters',...
     'Position',[(i-1)*5+1,2,3,3]);
hist2([pncp(ind),pshd(ind)],50,50)
title('ncp VS HS dist')
xlabel('ncp phase [5-15Hz]');
ylabel('spine-head distance phase [5-15Hz]');

i = i+1;
axes('Units',   'centimeters',...
     'Position',[(i-1)*5+1,2,3,3]);
hist2([pncp(ind),ppth(ind)],50,50)
xlabel('ncp phase [5-15Hz]');
ylabel('spine-head distance phase [5-15Hz]');
title('ncp Vs head pitch phase')

i = i+1;
axes('Units',   'centimeters',...
     'Position',[(i-1)*5+1,2,3,3]);
hist2([pncp(ind),prhm(ind)],50,50)
title('ncp Vs rhm phase')
xlabel('ncp phase [5-15Hz]');
ylabel('spine-head distance phase [5-15Hz]');

print(hfig,'-depsc2',fullfile(getenv('PROJECT'),'manuscripts/man2015-jgEd-MoCap/Figures/Figure_5',...
                     ['fig_ncpXrhm_phase.eps']))



subplot(1,n,i);i=i+1;
hist2([pncp(ind),circ_dist(ppth(ind),pshd(ind))],50,50)
title('pitch')



cxyz = xyz.copy;
cxyz.data = bsxfun(@minus,xyz.data,xyz(:,'fsu',:));
cxyz.filter('ButFilter',3,[5,15],'bandpass');

pcxyz = cxyz.copy;
pcxyz.data = sq(pcxyz(:,'hcom',:));
pcxyz = pcxyz.phase([5,15]);

subplot(1,n,i);i=i+1;
hist2([pncp(ind),circ_dist(ppth(ind),pshd(ind))],50,50)
title('pitch')

figure,sp=[];
for i = 1:3,
    sp(i) = subplot(1,3,i);
    hist2([pncp(ind),ButFilter(cxyz(ind,'hcom',i),3,[5,15]./(cxyz.sampleRate.*.5),'bandpass')],50,linspace(-1,1,50));
end


pcxyz = cxyz.copy;
pcxyz.data = sq(pcxyz(:,'hcom',:));
pcxyz = pcxyz.phase([5,15]);
figure,sp=[];
for i = 1:3,
    sp(i) = subplot(1,3,i);
    hist2([pncp(ind),pcxyz(ind,i)],50,50);
end
