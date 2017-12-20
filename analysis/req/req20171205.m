% req20171205 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: Determine the respiration frequency distribution during 
%               immobility
%  Bugs: NA


% Figure Settings ----------------------------------------------------------------------
OwnDir = '/storage/gravio/nextcloud/';
FigDir = 'MjgER2016/figures/supplementary';
% --------------------------------------------------------------------------------------

pitchReferenceTrial = 'Ed05-20140529.ont.all';

Trial = MTATrial.validate('Ed10-20140813.cof.all');

Trial = MTATrial.validate('Ed10-20140816.cof.all');
Trial = MTATrial.validate('Ed10-20140817.cof.gnd');
ncpChannel = 65;

Trial = MTATrial.validate('Ed05-20140528.cof.all');
Trial = MTATrial.validate('Ed05-20140529.ont.all');
ncpChannel = 2;

stc   = Trial.load('stc','msnn_ppsvd_raux');
xyz   = preproc_xyz(Trial,{'SPLINE_SPINE_HEAD_EQI'});
xyz.filter('ButFilter',3,2.4,'low');


pch = fet_HB_pitch(Trial);
map_to_reference_session(pch,Trial,pitchReferenceTrial);

% SPEED 
vh = xyz.vel({'spine_lower','head_front'},[1,2]);
vh.data = log10(abs(vh.data));
vh.data(~nniz(vh(:))) = nan;

% RHYTHMIC HEAD MOTION 
rhm = fet_rhm(Trial);
rhm.resample(xyz);
rhm.data(~nniz(rhm.data(:))) = nan;
% NORMALIZE rhm feature
rhm.data = nunity(rhm.data);

[rhms,fs,ts] = fet_rhm(Trial,[],'mtchglong',true);
% figure,imagesc(ts,fs,log10(rhms.data'));colormap('jet');axis('xy');
rhmf = rhms.copy();
fsInd = 2<fs&fs<13;
fsSub = fs(fsInd);
[~,rhmf.data] = max(log10(rhms(:,fsInd)),[],2);
rhmf.data = fsSub(rhmf.data);
rhmf.filter('ButFilter',3,2,'low');
rhmf.resample(xyz);
rhmp = rhms.copy;
rhmp.data = median(log10(rhmp(:,fs<6&fs<13)),2);
rhmp.resample(xyz);


% RESPIRATION 
ncp = fet_ncp(Trial,rhm,'mta',ncpChannel);
ncp.data = clip(ncp.data,-5e4,5e4);
%figure,plot(ncp.data);
ncp.filter('RectFilter',5,3);
ncp.filter('ButFilter',3,[1,20],'bandpass');
%hold on;plot(ncp.data);
ncp.data = nunity(ncp.data);
%hold on;plot(ncp.data*1000);
% FIND timepoints of inhalation
ncpPeaks = LocalMinima(ncp.data,8,-0.5);
ncpPeaks(ncpPeaks>size(ncp,1)) = [];
% COMPUTE instantaneous respiration rate (Hz)
ncpFreq = 1./((ncpPeaks-circshift(ncpPeaks,1)+circshift(ncpPeaks,-1)-ncpPeaks)/2/ncp.sampleRate);
% SMOOTH respiration rate (Hz)
ncpFreq = median(GetSegs(ncpFreq,circshift([1:size(ncpFreq,1)]',5),11));


pind = stc{'x+p+s',ncp.sampleRate};
pind = stc{'p',ncp.sampleRate};



ind = ncpPeaks;
sind = WithinRanges(ncpPeaks,pind.data);
ind = ncpPeaks(sind);

cs = jet(25);
edc = linspace(1,10,25);
cbin = discretize(ncpFreq(sind)',edc);
cbin(isnan(cbin)) = 25;
hfig = figure();
csc = cs(cbin,:);
csc(cbin==26,:) = repmat([1,1,1],[sum(cbin==26),1]);
scatter(pch(ind,3),rhmp(ind),5,csc,'filled')
xlim([-pi/2,pi/2]);
ylim([-10,-5]);
colormap(jet);
caxis([1,10]);
colorbar();

FigName = ['respiration_headpitch_X_rhmPow_C_ncpFreq',Trial.filebase];

print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));



cs = jet(25);
edc = linspace(-pi/2,pi/2,25);
cbin = discretize(pch(ind,3)',edc);
cbin(isnan(cbin)) = 25;
hfig = figure();
csc = cs(cbin,:);
csc(cbin==26,:) = repmat([1,1,1],[sum(cbin==26),1]);

figure
scatter(ncpFreq(sind)'+randn([sum(sind),1])/10,rhmf(ind),5,pch(ind,3),'filled')
colormap('jet');

figure
scatter(ncpFreq(sind)'+randn([sum(sind),1])/10,rhmf(ind),5,rhmp(ind),'filled')
colormap('jet');

figure
scatter(ncpFreq(sind)'+randn([sum(sind),1])/10,rhmf(ind),5,pch(ind,3),'filled')
colormap('jet');

figure
scatter(pch(ind,3),rhmf(ind),5,rhmp(ind),'filled')
colormap('jet');
xlabel('pitch');
ylabel('rhm freq');


figure
scatter(vh(ind,2),rhmf(ind),5,rhmp(ind),'filled')
colormap('jet');


figure
scatter(vh(ind,2),rhmf(ind),5,pch(ind,3),'filled')
colormap('jet');


daspect([1,1,1]);
line([0,12],[0,12],'Color','m');


pper = stc{'p'};
plper = cast(stc{'q'},'TimeSeries');
phper = cast(stc{'j'},'TimeSeries');
phlr = [];
pdur = [];
for p = pper.data'
    pls = sum(plper(p'));
    phs = sum(phper(p'));
    phlr(end+1) = pls/(phs+pls);
    pdur(end+1) = diff(p)./pper.sampleRate;
end

figure();
bar(0:0.05:1,histc(phlr,0:0.05:1),'histc');

figure();
plot(phlr'+randn([numel(phlr),1])/100,log10(pdur),'.');

% FIGURE head speed ------------------a----------------------------------------------
hfig = figure();
hfig.Position(3:4) = [1000,900];
hfig.PaperPositionMode = 'auto';
tind = stc{'t',ncp.sampleRate};
edx = linspace(-3,3,25);
edy = linspace(0,16,25);

% RESPIRATION during walk
pind = stc{'x',ncp.sampleRate};
subplot(331);
ind = WithinRanges(ncpPeaks,pind.data);
hist2([vh(ncpPeaks(ind),2),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);
title({'Respiration Frequency','during Walk'});
grid('on');
subplot(332);
ind = pind&tind.data;
ind = WithinRanges(ncpPeaks,ind.data);
hist2([vh(ncpPeaks(ind),2),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);
title({'Respiration Frequency','during Walk & Theta'});
grid('on');
subplot(333);
ind = pind-tind;
ind = WithinRanges(ncpPeaks,ind.data);
hist2([vh(ncpPeaks(ind),2),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);
title({'Respiration Frequency','during Walk - Theta'});
grid('on');

% RESPIRATION during pause
pind = stc{'p',ncp.sampleRate};
subplot(334);  
ind = WithinRanges(ncpPeaks,pind.data);
hist2([vh(ncpPeaks(ind),2),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);
title({'Respiration Frequency','during Pause'});
grid('on');
subplot(335);  
ind = pind&tind.data;
ind = WithinRanges(ncpPeaks,ind.data);
hist2([vh(ncpPeaks(ind),2),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);
title({'Respiration Frequency','during Pause & Theta'});
grid('on');
subplot(336);  
ind = pind-tind;
ind = WithinRanges(ncpPeaks,ind.data);
hist2([vh(ncpPeaks(ind),2),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);
title({'Respiration Frequency','during Pause - Theta'});
grid('on');

% RESPIRATION during sit
pind = stc{'s',ncp.sampleRate};
subplot(337);  
ind = WithinRanges(ncpPeaks,pind.data);
hist2([vh(ncpPeaks(ind),2),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);
title({'Respiration Frequency','during sit'});
grid('on');
subplot(338);  
ind = pind&tind.data;
ind = WithinRanges(ncpPeaks,ind.data);
hist2([vh(ncpPeaks(ind),2),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);;
title({'Respiration Frequency','during sit & Theta'});
grid('on');
subplot(339);
ind = pind-tind;
ind = WithinRanges(ncpPeaks,ind.data);
hist2([vh(ncpPeaks(ind),2),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);
title({'Respiration Frequency','during sit - Theta'});
grid('on');

% LABELS 
af(@(h)  xlabel(h,'log10 head speed'), findobj(hfig,'Type','Axes'));
af(@(h)  ylabel(h,'frequency'),        findobj(hfig,'Type','Axes'));

% SUPTITLE 
suptitle([Trial.filebase])

FigName = ['respiration_frequency_X_headSpeed',Trial.filebase];

print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

% END FIGURE head speed ----------------------------------------------------------------






% FIGURE body speed ----------------------------------------------------------------
hfig = figure();
hfig.Position(3:4) = [1000,900];
hfig.PaperPositionMode = 'auto';
tind = stc{'t',ncp.sampleRate};
edx = linspace(-3,3,25);
edy = linspace(0,16,25);

% RESPIRATION during walk
pind = stc{'x',ncp.sampleRate};
subplot(331);
ind = WithinRanges(ncpPeaks,pind.data);
hist2([vh(ncpPeaks(ind),1),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);
title({'Respiration Frequency','during Walk'});
grid('on');
subplot(332);
ind = pind&tind.data;
ind = WithinRanges(ncpPeaks,ind.data);
hist2([vh(ncpPeaks(ind),1),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);
title({'Respiration Frequency','during Walk & Theta'});
grid('on');
subplot(333);
ind = pind-tind;
ind = WithinRanges(ncpPeaks,ind.data);
hist2([vh(ncpPeaks(ind),1),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);
title({'Respiration Frequency','during Walk - Theta'});
grid('on');

% RESPIRATION during pause
pind = stc{'p',ncp.sampleRate};
subplot(334);  
ind = WithinRanges(ncpPeaks,pind.data);
hist2([vh(ncpPeaks(ind),1),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);
title({'Respiration Frequency','during Pause'});
grid('on');
subplot(335);  
ind = pind&tind.data;
ind = WithinRanges(ncpPeaks,ind.data);
hist2([vh(ncpPeaks(ind),1),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);
title({'Respiration Frequency','during Pause & Theta'});
grid('on');
subplot(336);  
ind = pind-tind;
ind = WithinRanges(ncpPeaks,ind.data);
hist2([vh(ncpPeaks(ind),1),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);
title({'Respiration Frequency','during Pause - Theta'});
grid('on');

% RESPIRATION during sit
pind = stc{'s',ncp.sampleRate};
subplot(337);  
ind = WithinRanges(ncpPeaks,pind.data);
hist2([vh(ncpPeaks(ind),1),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);
title({'Respiration Frequency','during sit'});
grid('on');
subplot(338);  
ind = pind&tind.data;
ind = WithinRanges(ncpPeaks,ind.data);
hist2([vh(ncpPeaks(ind),1),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);;
title({'Respiration Frequency','during sit & Theta'});
grid('on');
subplot(339);
ind = pind-tind;
ind = WithinRanges(ncpPeaks,ind.data);
hist2([vh(ncpPeaks(ind),1),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);
title({'Respiration Frequency','during sit - Theta'});
grid('on');

% LABELS 
af(@(h)  xlabel(h,'log10 body speed'), findobj(hfig,'Type','Axes'));
af(@(h)  ylabel(h,'frequency'),        findobj(hfig,'Type','Axes'));

% SUPTITLE 
suptitle([Trial.filebase])

FigName = ['respiration_frequency_X_bodySpeed',Trial.filebase];

print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

% END FIGURE body speed ----------------------------------------------------------------





% FIGURE head pitch ----------------------------------------------------------------
hfig = figure();
hfig.Position(3:4) = [1000,900];
hfig.PaperPositionMode = 'auto';
tind = stc{'t',ncp.sampleRate};
edx = linspace(-pi/2,pi/2,25);
edy = linspace(0,16,25);

% RESPIRATION during walk
pind = stc{'x',ncp.sampleRate};
subplot(331);
ind = WithinRanges(ncpPeaks,pind.data);
hist2([ang(ncpPeaks(ind),5,7,2),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);
title({'Respiration Frequency','during Walk'});
grid('on');
subplot(332);
ind = pind&tind.data;
ind = WithinRanges(ncpPeaks,ind.data);
hist2([ang(ncpPeaks(ind),5,7,2),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);
title({'Respiration Frequency','during Walk & Theta'});
grid('on');
subplot(333);
ind = pind-tind;
ind = WithinRanges(ncpPeaks,ind.data);
hist2([ang(ncpPeaks(ind),5,7,2),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);
title({'Respiration Frequency','during Walk - Theta'});
grid('on');

% RESPIRATION during pause
pind = stc{'p',ncp.sampleRate};
subplot(334);  
ind = WithinRanges(ncpPeaks,pind.data);
hist2([ang(ncpPeaks(ind),5,7,2),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);
title({'Respiration Frequency','during Pause'});
grid('on');
subplot(335);  
ind = pind&tind.data;
ind = WithinRanges(ncpPeaks,ind.data);
hist2([ang(ncpPeaks(ind),5,7,2),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);
title({'Respiration Frequency','during Pause & Theta'});
grid('on');
subplot(336);  
ind = pind-tind;
ind = WithinRanges(ncpPeaks,ind.data);
hist2([ang(ncpPeaks(ind),5,7,2),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);
title({'Respiration Frequency','during Pause - Theta'});
grid('on');

% RESPIRATION during sit
pind = stc{'s',ncp.sampleRate};
subplot(337);  
ind = WithinRanges(ncpPeaks,pind.data);
hist2([ang(ncpPeaks(ind),5,7,2),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);
title({'Respiration Frequency','during sit'});
grid('on');
subplot(338);  
ind = pind&tind.data;
ind = WithinRanges(ncpPeaks,ind.data);
hist2([ang(ncpPeaks(ind),5,7,2),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);;
title({'Respiration Frequency','during sit & Theta'});
grid('on');
subplot(339);
ind = pind-tind;
ind = WithinRanges(ncpPeaks,ind.data);
hist2([ang(ncpPeaks(ind),5,7,2),ncpFreq(ind)'+randn([sum(ind),1])],edx,edy);
title({'Respiration Frequency','during sit - Theta'});
grid('on');

% LABELS 
af(@(h)  xlabel(h,'head pitch (rad)'), findobj(hfig,'Type','Axes'));
af(@(h)  ylabel(h,'frequency'),        findobj(hfig,'Type','Axes'));

% SUPTITLE 
suptitle([Trial.filebase])

FigName = ['respiration_frequency_X_pitch',Trial.filebase];

print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

% END FIGURE head pitch ----------------------------------------------------------------



% CLUSTER 
[ys,fs,ts] = fet_ncp(Trial,[],'mtchglong',32,true);

hfig = figure();  
% Walk & Theta
pind = stc{'x',ncp.sampleRate};
ind = pind&tind.data;
ind = WithinRanges(ncpPeaks,ind.data);

subplot(231);
plot(vh(ncpPeaks(ind),1),ncpFreq(ind)'+randn([sum(ind),1]),'.');
cpnts = ClusterPP(gcf);
copyobj(get(gcf,'Children'),hfig);
set(


cind = double(ind);
cind(ind) = cpnts';

figure();hold('on');
for c = unique(cind)'
    if c==0,continue,end;
    mys = log10(nanmean(ys(round(ncpPeaks(cind==c)./ncp.sampleRate.*ys.sampleRate),:)));    
    errorbar(fs,mys,...
             log10(prctile(ys(round(ncpPeaks(cind==c)./ncp.sampleRate.*ys.sampleRate),:),5))-mys,...
             log10(prctile(ys(round(ncpPeaks(cind==c)./ncp.sampleRate.*ys.sampleRate),:),95))-mys);
end
legend({'pauseHpitch','pauseLpitch','pauseLResp'})


pind = stc{'p',ncp.sampleRate};
%ind = pind&tind.data;
ind = pind-tind.data;
ind = WithinRanges(ncpPeaks,ind.data);



% Walk - Theta
pind = stc{'x',ncp.sampleRate};
ind = pind-tind.data;
ind = WithinRanges(ncpPeaks,ind.data);

figure();  plot(vh(ncpPeaks(ind),1),ncpFreq(ind)'+randn([sum(ind),1]),'.');
cpnts = ClusterPP(gcf);
copyax('on');

copyax('off');
cind = double(ind);
cind(ind) = cpnts';

subplot(231);hold('on');
for c = unique(cind)'
    if c==0,continue,end;
    mys = log10(nanmean(ys(round(ncpPeaks(cind==c)./ncp.sampleRate.*ys.sampleRate),:)));    
    errorbar(fs,mys,...
             log10(prctile(ys(round(ncpPeaks(cind==c)./ncp.sampleRate.*ys.sampleRate),:),5))-mys,...
             log10(prctile(ys(round(ncpPeaks(cind==c)./ncp.sampleRate.*ys.sampleRate),:),95))-mys);
end
legend({'pauseHpitch','pauseLpitch','pauseLResp'})


pind = stc{'p',ncp.sampleRate};
%ind = pind&tind.data;
ind = pind-tind.data;
ind = WithinRanges(ncpPeaks,ind.data);
