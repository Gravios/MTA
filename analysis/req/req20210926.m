% req20210926
%    Tags: network state interneurons immobility 
%    Status: Active
%    Type: Analysis
%    Author: Justin Graboski
%    Final_Forms: NA
%    Project: General
%
%
% Conclusions
%  

% Interneuron pca component -> foward egocentric phase precession
% rad/Lm theta power proxy  -> foward egocentric phase precession
% intershank proxy power difference -> foward egocentric phase precession
% intershank proxy power difference -> velocity
% intershank proxy power difference -> rad/lm power ratio



% LOAD Session data --------------------------------------------------------------------------------

configure_default_args();

MjgER2016_load_data();

% LOAD behavioral data
% $$$ [pfd ,tags ,eigVec, eigVar, eigScore, validDims,...
% $$$  unitSubsets, unitIntersection, zrmMean, zrmStd] = req20180123_ver5(Trials);
% $$$ numComp = size(eigVec{1},2);
% $$$ pfindex = 1;

% LOAD behavioral scores
% $$$ MjgER2016_load_bhv_erpPCA_scores();

trialId = 20;

% LOAD units  
Trial = Trials{trialId};    % jg05-20120312.cof.all
unitsSubset = Trial.spk.get_unit_set(Trial,'pyramidal');
unitsPyr    = Trial.spk.get_unit_set(Trial,'pyramidal');
unitsInt =    Trial.spk.get_unit_set(Trial,'interneurons');

% SORT pyramidal cells by spatial bins
pft = pfs_2d_theta(Trial,unitsPyr);
[mxr,mxp]= pft.get_max_rate(unitsPyr);
mxd = sqrt(sum(mxp.^2,2));
mxn = mxp./mxd;
mxa = atan2(mxn(:,2),mxn(:,1));
mxg = nan(size(mxd));
mxg(mxd<=175) = 1;
abins = linspace(-pi,pi,9);
for a = 1:numel(abins)-1
    mxg(mxd>175& mxa>abins(a) & mxa<=abins(a+1)) = a+1;
end    
[mxv,mxi] = sort(mxg,'descend');
unitsPyr = unitsPyr(mxi);
ucolors = hsv(8);
unitsPyrColor = zeros([size(mxg,1),3]);
unitsPyrColor(mxv>1,:) = ucolors(mxv(mxv>1)-1,:);
unitsIntColor = repmat([0,0,0],[numel(unitsInt),1]);

Trial.load('nq');


sampleRate = 30; % Hertz
%sampleRate = 250; % Hertz
%spikeWindow = 0.025; % Seconds
%halfSpkWindow = 0.012; % Seconds
% LOAD behavioral variables
xyz = preproc_xyz(Trial,'trb',sampleRate);
fxyz = filter(xyz.copy(),'ButFilter',3,14,'low');
vxy = vel(filter(xyz.copy(),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2]);
stc = Trial.stc.copy();
ang = create(MTADang,Trial,fxyz);
pch = fet_HB_pitchB(Trial,sampleRate);

bfs      = cf(@(t,u)  compute_bhv_ratemaps(t,u),          Trials, units);

% LOAD bfs erpPCA
[eigVecs, eigScrs, eigVars, unitSubset, validDims, zrmMean, zrmStd] = ...
                    compute_bhv_ratemaps_erpPCA(bfs, units, [], [], false);
numComp = size(eigVecs,2);

% LINEAR shank
Trial.lfp.filename = [Trial.name,'.lfp'];
lfp = Trial.load('lfp',65:96);
% BUZ64 shank
Trial.lfp.filename = [Trial.name,'.lfp'];
clfp = Trial.load('lfp',57:64);

% FILTER lfp <- 30Hz lowpass
fclfp = filter(copy(clfp),'ButFilter',4,30,'low');


% PLOT example lfp traces of proxy and different layers
figure,
hold('on');
hold('on');
shift = -[1:8];
for c = 1:8
plot([1:size(fclfp,1)]./fclfp.sampleRate,unity(clfp(:,c))-5+shift(c))
end
plot([1:size(fclfp,1)]./fclfp.sampleRate,unity(lfp(:,4))+5,'k')
plot([1:size(fclfp,1)]./fclfp.sampleRate,unity(diff(fclfp(:,[1,8]),1,2)),'b')
plot([1:size(fclfp,1)]./fclfp.sampleRate,unity(lfp(:,18)),'m');
legend({'c1','c2','c3','c4','c5','c6','c7','c8','pyr','env','lm'});


% COMPUTE spectrum <- BUZ64 shank, between 1 and 35 Hz
elfp = copy(fclfp);
elfp.resample(250);
specArgsTheta = struct('nFFT',2^8,...
                  'Fs',  elfp.sampleRate,...
                  'WinLength',2^7,...
                  'nOverlap',2^7*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[1,35]);
[eys,efs,ets] = fet_spec(Trial,elfp,[],false,[],specArgsTheta);

%%%<<< UTILITY - in progress
% FIND theta power center on shank
etpow = eys.copy();
etpow.data = 1./sq(log10(mean(eys(:,efs>5&efs<12,:),2)));
figure();
imagesc(ets,1:8,bsxfun(@rdivide,etpow.data',sum(etpow.data')));
caxis([0.12,0.13]);
colormap('jet');
Lines(nonzeros(get([stc{'t',1}],'data')),[],'m');
%%%>>>

elfp = copy(clfp);
elfp.data = [unity(unity(flfp(:,1))-unity(diff(fclfp(:,[1,8]),1,2))),...
             unity(diff(fclfp(:,[1,8]),1,2)),...
             unity(unity(diff(fclfp(:,[3,8]),1,2))-unity(diff(fclfp(:,[1,6]),1,2)))];
specArgsTheta = struct('nFFT',2^11,...
                  'Fs',  lfp.sampleRate,...
                  'WinLength',2^10,...
                  'nOverlap',2^10*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[1,35]);
   
[cys,cfs,cts] = fet_spec(Trial,elfp,[],false,[],specArgsTheta);

%%%<<< FIGURE - example of RAD/LM theta proxy
states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit'};
stateColors = 'krggbbmy';

[hfig,fig,fax,sax] = set_figure_layout(figure(323232),'A4','landscape','centimeters',25,1.75,0.3,0.3);
% CREATE subplot axes
sax = place_subplot(hfig,sax,[1, 0, 1, 0],fig);
    hold(sax(end),'on');
    imagesc(lts,lfs,log10(lys(:,:,4))');     colormap('jet');
    axis('xy');                              caxis([3.25,6.75]);
    axis(sax(end),'tight');
    title('CA1 Spectral Time Series');
    ylabel({'CA1:PYR','Hz'});
sax = place_subplot(hfig,sax,[2, 0, 1, 0],fig);
    hold(sax(end),'on');
    imagesc(lts,lfs,log10(lys(:,:,13))');    colormap('jet');
    axis('xy');                              caxis([3.25,6.75]);
    axis(sax(end),'tight');    
    sax(end).XTickLabel = {};
    ylabel({'CA1:RAD','Hz'});
sax = place_subplot(hfig,sax,[3, 0, 1, 0],fig);
    hold(sax(end),'on');
    imagesc(lts,lfs,log10(lys(:,:,21))');    colormap('jet');
    axis('xy');                              caxis([3.25,6.75]);
    axis(sax(end),'tight');    
    ylabel({'CA1:LM','Hz'});
sax = place_subplot(hfig,sax,[4, 0, 1, 0],fig);
    hold(sax(end),'on');
    imagesc(cts,cfs,log10(cys(:,:,2))');     colormap('jet');
    axis('xy');                              caxis([-0.5,2.5]);
    axis(sax(end),'tight');    
    ylabel({'CA1:PROXY','Hz'});
sax = place_subplot(hfig,sax,[5, 0, 1, 0],fig);
    hold(sax(end),'on');
    plotSTC(Trial.stc,1,'text',states,stateColors);
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
    ylim([1,9]);
    xlabel('Time (s)');

af(@(a) set(a,'XTickLabel',{}), sax(1:end-1));
linkaxes(sax,'x');

sax = place_subplot(hfig,sax,[6, -5, 1, 0],fig,5/fig.subplot.width,5/fig.subplot.height);
    hold(sax(end),'on');
    imagesc(crossChanThetaPowCorr');       colormap('jet');
    axis(sax(end),'ij');                   caxis([0,1]);
    axis(sax(end),'tight');
    sax(end).XTick = [3,7,12,18,23,28];
    sax(end).YTick = [3,7,12,18,23,28];    
    sax(end).XTickLabel = {'ori','pyr','rad','lm','dg','hil'};
    sax(end).YTickLabel = {'ori','pyr','rad','lm','dg','hil'};
    cax = colorbar(sax(end));
    cax.Units = 'centimeters';
    cax.Position(1) = sum(sax(end).Position([1,3]))+0.1;
    box(sax(end),'on');
sax = place_subplot(hfig,sax,[6, 0.15, 1, 0],fig,5/fig.subplot.width,0.25/fig.subplot.height);
    hold(sax(end),'on');
    imagesc(proxChanThetaPowCorr);         colormap('jet');
    axis(sax(end),'ij');                   caxis([0,1]);
    axis(sax(end),'tight');    
    ylabel('proxy','Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle');
    sax(end).XTickLabel = {};
    sax(end).YTickLabel = {};
    box(sax(end),'on');
    title(sax(end),{'Cross Channel','Theta Power Correlation'});

%%%>>>

%%%<<< FIGURE correlation of csd and rad/lm theta power as function of theta phase
pyrCSD = copy(lfp);
pyrCSD.data = unity(sum(fclfp(:,[1,8]))-2*fclfp(:,4));
csdLim = [-3,3];
csdBin = linspace([csdLim,36]);
csdCtr = mean([csdBin(1:end-1);csdBin(2:end)]);
csdInd = discretize(pyrCSD.data,csdBin);

plfp = copy(lfp);
plfp.data = diff(fclfp(:,[1,8]),1,2);
phzProx = plfp.phase([4,12]);
phzLim = [-pi,pi];
phzBin = linspace([phzLim,37]);
phzCtr = mean([phzBin(1:end-1);phzBin(2:end)]);
phzInd = discretize(phzProx.data,phzBin);

tper = [Trial.stc{'t-s-m'}];
tper.cast('TimeSeries');
tper.resample(htpow);
% $$$ wper = [Trial.stc{'w&t'}];
% $$$ wper.cast('TimeSeries',xyz);
% $$$ pper = [Trial.stc{'p&t'}];
% $$$ pper.cast('TimeSeries',xyz);
% $$$ sper = [Trial.stc{'s-t'}];
% $$$ sper.cast('TimeSeries',xyz);

figure,
nind = logical(tper.data) & nniz(phzInd) & nniz(pyrCSD);
out = accumarray([phzInd(nind)],pyrCSD(nind),[numel(phzCtr),1],@mean);
plot(out)

htpow =copy(tpow);
htpow.resample(250);

tper = [Trial.stc{'t-s-m'}];
tper.cast('TimeSeries');
tper.resample(htpow);

xmin = LocalMinima(abs(circ_dist(phzProx.data,phzCtr(18))),70,0.1);
vq = interp1(xmin./pyrCSD.sampleRate,pyrCSD.data(xmin),[1:size(htpow)]./htpow.sampleRate)';
nind = logical(tper.data) & nniz(htpow(:,14)) & nniz(vq);
figure,plot(htpow(nind,14),vq(nind)','.')

for p = 1:numel(phzCtr)
    xmin = LocalMinima(abs(circ_dist(phzProx.data,phzCtr(p))),60,pi/2);
    vq = interp1(xmin./pyrCSD.sampleRate,pyrCSD.data(xmin),[1:size(htpow)]./htpow.sampleRate)';
    nind = logical(tper.data) & nniz(htpow(:,14)) & nniz(vq);
for c = 1:32
ccc(c,p) = corr(htpow(nind,c)./htpow(nind,11),vq(nind));
end
end
figure,imagesc(phzCtr,1:32,ccc)
colormap('jet');
%%%>>>

% Correlation of lm rad ratio with proxy
% SCATTER PLOT - proxy vs layers
ind = [stc{'w+p&t'}];
figure();
subplot(131);
plot(mean(log10(lys(ind,lfs>6&lfs<12,7)),2),mean(log10(cys(ind,cfs>6&cfs<12,2)),2),'.');
ylabel('proxy theta power');
xlabel('CA1:pyr theta power');
subplot(132);
plot(mean(log10(lys(ind,lfs>6&lfs<12,12)),2),mean(log10(cys(ind,cfs>6&cfs<12,2)),2),'.');
ylabel('proxy theta power');
xlabel('CA1:rad theta power');
subplot(133);
plot(mean(log10(lys(ind,lfs>6&lfs<12,21)),2),mean(log10(cys(ind,cfs>6&cfs<12,2)),2),'.');
ylabel('proxy theta power');
xlabel('CA1:lm theta power');



% COMPUTE Inter channel theta power correlation
for c1 = 1:32
    for c2 = 1:32
    lp = mean(log10(lys(ind,lfs>6&lfs<12,c1)),2);
    tp = mean(log10(lys(ind,lfs>6&lfs<12,c2)),2);
    nind = nniz(lp)&nniz(tp);
    crossChanThetaPowCorr(c1,c2) = corr(lp(nind),tp(nind));
    end
end

proxChanThetaPowCorr = [];
tp = mean(log10(cys(ind,cfs>6&cfs<12,2)),2);%+mean(log10(cys(ind,cfs>6&cfs<12,2)),2);
for c = 1:32
    lp = mean(log10(lys(ind,lfs>6&lfs<12,c)),2);
    nind = nniz(lp)&nniz(tp);
    proxChanThetaPowCorr(c) = corr(lp(nind),tp(nind));
end

hfig = figure();
sax = gobjects([1,0]);
sax(end+1) = axes(hfig,'Position',[0.15,0.15,0.6,0.7]);
imagesc(ccc');
colormap('jet');
caxis([0,1]);
sax(end).XTickLabel = {};
sax(end).YTickLabel = {};
xlabel(sax(end),'ori          pyr             rad             lm           dg                     ');
ylabel(sax(end),'            dg           lm             rad             pyr          ori');
sax(end+1) = axes(hfig,'Position',[0.15,0.875,0.6,0.025]);
imagesc(rcc);
colormap('jet');
caxis([0,1]);
sax(end).XTickLabel = {};
sax(end).YTickLabel = {};
cax = colorbar(sax(1));
cax.Position(1) = sum(sax(2).Position([1,3]))+0.005;
ylabel(cax,'Theta Power Correlation');



% CROSS-SPECTRA of buz64 shank top and bottom channels
elfp = copy(clfp);
elfp.data = [fclfp(:,[1,8])];
specArgsTheta = struct('nFFT',2^10,...
                  'Fs',  lfp.sampleRate,...
                  'WinLength',2^9,...
                  'nOverlap',2^9*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[1,35]);
[ccys,ccfs,ccts,ccphi] = fet_spec(Trial,elfp,[],false,[],specArgsTheta,[],true);

figure,
subplot(411);
imagesc(ccts,ccfs,ccphi(:,:,1,2)');axis('xy');
colormap(gca(),'hsv');
subplot(412);
imagesc(cts,cfs,log10(lys(:,:,13))');axis('xy');
colormap(gca(),'jet');
subplot(413);
imagesc(cts,cfs,log10(lys(:,:,21))');axis('xy');
colormap(gca(),'jet');
subplot(414);
hold('on');
plotSTC(Trial.stc,1,'text',states,stateColors);
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
ylim([1,9]);
linkaxes(findobj(gcf(),'Type','Axes'),'x');


proxPhzDiff = copy(ccys);
proxPhzDiff.data = ccphi(:,:,1,2);
proxPhzDiff.resample(cys);


% CORRELATION between power and phase shift
bcorr = [];
breg = proxPhzDiff(ind,7);
for c = 1:32
    lrbratio = mean(log10(lys(ind,lfs>6&lfs<12,c)),2);
    nind = nniz(breg)&nniz(lrbratio);
    bcorr(c) = corr(breg(nind),lrbratio(nind));        
end
%figure,imagesc(1:32,ccfs(2:2:numel(ccfs)),bpcorr),axis('xy');
figure,plot(1:32,bcorr)


refChan = 7; 
bpcorr = [];
breg = proxPhzDiff(ind,7);
for c = 1:32
    %if c == refChan, continue;end;
    %for f = 2:2:numel(ccfs)
    for c2 = c+1:32
        %lrbratio = log2(mean(log10(lys(ind,lfs>6&lfs<12,c)),2)./mean(log10(lys(ind,lfs>6&lfs<12,c2)),2));
        lrbratio = log2(mean(log10(lys(ind,lfs>6&lfs<12,c)),2)./mean(log10(lys(ind,lfs>6&lfs<12,c2)),2));        

        nind = nniz(breg)&nniz(lrbratio);
        %bpcorr(f/2,c) = corr(breg(nind),lrbratio(nind));
        bpcorr(c,c2) = corr(breg(nind),lrbratio(nind));        
    end
end
%figure,imagesc(1:32,ccfs(2:2:numel(ccfs)),bpcorr),axis('xy');
figure,imagesc(1:32,1:32,bpcorr)

proxPhzDiff.resample(cys);
proxPhzDiff.resample(vxy);




%%%<<< UTILITY - segmentation of theta state and frequency estimation
rlcys = copy(cys);
rlcys.data = RectFilter(RectFilter(RectFilter(RectFilter(RectFilter(RectFilter(RectFilter(RectFilter(RectFilter(RectFilter(log10(cys(:,:,2)),3,1)',3,1)',3,1)',3,1)',3,1)',3,1)',3,1)',3,1)',3,1)',3,1)';

[mfv,mfi] = max(bsxfun(@rdivide,rlcys.data(:,cfs<12),sum(rlcys(:,cfs<12),2)),[],2);

figure,
subplot(411);
    hold('on');
    imagesc(lts,lfs,log10(lys(:,:,19))');
    axis('xy');
    colormap('jet');
    plot(cts,cfs(mfi));;
    caxis([4,7])
subplot(412);
    % $$$ imagesc(cts,cfs,rlcys.data');
    % $$$ axis('xy');
    % $$$ colormap('jet');
    % $$$ caxis([-1,2.5])
    %plot([1:size(dufrInt,1)]./dufrInt.sampleRate,diff(out(:,2:3),1,2));
    plot([1:size(dufrInt,1)]./dufrInt.sampleRate,(-RectFilter(out,3,5)./200).^(-1/2))
    %plot([1:size(dufrInt,1)]./dufrInt.sampleRate,RectFilter(out,21,3));
    legend({'1','2','3'})
subplot(413);
    hold('on');
    imagesc(cts,cfs,bsxfun(@rdivide,rlcys.data,sum(rlcys(:,cfs<12),2))');;
    plot(cts,cfs(mfi));;
    Lines([],5,'k');
    caxis([0.0,0.07])
    axis('xy');
    %caxis([-1,2]);
    colormap('jet');
subplot(414);
    hold('on');
    plotSTC(Trial.stc,1,'text',states,stateColors);
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
    %plot(lts,(mean(log10(lys(:,lfs<15,7)),2)-2)*5-10);
    plot(cts,RectFilter((mean(log10(cys(:,cfs<15,2)),2)-2).*3+5,5,3));
    plot(cts,RectFilter((mean(log10(cys(:,cfs<5,2)),2)-2).*3+5,5,3));
    Lines([],4,'k');
    % $$$ plot(cts,mean(log10(cys(:,cfs>6&cfs<12,2)),2)./mean(log10(cys(:,cfs<3,2)),2)+2);
    % $$$ plot(cts,mean(log10(cys(:,cfs>6&cfs<12,2)),2)./mean(log10(cys(:,cfs<3|(cfs>12&cfs<15),2)),2)+2);
% $$$ plot(cts,mean(log10(cys(:,cfs>6&cfs<12,2))+2,2));
ylim([1,9]);
linkaxes(findobj(gcf(),'Type','Axes'),'x');


figure,
hold('on');
%plot(out./200)
plot((-RectFilter(out,3,5)./200).^(-1/2))
%%%>>>


tper = [Trial.stc{'t-s-m'}];
tper.cast('TimeSeries',xyz);

wper = [Trial.stc{'w&t'}];
wper.cast('TimeSeries',xyz);

pper = [Trial.stc{'p&t'}];
pper.cast('TimeSeries',xyz);

sper = [Trial.stc{'s-t'}];
sper.cast('TimeSeries',xyz);

figure,
subplot(131);
hist2((-RectFilter(out(~logical(tper.data),[2,3]),3,5)./200).^(-1/1.5),100,linspace([0,10,100]))
subplot(132);
hist2((-RectFilter(out(logical(tper.data),[2,3]),3,5)./200).^(-1/1.5),100,linspace([0,10,100]))
subplot(133);
hist2((-RectFilter(out(logical(tper.data),[2,1]),3,5)./200).^(-1/1.5),100,linspace([0,10,100]))

figure
hist2((-RectFilter(out(logical(tper.data),[2,3]),3,5)./200).^(-1/1.5),linspace([2,8,30]),linspace([0,5,30]))

intRunLim = [2,9];
intRunBin = linspace([intRunLim,10]);
intRunCtr = mean([intRunBin(1:end-1);intRunBin(2:end)]);
intRunInd = discretize((-RectFilter(out(:,[2]),3,1)./200).^(-1/1.5),intRunBin);
%intRunInd = discretize((-RectFilter(out(:,[2]),3,5)./200).^(-1/1.5),intRunBin);
intSpwLim = [1,9];
intSpwBin = linspace([intSpwLim,10]);
intSpwCtr = mean([intSpwBin(1:end-1);intSpwBin(2:end)]);
intSpwInd = discretize((-RectFilter(out(:,[1]),3,1)./200).^(-1/1.5),intSpwBin);
%intSpwInd = discretize((-RectFilter(out(:,[1]),3,5)./200).^(-1/1.5),intSpwBin);

lvxyLim = [-2,2];
lvxyBin = linspace([lvxyLim,10]);
lvxyCtr = mean([lvxyBin(1:end-1);lvxyBin(2:end)]);
lvxyInd = discretize(lvxy(:,2),lvxyBin);


efreq = copy(rlcys);
efreq.data = cfs(mfi);
efreq.resample(xyz);
efreqLim =[5,11];
efreqBin = linspace([efreqLim,12]);
efreqInd = discretize(efreq.data,efreqBin);

dcind = false([size(xyz,1),1]);
%dcind(dc.ind) = true;
dcind(dc.ind) = dc.ucnt>3&sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))<350;

hba = nan([size(xyz,1),1]);
hba(dc.ind) = dc.hbang.data;

%nind = nniz(intSpwInd)&nniz(intRunInd)&nniz(erf)&logical(tper.data)&dcind&hba>0.1;
%nind = nniz(intSpwInd)&nniz(intRunInd)&nniz(erf)&logical(sper.data)&dcind;
nind = nniz(intSpwInd)&nniz(intRunInd)&nniz(erf)&logical(wper.data|pper.data)&dcind;
nind = nniz(lvxyInd)&nniz(intSpwInd)&nniz(intRunInd)&nniz(erf)&logical(pper.data)&dcind;
mout = accumarray([intRunInd(nind),intSpwInd(nind)],erf(nind),[numel(intRunBin)-1,numel(intSpwBin)-1],@mean);
sout = accumarray([intRunInd(nind),intSpwInd(nind)],erf(nind),[numel(intRunBin)-1,numel(intSpwBin)-1],@std);
figure,
subplot(131);
imagesc(intRunCtr,intSpwCtr,mout');
colormap('jet');
axis('xy');
%caxis([-20,20])
%caxis([-100,100])
caxis([-20,80])
%caxis([-20,200])
%caxis([0,200])
subplot(132);
imagesc(intRunCtr,intSpwCtr,sout');
colormap('jet');
axis('xy');
caxis([0,200])
subplot(133);
hist2((-RectFilter(out(nind,[2,1]),3,1)./200).^(-1/1.5),intRunBin,intSpwBin);


rtpow = cys.copy();
%rtpow.data = mean(log10(cys(:,cfs>5&cfs<12,2)),2);
rtpow.data = nunity(mean(log10(cys(:,cfs>5&cfs<12,2)),2)) - nunity(mean(log10(lys(:,cfs>5&cfs<12,7)),2));
rtpow.resample(xyz);
%rtpowLim = [0.5,3];
%rtpowLim = [0.15,0.5];
rtpowLim = [-2,2];
rtpowBin = linspace([rtpowLim,10]);
rtpowCtr = mean([rtpowBin(1:end-1);rtpowBin(2:end)]);
rtpowInd = discretize(rtpow.data,rtpowBin);


tpowLim = [3.5,6.5];
tpowBin = linspace([tpowLim,20]);
tpowCtr = mean([tpowBin(1:end-1);tpowBin(2:end)]);
tpowInd = discretize(tpow(:,7),tpowBin);



%nind = nniz(intSpwInd)&nniz(intRunInd)&nniz(erf)&logical(tper.data)&dcind&hba>0.1;
%nind = nniz(intSpwInd)&nniz(intRunInd)&nniz(erf)&logical(sper.data)&dcind;
%nind = nniz(intSpwInd)&nniz(intRunInd)&nniz(erf)&logical(wper.data|pper.data)&dcind;
nind = nniz(rtpowInd)&nniz(intSpwInd)&nniz(intRunInd)&nniz(erf)&logical(pper.data)&dcind;
mout = accumarray([rtpowInd(nind),intSpwInd(nind)],erf(nind),[numel(rtpowBin)-1,numel(intSpwBin)-1],@mean);
sout = accumarray([rtpowInd(nind),intSpwInd(nind)],erf(nind),[numel(rtpowBin)-1,numel(intSpwBin)-1],@std);
figure,
subplot(131);
imagesc(rtpowCtr,intSpwCtr,mout');
colormap('jet');
axis('xy');
%caxis([-20,20])
%caxis([-100,100])
caxis([-30,100])
%caxis([-20,200])
%caxis([0,200])
subplot(132);
imagesc(rtpowCtr,intSpwCtr,sout');
colormap('jet');
axis('xy');
caxis([0,100])
%caxis([0,200])
subplot(133);
hist2([rtpow(nind,1),(-RectFilter(out(nind,[1]),3,1)./200).^(-1/1.5)],rtpowBin,intSpwBin);


%nind = nniz(intSpwInd)&nniz(intRunInd)&nniz(erf)&logical(tper.data)&dcind&hba>0.1;
%nind = nniz(intSpwInd)&nniz(intRunInd)&nniz(erf)&logical(sper.data)&dcind;
%nind = nniz(intSpwInd)&nniz(intRunInd)&nniz(erf)&logical(wper.data|pper.data)&dcind;
nind = nniz(tpowInd)&nniz(intSpwInd)&nniz(intRunInd)&nniz(erf)&logical(wper.data|pper.data)&dcind;
mout = accumarray([tpowInd(nind),intSpwInd(nind)],erf(nind),[numel(tpowBin)-1,numel(intSpwBin)-1],@mean);
sout = accumarray([tpowInd(nind),intSpwInd(nind)],erf(nind),[numel(tpowBin)-1,numel(intSpwBin)-1],@std);
figure,
subplot(131);
imagesc(tpowCtr,intSpwCtr,mout');
colormap('jet');
axis('xy');
%caxis([-20,20])
%caxis([-100,100])
caxis([-20,80])
%caxis([-20,200])
%caxis([0,200])
subplot(132);
imagesc(tpowCtr,intSpwCtr,sout');
colormap('jet');
axis('xy');
caxis([0,200])
subplot(133);
hist2([tpow(nind,1),(-RectFilter(out(nind,[1]),3,1)./200).^(-1/1.5)],tpowBin,intSpwBin);


%nind = nniz(intSpwInd)&nniz(intRunInd)&nniz(erf)&logical(tper.data)&dcind&hba>0.1;
%nind = nniz(intSpwInd)&nniz(intRunInd)&nniz(erf)&logical(sper.data)&dcind;
%nind = nniz(intSpwInd)&nniz(intRunInd)&nniz(erf)&logical(wper.data|pper.data)&dcind;
nind = nniz(tpowInd)&nniz(lvxyInd)&nniz(intSpwInd)&nniz(intRunInd)&nniz(erf)&logical(wper.data|pper.data)&dcind;
mout = accumarray([lvxyInd(nind),tpowInd(nind)],erf(nind),[numel(lvxyBin)-1,numel(tpowBin)-1],@mean);
sout = accumarray([lvxyInd(nind),tpowInd(nind)],erf(nind),[numel(lvxyBin)-1,numel(tpowBin)-1],@std);
figure,
subplot(131);
imagesc(lvxyCtr,tpowCtr,mout');
colormap('jet');
axis('xy');
%caxis([-20,20])
%caxis([-100,100])
caxis([-20,120])
%caxis([-20,200])
%caxis([0,200])
subplot(132);
imagesc(lvxyCtr,tpowCtr,sout');
colormap('jet');
axis('xy');
caxis([0,200])
subplot(133);
hist2([lvxy(nind,2),tpow(nind,7)],lvxyBin,tpowBin);

%nind = nniz(intSpwInd)&nniz(intRunInd)&nniz(erf)&logical(tper.data)&dcind&hba>0.1;
%nind = nniz(intSpwInd)&nniz(intRunInd)&nniz(erf)&logical(sper.data)&dcind;
%nind = nniz(intSpwInd)&nniz(intRunInd)&nniz(erf)&logical(wper.data|pper.data)&dcind;
nind = nniz(rtpowInd)&nniz(lvxyInd)&nniz(intSpwInd)&nniz(intRunInd)&nniz(erf)&logical(pper.data)&dcind;
mout = accumarray([lvxyInd(nind),rtpowInd(nind)],erf(nind),[numel(lvxyBin)-1,numel(rtpowBin)-1],@mean);
sout = accumarray([lvxyInd(nind),rtpowInd(nind)],erf(nind),[numel(lvxyBin)-1,numel(rtpowBin)-1],@std);
figure,
subplot(131);
imagesc(lvxyCtr,rtpowCtr,mout');
colormap('jet');
axis('xy');
%caxis([-20,20])
%caxis([-100,100])
caxis([-20,80])
%caxis([-20,200])
%caxis([0,200])
subplot(132);
imagesc(lvxyCtr,rtpowCtr,sout');
colormap('jet');
axis('xy');
caxis([0,200])
subplot(133);
hist2([lvxy(nind,2),rtpow(nind)],lvxyBin,rtpowBin);




mout = accumarray([lvxyInd(nind)],erf(nind,1),[numel(lvxyBin)-1,1],@mean);
sout = accumarray([lvxyInd(nind)],erf(nind,1),[numel(lvxyBin)-1,1],@std);
figure,
hold('on');
plot(lvxy(nind,2),erf(nind),'.');
plot(lvxyCtr,mout);
plot(lvxyCtr,mout+sout,'r');
plot(lvxyCtr,mout-sout,'r');



mout = accumarray([intSpwInd(nind)],lvxy(nind,1),[numel(intSpwBin)-1,1],@mean);
sout = accumarray([intSpwInd(nind)],lvxy(nind,1),[numel(intSpwBin)-1,1],@std);
figure,
hold('on');
plot([(-RectFilter(out(nind,[1]),3,1)./200).^(-1/1.5)],lvxy(nind,1),'.');
plot(intSpwCtr,mout);
plot(intSpwCtr,mout+sout,'r');
plot(intSpwCtr,mout-sout,'r');

mout = accumarray([intSpwInd(nind)],lvxy(nind,2),[numel(intSpwBin)-1,1],@mean);
sout = accumarray([intSpwInd(nind)],lvxy(nind,2),[numel(intSpwBin)-1,1],@std);
figure,
hold('on');
plot([(-RectFilter(out(nind,[1]),3,1)./200).^(-1/1.5)],lvxy(nind,2),'.');
plot(intSpwCtr,mout);
plot(intSpwCtr,mout+sout);
plot(intSpwCtr,mout-sout);

mout = accumarray([intSpwInd(nind)],erf(nind),[numel(intSpwBin)-1,1],@mean);
sout = accumarray([intSpwInd(nind)],erf(nind),[numel(intSpwBin)-1,1],@std);
figure,
hold('on');
plot([(-RectFilter(out(nind,[1]),3,1)./200).^(-1/1.5)],erf(nind),'.')
plot(intSpwCtr,mout);
plot(intSpwCtr,mout+sout,'r');
plot(intSpwCtr,mout-sout,'r');
ylim([-200,200])

mout = accumarray([intSpwInd(nind)],erl(nind),[numel(intSpwBin)-1,1],@mean);
sout = accumarray([intSpwInd(nind)],erl(nind),[numel(intSpwBin)-1,1],@std);
figure,
hold('on');
plot([(-RectFilter(out(nind,[1]),3,1)./200).^(-1/1.5)],erl(nind),'.')
plot(intSpwCtr,mout);
plot(intSpwCtr,mout+sout,'r');
plot(intSpwCtr,mout-sout,'r');
ylim([-200,200])

[R,P] = corr([[intSpwInd(nind)],erf(nind)]);
[R,P] = corr([[intRunInd(nind)],erf(nind)]);

nind = nniz(efreqInd)&nniz(intRunInd)&nniz(erf)&logical(pper.data)&dcind;
mout = accumarray([efreqInd(nind),intRunInd(nind)],erf(nind),[numel(efreqBin)-1,numel(intRunBin)-1],@mean);
sout = accumarray([efreqInd(nind),intRunInd(nind)],erf(nind),[numel(efreqBin)-1,numel(intRunBin)-1],@std);

figure,
subplot(131);
imagesc(mout');
colormap('jet');
axis('xy');
%caxis([-20,20])
%caxis([-100,100])
%caxis([-20,80])
%caxis([0,200])
subplot(132);
imagesc(sout');
colormap('jet');
axis('xy');
caxis([0,150])
subplot(133);
%hist2((-RectFilter(out(logical(tper.data),[2,1]),3,1)./200).^(-1/1.5),intRunBin,intRunBin);
hist2([efreq(nind),(-RectFilter(out(nind,[1]),3,1)./200).^(-1/1.5)],efreqBin,intRunBin);

figure,
subplot(121);
imagesc(rlcys(48242,1:2:24)'-rlcys(48242,1:2:24));
caxis([-1,1]);
colormap('jet');
subplot(122);
imagesc(rlcys(48046,1:2:24)'-rlcys(48046,1:2:24));
caxis([-1,1]);
colormap('jet');

sum(sum((rlcys(48046,1:2:24)'-rlcys(48046,1:2:24)).*(rlcys(48242,1:2:24)'-rlcys(48242,1:2:24))))

sum(sum((rlcys(48164,1:2:24)'-rlcys(48164,1:2:24)).*(rlcys(48242,1:2:24)'-rlcys(48242,1:2:24))))

mg = nan([size(rlcys,1),1]);
for tt = 48000:49200
mg(tt) = log10(sqrt(sum(sum(diff(diff((cys(tt,1:2:24)'-cys(tt,1:2:24)).*(cys(48242,1:2:24)'-cys(48242,1:2:24)),1,2),1,1).^2))));
end
plot(cts,mg,'m')


figure,
hold('on');
plot(cts,mean(log10(cys(:,cfs>6&cfs<12,2)),2)./mean(log10(cys(:,cfs<20,2)),2));
plot(cts,mean(log10(cys(:,cfs>6&cfs<12,2)),2)./mean(log10(cys(:,cfs<5,2)),2));
plot(cts,mean(log10(cys(:,cfs>6&cfs<12,2)),2));


% IS the phase difference between the top and bottom channels indicative of rad/lm tpow ratio
proxPhzDiff = copy(ccys);
proxPhzDiff.data = circ_mean(ccphi(:,ccfs>2&ccfs<10,1,2),[],2);
proxPhzDiff.resample(cys);

ind = [stc{'w&t'}];
breg = proxPhzDiff(ind);
lrbratio = mean(log10(lys(ind,lfs>6&lfs<12,14)),2)./mean(log10(lys(ind,lfs>6&lfs<12,21)),2);
figure,plot(breg,lrbratio,'.')
nind = nniz(breg)&nniz(lrbratio);
corr(breg(nind),lrbratio(nind))


lrbratio = mean(log10(lys(ind,lfs>6&lfs<12,14)),2)./mean(log10(lys(ind,lfs>6&lfs<12,21)),2);
breg = mean(log10(cys(ind,lfs>6&lfs<12,2)),2);
figure,plot(breg,lrbratio,'.')
nind = nniz(breg)&nniz(lrbratio);
corr(breg(nind),lrbratio(nind))


lrbratio = mean(log10(lys(ind,lfs>6&lfs<12,13)),2)./mean(log10(lys(ind,lfs>6&lfs<12,21)),2);
breg = mean(log10(cys(ind,lfs>6&lfs<12,2)),2);%-polyval(pc(c,:),bratio);
oreg = sqrt(mean(log10(cys(ind,lfs>6&lfs<12,3)),2));
figure,plot3(lrbratio,breg,oreg,'.');


lrbratio = mean(log10(lys(ind,lfs>6&lfs<12,13)),2);
breg = mean(log10(cys(ind,lfs>6&lfs<12,2)),2);%-polyval(pc(c,:),bratio);
figure,plot(breg,lrbratio,'.')
nind = nniz(breg)&nniz(lrbratio);
corr(breg(nind),lrbratio(nind))

figure,plot(bratio,lrbratio,'.')
figure,plot(bratio,breg,'.')
nind = nniz(bratio)&nniz(lrbratio);
corr(bratio(nind),lrbratio(nind))

% pyr power vs lrbratio
lrbratio = mean(log10(lys(ind,lfs>6&lfs<12,12)),2)./mean(log10(lys(ind,lfs>6&lfs<12,21)),2);
breg = mean(log10(lys(ind,lfs>6&lfs<12,21)),2);%-polyval(pc(c,:),bratio);
figure,plot(breg,lrbratio,'.');
nind = nniz(breg)&nniz(lrbratio);
corr(breg(nind),lrbratio(nind))

nind = nniz(breg)&nniz(bratio);
corr(breg(nind),bratio(nind))
figure,plot(breg,bratio,'.');
figure,plot(breg,mean(log10(lys(ind,lfs>6&lfs<12,14)),2)-mean(log10(lys(ind,lfs>6&lfs<12,21)),2),'.');
figure,plot(breg,mean(log10(lys(ind,lfs>6&lfs<12,21)),2),'.');


cphz = clfp.phase;

figure,
subplot(2,1,1);
hold('on');
plot([1:size(clfp,1)]./clfp.sampleRate,circ_dist(cphz(:,1),cphz(:,8)))
plot([1:size(clfp,1)]./clfp.sampleRate,unity(sum(fclfp(:,1:end-1)-fclfp(:,2:end),2)),'b')
subplot(2,1,2);
hold('on');
plotSTC(Trial.stc,1,'text',states,stateColors);
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
ylim([1,9]);
linkaxes(findobj(gcf(),'Type','Axes'),'x');

dcphz = copy(cphz);
dcphz.data = circ_dist(cphz(:,1),cphz(:,8));
figure
hold('on');
ind = [stc{'lloc'}];
histogram(dcphz(ind),linspace(-pi,pi,100),'EdgeColor','none');
ind = [stc{'hloc'}];
histogram(dcphz(ind),linspace(-pi,pi,100),'EdgeColor','none');


figure,
for c = 1:3,
subplot(4,1,c);
imagesc(cts,cfs,(bsxfun(@rdivide,log10(cys(:,:,c)),sum(log10(cys(:,:,c)),2)))');
axis('xy');
end
colormap('jet');
subplot(4,1,4);
hold('on');
plotSTC(Trial.stc,1,'text',states,stateColors);
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
ylim([1,9]);
linkaxes(findobj(gcf(),'Type','Axes'),'x');





phz = load_theta_phase(Trial,...
                       lfp.sampleRate,...
                       sessionList(trialId).thetaRefGeneral,...
                       phzCorrection(trialId));


spk = Trial.load('spk',Trial.lfp.sampleRate,'',[unitsPyr,unitsInt]);
pyr = Trial.load('spk',Trial.lfp.sampleRate,'',unitsPyr);
pyrt = Trial.load('spk',Trial.lfp.sampleRate,'theta-sit-groom',unitsPyr);
int = Trial.load('spk',Trial.lfp.sampleRate,'theta-sit-groom',unitsInt);
% $$$ non = Trial.load('spk',Trial.lfp.sampleRate,'',unitsNon);


% ORDER interneurons based on theta phase preference
intPhzPref = zeros([numel(unitsInt),1]);
for ii = 1:numel(unitsInt),
    intPhzPref(ii) = circ_mean(phz(int(unitsInt(ii))));
end
intPhzPref(intPhzPref<0) = intPhzPref(intPhzPref<0)+2*pi;
[mpv,mpi] = sort(intPhzPref,'descend');
unitsInt = unitsInt(mpi);





% Using Only pyramidal cells to determine theta phase preference of interneurons
figure();
for intInd= 1:numel(unitsInt),
ires = int(unitsInt(intInd));
[mccg,tbins] = CCG([pyrt.res;ires],[ones(size(pyrt.res));2*ones(size(ires))],5,2^7,lfp.sampleRate,[1,2]);

x = tbins;
y = RectFilter(mccg(:,1,2)-mean(mccg(:,1,2)),5)'-RectFilter(mccg(:,1,2)-mean(mccg(:,1,2)),11)';
yr = (max(y)-min(y));                               % Range of ?y?
yz = y-max(y)+(yr/2);
zx = x(yz .* circshift(yz,[0 1]) <= 0);     % Find zero-crossings
per = 2*mean(diff(zx));                     % Estimate period

dy = diff(y);
dzx = dy(yz .* circshift(yz,[0 1]) <= 0);
pinds = [(find(zx<0,1,'last')),(find(zx>0,1,'first'))];
ccgphz = (zx(pinds(find(dzx(pinds)>0,1,'first'))));

ym = mean(y);                               % Estimate offset
p = polyfit(tbins((((numel(y)-1)/2)+2):end),mean([fliplr(y(1:((numel(y)-1)/2)));y((((numel(y)-1)/2)+2):end)]),1);
[xp,xpi] = min(abs(zx));
xp = xp.*sign(zx(xpi));
b = [ per;  ccgphz; 150; 1];

fit = @(b,x) b(4).*exp(-0.5.*(x.^2./b(3).^2)).*  ...
              ((cos(2*pi*(x + b(2))./b(1) )));    % Function to fit
fcn = @(b) sum((fit(b,x) - y).^2);                              % Least-Squares cost function
s(:,intInd) = fminsearch(fcn, b,struct('MaxIter',10000)) ;                      % Minimise Least-Squares
xp = linspace(min(x),max(x));
subplot(4,4,intInd);
%figure
hold('on');
plot(x,s(4,intInd).*exp(-0.5.*(x.^2./s(3,intInd).^2)).*cos(2.*pi.*(x+s(2,intInd))./s(1,intInd)));
plot(x,y);
Lines([],0,'r');
Lines(0,[],'r');
title(num2str(unitsInt(intInd)));
end

mpvp = 2*pi*s(2,:)'./s(1,:)';
mpvp(mpvp<-pi) = mpvp(mpvp<-pi)+2*pi;
mpvp(mpvp<-pi/2) = mpvp(mpvp<-pi/2)+pi;
mpvp = -(mpvp+pi/2)+2*pi
% $$$ mpvp = mod(-mpvp-pi/2,2*pi);
% $$$ mpvp(mpvp<pi/2) = mpvp(mpvp<pi/2)+pi;
figure,
%plot([mpv,mod(-mpvp'-pi/2,2*pi)],'-+')
plot([mpv,mpvp],'-+');
Lines([],-pi,'r');
Lines([],-pi/2,'r');
Lines([],pi,'r');
Lines([],pi/2,'r');
Lines([],3*pi/2,'r');
Lines([],-3*pi/2,'r');




% Decoding stuff
ufr = Trial.load('ufr',xyz,[],unitsPyr,spikeWindow,'boxcar');
pfs = pfs_2d_states(Trial,unitsPyr,[],{'loc&theta'});
ratemaps = [];
for unit = unitsPyr
    ratemaps(:,end+1) = pfs{1}.data.rateMap(:,pfs{1}.data.clu==unit);
end
mask = create_tensor_mask(pfs{1}.adata.bins);
ratemaps(~mask(:),:) ==nan;

bhvRatemaps = [];
for unit = unitsPyr
    bhvRatemaps(:,end+1) = bfs{trialId}.data.rateMap(:,bfs{trialId}.data.clu==unit);
    bhvRatemaps(~validDims,end) = nan;
    bhvRatemaps(validDims&~nniz(bhvRatemaps),end) = 0.01;
end




% INTERNEURON theta phase preference rate ratio
ufrInt = Trial.load('ufr',lfp,spk,unitsInt,0.2,'gauss');
fufrInt = filter(copy(ufrInt),'ButFilter',4,1.2,'low');


lvxy = copy(vxy);
lvxy.data(lvxy.data<=0.0001) = 0.0001;
lvxy.data = log10(lvxy.data);

figure();
hist2(lvxy(nniz(lvxy),:),linspace(-3.5,1.8,20),linspace(-3.5,1.8,20));

velBins = linspace(-2.5,1.8,20);
velInds = discretize(lvxy(:,2),velBins);

stcm = stc2mat(Trial.stc,xyz,{'theta','rear','loc','pause','sit'});

figure();
hist2(lvxy(nniz(lvxy) & stcm(:,4)==4,:),linspace(-3.5,1.8,20),linspace(-3.5,1.8,20));


dufrInt = Trial.load('ufr',xyz,spk,unitsInt,0.12,'gauss');
dufrInt = resample(copy(ufrInt),xyz);

figure();
for vind = 1:20;
subplot(4,5,vind);
imagesc((bsxfun(@minus,dufrInt(velInds==vind & stcm(:,4)==4,:),median(dufrInt(velInds==vind & stcm(:,4)==4,:)))'*bsxfun(@minus,dufrInt(velInds==vind & stcm(:,4)==4,:),median(dufrInt(velInds==vind & stcm(:,4)==4,:))))./sum(velInds==vind & stcm(:,4)==4));
end

ind = (stcm(:,4)==4|stcm(:,5)==5) & stcm(:,1)~=1;
%ind = (stcm(:,3)==3|stcm(:,5)==5) & stcm(:,1)~=1;
suI = dufrInt(ind,:);
svI = lvxy(ind,:);
nClust = 3;
[idx,C,SUMD,D] = kmeans(suI,nClust,'Distance','correlation');


kCov = zeros(size(dufrInt,2),size(dufrInt,2),nClust);
for vind = 1:nClust,
    kCov(:,:,vind) = (bsxfun(@minus,suI(idx==vind,:),median(suI(idx==vind,:)))' ...
                       *bsxfun(@minus,suI(idx==vind,:),median(suI(idx==vind,:)))) ...
                      ./sum(idx==vind);
end

figure,
for vind = 1:nClust,
    subplot2(nClust,2,vind,1);
        imagesc(kCov(:,:,vind));
    subplot2(nClust,2,vind,2);
    histogram(svI(idx==vind,2),velBins);
end

out = [];
for i =  1:nClust
    out(:,i) = -.5*log(det(kCov(:,:,i)))...
        -.5*dot((bsxfun(@minus,dufrInt(:,:),median(suI(idx==vind,:)))/kCov(:,:,i))',...
                 bsxfun(@minus,dufrInt(:,:),median(suI(idx==vind,:)))')';
end


figure,
subplot(211);
%plot([1:size(dufrInt,1)]./dufrInt.sampleRate,diff(out(:,2:3),1,2));
plot([1:size(dufrInt,1)]./dufrInt.sampleRate,out);
legend({'1','2','3'})
subplot(212);
hold('on');
plotSTC(Trial.stc,1,'text',states,stateColors);
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
ylim([1,9]);
linkaxes(findobj(gcf(),'Type','Axes'),'x');


% NEXT hmmg
% mean isi of sliding window
res = spk(43);


mccg = [];
t = 5250000;
for t = 5250000:5260000
[mccg(t-5250000+1,:),tbins] = CCG(res,1,5,60,1250,1,'hz',[t-60,t+60]);
end

figure,imagesc([5250000:5260000]./1250,tbins,mccg');
figure,plot([5250000:5260000]./1250,sum(mccg,2,'omitnan'));
figure,plot([5250000:5260000]./1250,dufrInt(round((5250000:5260000)./1250*xyz.sampleRate),1));


figure();
subplot(311);
hold('on');
% $$$ plot([1:size(ufrInt,1)]/lfp.sampleRate,sum(ufrInt(:,mpv>2&mpv<3),2)./sum(mpv>2&mpv<3));
% $$$ plot([1:size(ufrInt,1)]/lfp.sampleRate,sum(ufrInt(:,mpv>3.1),2)./sum(mpv>3.1));
% $$$ plot([1:size(ufrInt,1)]/lfp.sampleRate,sum(fufrInt(:,mpv>2&mpv<3),2)./sum(mpv>2&mpv<3));
% $$$ plot([1:size(ufrInt,1)]/lfp.sampleRate,sum(fufrInt(:,mpv>3.1),2)./sum(mpv>3.1));
plot([1:size(dufrInt,1)]./dufrInt.sampleRate,out);
subplot(312);
imagesc(cts,cfs,log10(cys(:,:,2))');
axis('xy');
colormap('jet');
subplot(313);
hold('on');
plotSTC(Trial.stc,1,'text',states,stateColors);
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
ylim([1,9]);
linkaxes(findobj(gcf(),'Type','Axes'),'x');

figure,plot(log2(out(:,1)./out(:,2)),log2(out(:,1)./out(:,3)),'.');

ind = (stcm(:,4)==4|stcm(:,5)==5) & stcm(:,1)~=1;
ind = (stcm(:,5)==5) & stcm(:,1)~=1;
ind = (stcm(:,3)==3|stcm(:,4)==4) & stcm(:,1)==1;
figure,
scatter(log2(out(ind,1)./out(ind,2)),tpow(ind,7),10,erd(ind,1),'Filled');
colormap('jet');


intFet = log2(out(:,3)./out(:,1));

icrLim = [-3,4];
icrBin = linspace([icrLim,30]);
icrCtr = mean([icrBin(1:end-1);icrBin(2:end)]);
icrInd = discretize(intFet,icrBin);
%figure,plot(lwp(:),log2(out(:,2)./out(:,3)),'.');


nind = nniz([lwpInd,icrInd,erd]) & ind;
figure,
subplot(141);
hist2([lwp(nind),intFet(nind)],lwpBin,icrBin);
colorbar();
caxis([0,500])
subplot(142);
tpout = accumarray([lwpInd(nind),icrInd(nind)],lvxy(nind,2),[numel(lwpCtr),numel(icrCtr)],@mean);
imagesc(lwpCtr,icrCtr,tpout');  colorbar();        caxis([-2,0.5]);
axis('xy');                     colormap('jet');
subplot(143);
tpout = accumarray([lwpInd(nind),icrInd(nind)],erd(nind,1),[numel(lwpCtr),numel(icrCtr)],@mean);
imagesc(lwpCtr,icrCtr,tpout');  colorbar();
axis('xy');                      colormap('jet');
subplot(144);
tpout = accumarray([lwpInd(nind),icrInd(nind)],erf(nind,1),[numel(lwpCtr),numel(icrCtr)],@mean);
imagesc(lwpCtr,icrCtr,tpout');  colorbar();
axis('xy');                      colormap('jet');

ForAllSubplots('Lines([],0,''r'');Lines(0,[],''r'');')

% Comparison of spectra from depth recordings in the hippocampal CA1
% Preprocess - remove movement artifact via ica (weiwei)
% Alginment - Comparison of the spectra requires that anatomical locations are matched/mapped to the 
%             selected channels.
%    Method - CSD -> blocked and averaged over time to estimate the location of the primary theta 
%             sources and sinks.
%             - Function: csd
%    Method - ripple reversal point -> blocked and averaged over time to estimate the location of the
%             pyramidal layer.
%    

step = 1;
ch=[step+1:32-step];
mcsd = resample(lfp.copy(),60);
mcsd.data = mcsd(:,ch+step) - 2*mcsd(:,ch) + mcsd(:,ch-step);
mcsd.data = -mcsd.data/(step^2);


rlfp = resample(lfp.copy(),30);
bcsd = copy(rlfp);
[bcsd.data,nts,nch,dcsd] = CurSrcDns(bcsd.data,[],'c',[],[],60,2);

figure();
subplot(311);
    %imagesc(nts./960,1:size(bcsd,2),bcsd.data');
    imagesc([1:size(dcsd,1)]./30,4:29,dcsd');
    Lines([],3,'m');
    Lines([],9,'m');    
    Lines([],17,'m');
% $$$     Lines([],217/30*3,'m');
% $$$     Lines([],217/30*9,'m');    
% $$$     Lines([],217/30*17,'m');
    caxis([-2000,2000]);
subplot(312);    
    imagesc(lts,1:32,sq(mean(log10(lys(:,lfs>5&lfs<12,:)),2))');
    axis(sax(1),'ij');
    colormap('jet');
    caxis([4,6.7]);
subplot(313);
    hold('on');
    plotSTC(Trial.stc,1,'text',states,stateColors);
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
    ylim([1,9]);
    plot([1:size(errdist,1)]./30,errdist./50,'m')
linkaxes(findobj(gcf(),'Type','Axes'),'x');




%%%<<< LINEAR 32 lfp spectra -----------------------------------------------------------------------
specArgsTheta = struct('nFFT',2^11,...
                  'Fs',  lfp.sampleRate,...
                  'WinLength',2^10,...
                  'nOverlap',2^10*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[1,35]);
   
[lys,lfs,lts] = fet_spec(Trial,lfp,[],[],[],specArgsTheta);

% COMPUTE Inter channel theta power correlation
crossChanThetaPowCorr = nan([32,32]);
for c1 = 1:32
    for c2 = 1:32
    lp = mean(log10(lys(ind,lfs>6&lfs<12,c1)),2);
    tp = mean(log10(lys(ind,lfs>6&lfs<12,c2)),2);
    nind = nniz(lp)&nniz(tp);
    crossChanThetaPowCorr(c1,c2) = corr(lp(nind),tp(nind));
    end
end


% COMPUTE Inter channel theta power correlation
mavgCrossChanThetaPowCorr = nan([100,32,32]);
for tind = 51:100
for c1 = 1:32
    for c2 = 1:32
    lp = mean(log10(lys(tind-50:tind+50,lfs>6&lfs<12,c1)),2);
    tp = mean(log10(lys(tind-50:tind+50,lfs>6&lfs<12,c2)),2);
    nind = nniz(lp)&nniz(tp);
    mavgCrossChanThetaPowCorr(tind,c1,c2) = corr(lp(nind),tp(nind));
    end
end
end

figure,
subplot(141);
imagesc(sq(mavgCrossChanThetaPowCorr(51,:,:))),colorbar();
subplot(142);
imagesc(sq(mavgCrossChanThetaPowCorr(74,:,:))),colorbar();
subplot(143);
imagesc(sq(mavgCrossChanThetaPowCorr(25000,:,:))),colorbar();
subplot(144);
imagesc(sq(mavgCrossChanThetaPowCorr(50000,:,:))),colorbar();


proxChanThetaPowCorr = [];
tp = mean(log10(cys(ind,cfs>6&cfs<12,2)),2);%+mean(log10(cys(ind,cfs>6&cfs<12,2)),2);
for c = 1:32
    lp = mean(log10(lys(ind,lfs>6&lfs<12,c)),2);
    nind = nniz(lp)&nniz(tp);
    proxChanThetaPowCorr(c) = corr(lp(nind),tp(nind));
end

ind = [stc{'w&t-m-s'}];
ind = [stc{'p&t-m-s'}];
proxChanThetaRatioCorr = [];
tp = mean(log10(cys(ind,cfs>6&cfs<12,2)),2);%+mean(log10(cys(ind,cfs>6&cfs<12,2)),2);
for c1 = 1:32
    for c2 = c1+1:32    
    lp = mean(log10(lys(ind,lfs>6&lfs<12,c1)),2)./mean(log10(lys(ind,lfs>6&lfs<12,c2)),2);
    nind = nniz(lp)&nniz(tp);
    proxChanThetaRatioCorr(c1,c2) = corr(lp(nind),tp(nind));
    end
end

figure,
imagesc(proxChanThetaRatioCorr);
colormap('jet');

%%%>>>

%%%<<< COMPUTE lfp features ------------------------------------------------------------------------
flys = lys.copy();
flys.filter('ButFilter',4,1,'low');
flys.data(flys.data<0) = 0.001;

% jg05-20120310
chanOri = 2;
chanPyr = 6;
chanRad = 9;
chanLm = 17;

% jg05-20120312
chanOri = 3;
chanPyr = 7;
chanRad = 13;
chanLm = 19;


% plot theta power for CA1{ori,pyr,rad,lm} and DG
figure();
sax = tight_subplot(3,1,[0.01,0.01],[0.1,0.1],[0.1,0.1]);
axes(sax(1));
imagesc(lts,1:32,sq(mean(log10(lys(:,lfs>5&lfs<12,:)),2))');
axis(sax(1),'ij');
colormap('jet');
caxis([4,6.7]);
%caxis([-3,3]);
axes(sax(2));
imagesc(lts,1:32,sq(mean(log10(lys(:,lfs>20&lfs<28,:)),2))');
axis(sax(2),'ij');
colormap('jet');
%caxis([-3,3]);
caxis([4,6.7]);
axes(sax(3));
plotSTC(Trial.stc,1,'text',states,stateColors);
ylim(sax(end),[0,9]);
linkaxes(sax,'x');



% CA1LM Theta power 
thp = lys.copy();
thp.data = mean(lys(:,lfs>5&lfs<10,chanLm),2)./mean(lys(:,lfs<4|(lfs<13&lfs>10),chanLm),2);
%thp.data(nniz(thp)) = ButFilter(thp(nniz(thp)),4,[0.5]./(lys.sampleRate./2),'low');
thp.resample(vxy);

% CA1 rad/lm theta power ratio
thr = lys.copy();
thr.data = mean(log10(lys(:,lfs>6&lfs<12,chanRad)),2)./mean(log10(lys(:,lfs>6&lfs<12,chanLm)),2);
thr.resample(vxy);

% CA1 pyr/lm theta power ratio
thl = lys.copy();
thl.data = mean(log10(lys(:,lfs>6&lfs<12,chanPyr)),2)./mean(log10(lys(:,lfs>6&lfs<12,chanLm)),2);
thl.resample(vxy);

% CA1 pyr/lm theta power ratio
% $$$ thpPyr = lys.copy();
% $$$ thpPyr.data = mean(lys(:,lfs>5&lfs<10,2),2)./mean(lys(:,lfs<3|(lfs<13&lfs>10),2),2);
% $$$ thpPyr.data(nniz(thpPyr)) = ButFilter(thpPyr(nniz(thpPyr)),4,[0.5]./(lys.sampleRate./2),'low');
% $$$ thpPyr.resample(vxy);
% $$$ thpRad = lys.copy();
% $$$ thpRad.data = mean(lys(:,lfs>5&lfs<10,chanRad),2)./mean(lys(:,lfs<3|(lfs<13&lfs>10),chanRad),2);
% $$$ thpRad.data(nniz(thpRad)) = ButFilter(thpRad(nniz(thpRad)),4,[0.5]./(lys.sampleRate./2),'low');
% $$$ thpRad.resample(vxy);

% MEAN low frequency power <15Hz or <20Hz
lwp = lys.copy();
lwpFreq = 15;
lwpMean = mean(log10(lys(stc{'x'},lfs<lwpFreq,chanPyr)),2);
lwpMean = mean(lwpMean(nniz(lwpMean)));
lwpStd = mean(log10(lys(stc{'x'},lfs<lwpFreq,chanPyr)),2);
lwpStd = std(lwpStd(nniz(lwpStd)));
lwp.data = nunity(mean(log10(lys(:,lfs<lwpFreq,chanPyr)),2),[],...
                  lwpMean,lwpStd);
lwp.resample(vxy);
%%%>>>

figure();
subplot(211);
hold('on');
plot([1:size(ufrInt,1)]/lfp.sampleRate,sum(ufrInt(:,mpv>2&mpv<3),2)./sum(mpv>2&mpv<3));
plot([1:size(ufrInt,1)]/lfp.sampleRate,sum(ufrInt(:,mpv>3.1),2)./sum(mpv>3.1));
plot([1:size(ufrInt,1)]/lfp.sampleRate,sum(fufrInt(:,mpv>2&mpv<3),2)./sum(mpv>2&mpv<3));
plot([1:size(ufrInt,1)]/lfp.sampleRate,sum(fufrInt(:,mpv>3.1),2)./sum(mpv>3.1));
plot(lts,thp,'r')
plot(lts,lwp,'m')
subplot(212);
hold('on');
plotSTC(Trial.stc,1,'text',states,stateColors);
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
ylim([1,9]);
linkaxes(findobj(gcf(),'Type','Axes'),'x');

intFet = fufrInt.copy();
intFet.data = sum(fufrInt(:,mpv>3&mpv<4),2)./sum(mpv>3&mpv<4)-sum(fufrInt(:,mpv>4.1),2)./sum(mpv>4.1);
intFet.resample(vxy);

norm = 'xprob';
norm = '';


figure,
ind = stc{'s'};
subplot(221);
hist2([log10(vxy(ind,1)),lwp(ind,1)],linspace(-3,1,50),linspace(-3,3,50),norm);
subplot(222);
hist2([log10(vxy(ind,2)),lwp(ind,1)],linspace(-3,1,50),linspace(-3,3,50),norm);
subplot(223);
hist2([log10(vxy(ind,1)),thpLm(ind,1)],linspace(-3,1,50),linspace(0,5,50),norm);
subplot(224);
hist2([log10(vxy(ind,2)),thpLm(ind,1)],linspace(-3,1,50),linspace(0,5,50),norm);


figure,
ind = stc{'s'};
subplot(221);
hist2([log10(vxy(ind,1)),lwp(ind,1)],linspace(-3,1,50),linspace(-3,3,50),norm);
subplot(222);
hist2([log10(vxy(ind,2)),lwp(ind,1)],linspace(-3,1,50),linspace(-3,3,50),norm);
subplot(223);
hist2([thpLm(ind,1),lwp(ind,1)],linspace(0,5,50),linspace(-3,3,50),norm);
subplot(224);
hist2([log10(vxy(ind,2)),thpLm(ind,1)],linspace(-3,1,50),linspace(0,5,50),norm);


dc = accumulate_decoding_vars(Trial,                               ...
                              units{trialId},                      ...
                              sessionList(trialId).thetaRefGeneral,...
                              phzCorrection(trialId),              ...
                              headRotation{trialId},               ...
                              hbangCorrection{trialId});




thpLmBin = linspace(0,5,20);
thpLmInd = discretize(thpLm.data,thpLmBin);

vxyBin = linspace(-2.5,2,20);
vxyInd = discretize(log10(vxy(:,2)),vxyBin);

lwpLim = [-4,4];
lwpBin = linspace([lwpLim,30]);
lwpCtr = mean([lwpBin(1:end-1);lwpBin(2:end)]);
lwpInd = discretize(lwp.data,lwpBin);

thlBin = linspace(0.7,0.925,30);
thlInd = zeros(size(thl));
thlInd(nniz(thl)) = discretize(thl(nniz(thl)),thlBin);


chanRad = 13;chanLm = 19;
chanRad = 3;chanLm = 10;
chanRad = 4;chanLm = 7;
%chanRad = 3;chanLm = 12; thrBin = linspace(0.825,1.1,30);
chanRad = 8;chanLm = 22; thrBin = linspace(0.775,0.95,30);
%chanRad = 14;chanLm = 22; thrBin = linspace(0.825,1.1,30);
thr = lys.copy();
thr.data = mean(log10(lys(:,lfs>6&lfs<12,chanRad)),2)./mean(log10(lys(:,lfs>6&lfs<12,chanLm)),2);
thr.resample(vxy);
thrInd = zeros(size(thr));
thrInd(nniz(thr)) = discretize(thr(nniz(thr)),thrBin);

chanRad = 13;chanLm = 19;
chanRad = 3;chanLm = 10;
%chanRad = 13;chanLm = 21;thlBin = linspace(0.84,1,30);
chanRad = 4;chanLm = 7;  thlBin = linspace(0.925,1.1,30);
%chanRad = 8;chanLm = 22; thlBin = linspace(0.775,0.95,30);
thl = lys.copy();
thl.data = mean(log10(lys(:,lfs>6&lfs<12,chanRad)),2)./mean(log10(lys(:,lfs>6&lfs<12,chanLm)),2);
thl.resample(vxy);
thlInd = zeros(size(thl));
thlInd(nniz(thl)) = discretize(thl(nniz(thl)),thlBin);

ppdlim = [-2,0];
ppdBin = linspace([ppdlim,30]);
ppdInd = zeros(size(proxPhzDiff));
ppdInd(nniz(proxPhzDiff)) = discretize(proxPhzDiff(nniz(proxPhzDiff)),ppdBin);


%xFet = lwp(:,1);   xInd = lwpInd; xBin = lwpBin;
%xFet = thl(:,1);   xInd = thlInd; xBin = thlBin;
xFet = proxPhzDiff(:,1);   xInd = ppdInd; xBin = ppdBin;
yFet = thr(:,1);   yInd = thrInd; yBin = thrBin;
errdist = nan([size(xyz,1),1]);
%errdist(dc.ind) = sqrt(sum(dc.ecom(:,1:2).^2,2)); edlim = [0,400];
%errdist(dc.ind) = dc.esax(:,1);;
errdist(dc.ind) = dc.ecom(:,1);; edlim = [-100,100];
%errdist(dc.ind) = dc.ecom(:,2);; edlim = [-100,100];


indGrps = {dc.stcm(:,8)==8&dc.stcm(:,1)==1&dc.ucnt>2,...
           dc.stcm(:,8)==8&dc.stcm(:,1)~=1&dc.ucnt>2,...
           (dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)==1&dc.ucnt>2,...
           (dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)~=1&dc.ucnt>2,...
           (dc.stcm(:,6)==6|dc.stcm(:,4)==4)&dc.stcm(:,1)==1&dc.ucnt>2,...
           (dc.stcm(:,6)==6|dc.stcm(:,4)==4)&dc.stcm(:,1)~=1&dc.ucnt>2};
indLbls = {'sit theta','sit not theta',...           
           'loc',      'loc not theta',...
           'pause',    'pause not theta'};

figure
sax = tight_subplot(numel(indGrps),3,[0.05,0.05],[0.05,0.05],[0.05,0.05]);
for s = 1:numel(indGrps)
ind = false([size(xyz,1),1]);
ind(dc.ind(indGrps{s})) = true;
ind(~nniz(xInd)|~nniz(yInd)) = false;
ind = find(ind);
outc = hist2([xFet(ind,1),yFet(ind,1)],xBin,yBin,'');
%out = accumarray([xInd(ind),yInd(ind)],vxy(ind,2),[numel(xBin)-1,numel(yBin)-1],@mean);
out = accumarray([xInd(ind),yInd(ind)],errdist(ind,1),[numel(xBin)-1,numel(yBin)-1],@mean);
out(outc(:)<40) = nan;
%outs = accumarray([xInd(ind),yInd(ind)],vxy(ind,2),[numel(xBin)-1,numel(yBin)-1],@std);
outs = accumarray([xInd(ind),yInd(ind)],errdist(ind,1),[numel(xBin)-1,numel(yBin)-1],@std);
outs(outc(:)<40) = nan;
axes(sax((s-1)*3+1));
    set(pcolor(mean([xBin(2:end);xBin(1:end-1)]),mean([yBin(2:end);yBin(1:end-1)]),out'),'EdgeColor','none');
    axis('xy');
    colormap('jet');
    colorbar();
    %caxis([0,600]);    
    %caxis([0,300]);    
    caxis(edlim);
    %caxis([0,30]);    
    %Lines([],0.925,'r');
% $$$     Lines(0,[],'r');    
    %line(xBin([1,end]),yBin([1,end]),'Color','m');
axes(sax((s-1)*3+2));    
    imagesc(xBin,yBin,outs');
    axis('xy');
    colormap('jet');
    colorbar();
    caxis([0,300]);    
    %caxis([0,150]);        
    %caxis([0,30]);        
    %line(xBin([1,end]),yBin([1,end]),'Color','m');    
axes(sax((s-1)*3+3));    
    set(pcolor(mean([xBin(2:end);xBin(1:end-1)]),mean([yBin(2:end);yBin(1:end-1)]),outc'),'EdgeColor','none');
    %Lines([],0.925,'r');
% $$$     Lines(0,[],'r');    
%line(xBin([1,end]),yBin([1,end]),'Color','m');    
    title(indLbls{s});
end

% automatically find vetor and midpoint in pause to transform along the direction with largest change
s = find(~cellfun(@isempty,regexp(indLbls,'^pause$')));
%s = find(~cellfun(@isempty,regexp(indLbls,'^loc$')));
% FILTER inds for state and other criteria specified in indGrps
ind = false([size(xyz,1),1]);
ind(dc.ind(indGrps{s})) = true;
ind(~nniz(xInd)|~nniz(yInd)) = false;
ind = find(ind);
% COMPUTE count and conditional mean
outc = hist2([xFet(ind,1),yFet(ind,1)],xBin,yBin,'');
out = accumarray([xInd(ind),yInd(ind)],errdist(ind,1),[numel(xBin)-1,numel(yBin)-1],@mean);
out(outc(:)<30) = nan;
% COMPUTE xy grid
xBinCntr = mean([xBin(1:end-1);xBin(2:end)])';
yBinCntr = mean([yBin(1:end-1);yBin(2:end)])';
binGrid = cell([1,2]);
[binGrid{:}] = ndgrid(xBinCntr,yBinCntr);
% COMPUTE xy center
xyCntr = mean([binGrid{1}(outc(:)<30),binGrid{2}(outc(:)<30)]);

binGrid = {binGrid{1}-xyCntr(1),binGrid{2}-xyCntr(2)};

% COMPUTE maximum variance of conditional mean as a function of grid rotation around xyCntr
% ROTATE grid points and compute
gridDiagDist = sqrt(sum(([binGrid{1}(1,1),binGrid{2}(1,1)]-[binGrid{1}(2,2),binGrid{2}(2,2)]).^2));
rotate_grid = @(g,a) multiprod(g,[cos(a),-sin(a);sin(a),cos(a)],[2],[1,2]);

alphaSteps = -3*pi/2:0.01:3*pi/2;
xyrVar = nan(size(alphaSteps));
for alpha = 1:numel(alphaSteps)
tempGrid = rotate_grid([binGrid{1}(:),binGrid{2}(:)],alphaSteps(alpha));
tempInds = discretize(tempGrid(:,1),...
                      (min(tempGrid(:,1))-gridDiagDist) ...
                      : gridDiagDist ...
                      : max(tempGrid(:,1))+gridDiagDist);
tempOut = accumarray(tempInds,out(:),[max(tempInds),1],@(x) mean(x,'omitnan'));
xyrVar(alpha) = sum((tempOut(nniz(tempOut))-mean(tempOut(nniz(tempOut)))).^2)./sum(nniz(tempOut));
end

xyrVarSmth = RectFilter(xyrVar',31,3);
figure();
hold('on');
plot(alphaSteps,xyrVar);
plot(alphaSteps,xyrVarSmth);

mind = LocalMinima(-xyrVarSmth,31,0);

bestRotation = alphaSteps(mind(find(abs(alphaSteps(mind))-min(abs(alphaSteps(mind)))==0)));


plot(alphaSteps,RectFilter(xyrVar',61));

tempGrid = rotate_grid([binGrid{1}(:),binGrid{2}(:)],bestRotation);
figure();
surf(reshape(tempGrid(:,1),[29,29]),reshape(tempGrid(:,2),[29,29]),out);
colormap('jet');

thn = copy(thr);
thn.data = rotate_grid(bsxfun(@minus,[thl.data,thr.data],xyCntr),bestRotation);
thn.data = thn.data(:,1);




% $$$ x = tbins;
% $$$ y = RectFilter(mccg(:,1,2)-mean(mccg(:,1,2)),5)'-RectFilter(mccg(:,1,2)-mean(mccg(:,1,2)),11)';
% $$$ yr = (max(y)-min(y));                               % Range of ?y?
% $$$ yz = y-max(y)+(yr/2);
% $$$ zx = x(yz .* circshift(yz,[0 1]) <= 0);     % Find zero-crossings
% $$$ per = 2*mean(diff(zx));                     % Estimate period
% $$$ dy = diff(y);
% $$$ dzx = dy(yz .* circshift(yz,[0 1]) <= 0);
% $$$ pinds = [(find(zx<0,1,'last')),(find(zx>0,1,'first'))];
% $$$ ccgphz = (zx(pinds(find(dzx(pinds)>0,1,'first'))));
% $$$ 
% $$$ ym = mean(y);                               % Estimate offset
% $$$ p = polyfit(tbins((((numel(y)-1)/2)+2):end),mean([fliplr(y(1:((numel(y)-1)/2)));y((((numel(y)-1)/2)+2):end)]),1);
% $$$ [xp,xpi] = min(abs(zx));
% $$$ xp = xp.*sign(zx(xpi));
% $$$ b = [ per;  ccgphz; 150; 1];
% $$$ 
% $$$ fit = @(b,x) b(4).*exp(-0.5.*(x.^2./b(3).^2)).*  ...
% $$$               ((cos(2*pi*(x + b(2))./b(1) )));    % Function to fit
% $$$ fcn = @(b) sum((fit(b,x) - y).^2);                              % Least-Squares cost function
% $$$ s(:,intInd) = fminsearch(fcn, b,struct('MaxIter',10000)) ;                      % Minimise Least-Squares
% $$$ xp = linspace(min(x),max(x));






% Manually find vetor and midpoint in pause to transform along the direction with largest change
thd = copy(thr);
thd.data = bsxfun(@minus,[thl.data,thr.data],[1.019,0.8565]);
([0.9763,0.9047]-[1.019,0.8565])./
1.055,0.8082
figure,
hold('on');
histogram(thr([Trial.stc{'lloc'}]),thrBin,'EdgeColor','none')
histogram(thr([Trial.stc{'hloc'}]),thrBin,'EdgeColor','none')
histogram(thr([Trial.stc{'rear'}]),thrBin,'EdgeColor','none')

figure,
hold('on');

histogram(thr([Trial.stc{'hpause&theta'}]),thrBin,'EdgeColor','none')
histogram(thr([Trial.stc{'rear-theta'}]),thrBin,'EdgeColor','none')
histogram(thr([Trial.stc{'rear&theta'}]),thrBin,'EdgeColor','none')
legend({'lp','hp','r','rt'})




ind = false([size(xyz,1),1]);
%ind(dc.ind(dc.stcm(:,8)==8&dc.stcm(:,1)~=1)) = true;
ind(dc.ind(dc.stcm(:,2)==2&dc.stcm(:,1)==1)) = true;
ind(dc.ind((dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)==1)) = true;
ind(dc.ind((dc.stcm(:,4)==4|dc.stcm(:,6)==6)&dc.stcm(:,1)==1)) = true;
%ind(dc.ind(dc.stcm(:,3)==3&dc.stcm(:,1)==1)) = true;
figure,plot(tpow(ind,1),errdist(ind,1),'.')

figure();
for c = 1:32
    subplot(6,6,c);
    hist2([tpow(ind,c),errdist(ind,1)],linspace(4,7,40),linspace(0,300,40),'yprob');
    %caxis([0,0.2])    
end

figure();
for c = 1:32
    subplot(6,6,c);
    hist2([tpow(ind,c)./tpow(ind,27),errdist(ind,1)],linspace(0.75,1.2,50),linspace(0,800,40),'xprob');
    caxis([0,0.1])    
end


errdistBin = linspace(0,800,21);
errdistCtr = mean([errdistBin(1:end-1);errdistBin(2:end)]);
errdistInd = discretize(errdist,errdistBin);
tpowm = [];
tpows = [];
b = 1;
for band = lfs'
tpow = copy(lys);
tpow.data = sq(mean(log10(lys(:,lfs==band,:)),2));
tpow.resample(xyz);
for bin = 1:20,
tpowm(:,bin,b) = mean(tpow(ind&errdistInd==bin,:),'omitnan');
tpows(:,bin,b) = std(tpow(ind&errdistInd==bin,:),'omitnan');
end
b = b+1;
end

figure();
for b = 1:size(tpowm,3)
subplot2(4,size(tpowm,3),1,b);
    imagesc(errdistCtr,1:32,tpowm(:,:,b));    
subplot2(4,size(tpowm,3),2,b);
    imagesc(errdistCtr,1:32,bsxfun(@rdivide,tpowm(:,:,b),sum(tpowm(:,:,b),1)));
subplot2(4,size(tpowm,3),3,b);
    imagesc(errdistCtr,1:32,bsxfun(@rdivide,tpowm(:,:,b),sum(tpowm(:,:,b),2)));
subplot2(4,size(tpowm,3),4,b);    
    imagesc(errdistCtr,1:32,tpows(:,:,b));
end


bins = 1:20;
for b = 1:size(tpowm,3)
    for chan = 1:32
        try
    p = polyfit(bins(~isnan(tpowm(chan,:,b)))',tpowm(chan,~isnan(tpowm(chan,:,b)),b)',1);        
    tpowc(chan,b) = p(1);
    end
    end
end

figure,
imagesc(lfs,1:32,tpowc)
colormap('jet');
caxis(max(abs(caxis)).*[-1,1]);

c = 9;
f = 5;
figure()
    plot(errdistCtr,tpowm(c,:,f),'-+')
    hold('on');
    plot(errdistCtr,tpowm(c,:,f)+tpows(c,:,f),'-+r')
    plot(errdistCtr,tpowm(c,:,f)-tpows(c,:,f),'-+r')

c = 9;

figure()
hold('on');
for f = 1:numel(lfs);
    plot(errdistCtr,tpowm(c,:,f),'-')
end

    

figure();
hist2([thr(ind),errdist(ind)],linspace(0.86,0.96,50),linspace(-300,300,40),'yprob');
caxis([0,0.08]);
Lines([],50,'r');


tpow = copy(lys);
tpow.data = sq(mean(log10(lys(:,lfs>5&lfs<12,:)),2));
tpow.resample(xyz);

gpow = copy(lys);
gpow.data = sq(mean(log10(lys(:,lfs>30,:)),2));
gpow.resample(xyz);


% THR vs ERRDIST
ind = false([size(xyz,1),1]);
ind(dc.ind(dc.stcm(:,8)==8&dc.stcm(:,1)~=1&dc.ucnt>2&dc.ucnt>2)) = true;
ind(dc.ind(dc.stcm(:,2)==2&dc.stcm(:,1)==1)) = true;
ind(dc.ind((dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)==1 ...
           &dc.ucnt>2&sqrt(sum(xyz(dc.ind,'hcom',[1,2]).^2,3))<350)) = true;
ind(dc.ind((dc.stcm(:,4)==4|dc.stcm(:,6)==6)&dc.stcm(:,1)==1&dc.ucnt>2)) = true;
ind(dc.ind((dc.stcm(:,6)==6)&dc.stcm(:,1)==1&dc.ucnt>2&sqrt(sum(xyz(dc.ind,'hcom',[1,2]).^2,3))<350)) = true;
ind(dc.ind((dc.stcm(:,4)==4)&dc.stcm(:,1)==1)) = true;
%ind(dc.ind(dc.stcm(:,3)==3&dc.stcm(:,1)==1)) = true;

chanRad = 13;chanLm = 19;
chanRad = 3;chanLm = 10;
chanRad = 4;chanLm = 7;
chanRad = 8;chanLm = 22;
thr = lys.copy();
thr.data = mean(log10(lys(:,lfs>6&lfs<12,chanRad)),2)./mean(log10(lys(:,lfs>6&lfs<12,chanLm)),2);
thr.resample(vxy);


norm = 'yprob';
%norm = '';
clim = [0,0.15];
%clim = 'auto';
thrlim = [0.75,0.95];
thrlim = [0.75,1.1];
thrlim = [0.95,1.075];
%thrlim = [0.875,0.98];
tpowRngLm = [5.5,7];
tpowRngLm = [4,5.7];
tpowRngRad = [4.5,6.5];
%tpowRngRad = [4.5,6.5];
efet = 'ecom';
figure();
errdist = nan([size(xyz,1),1]);
errdist(dc.ind) = sqrt(sum(dc.(efet)(:,1:2).^2,2));
subplot2(3,3,1,1);
    hist2([thr(ind),errdist(ind)],linspace([thrlim,30]),linspace(0,300,20),norm);
    [p,s]  = polyfit(thr(ind),errdist(ind),1);
    [r,pval] = corr(thr(ind),errdist(ind));
    line(thrlim,polyval(p,thrlim),'color','r');
    caxis(clim);
    title(['R^2 = ',num2str(round(r,3)),' pval = ' num2str(round(pval,3))]);
subplot2(3,3,1,2);
    hist2([tpow(ind,chanRad),errdist(ind)],linspace([tpowRngRad,30]),linspace(0,300,20),norm);
    [p,s]  = polyfit(tpow(ind,chanRad),errdist(ind),1);
    [r,pval] = corr(tpow(ind,chanRad),errdist(ind));
    line(tpowRngRad,polyval(p,tpowRngRad),'color','r');
    caxis(clim);
    title(['R^2 = ',num2str(round(r,3)),' pval = ' num2str(round(pval,3))]);
subplot2(3,3,1,3);
    hist2([tpow(ind,chanLm),errdist(ind)],linspace([tpowRngLm,30]),linspace(0,300,20),norm);
    [p,s]  = polyfit(tpow(ind,chanLm),errdist(ind),1);
    [r,pval] = corr(tpow(ind,chanLm),errdist(ind));
    line(tpowRngLm,polyval(p,tpowRngLm),'color','r');
    caxis(clim);
    title(['R^2 = ',num2str(round(r,3)),' pval = ' num2str(round(pval,3))]);
errdist = nan([size(xyz,1),1]);    
errdist(dc.ind) = dc.(efet)(:,1);;
subplot2(3,3,2,1);
    hist2([thr(ind),errdist(ind)],linspace([thrlim,30]),linspace(-200,200,20),norm);
    [p,s]  = polyfit(thr(ind),errdist(ind),1);
    [r,pval] = corr(thr(ind),errdist(ind));
    line(thrlim,polyval(p,thrlim),'color','r');
    caxis(clim);
    title(['R^2 = ',num2str(round(r,3)),' pval = ' num2str(round(pval,3))]);
subplot2(3,3,2,2);
    hist2([tpow(ind,chanRad),errdist(ind)],linspace([tpowRngRad,30]),linspace(-200,200,20),norm);
    [p,s]  = polyfit(tpow(ind,chanRad),errdist(ind),1);
    [r,pval] = corr(tpow(ind,chanRad),errdist(ind));
    line(tpowRngRad,polyval(p,tpowRngRad),'color','r');
    caxis(clim);
    title(['R^2 = ',num2str(round(r,3)),' pval = ' num2str(round(pval,3))]);
subplot2(3,3,2,3);
    hist2([tpow(ind,chanLm),errdist(ind)],linspace([tpowRngLm,30]),linspace(-200,200,20),norm);
    [p,s]  = polyfit(tpow(ind,chanLm),errdist(ind),1);
    [r,pval] = corr(tpow(ind,chanLm),errdist(ind));
    line([tpowRngLm],polyval(p,[tpowRngLm]),'color','r');
    caxis(clim);
    title(['R^2 = ',num2str(round(r,3)),' pval = ' num2str(round(pval,3))]);
errdist = nan([size(xyz,1),1]);    
errdist(dc.ind) = dc.(efet)(:,2);;
subplot2(3,3,3,1);
    hist2([thr(ind),errdist(ind)],linspace([thrlim,30]),linspace(-200,200,20),norm);
    [p,s]  = polyfit(thr(ind),errdist(ind),1);
    [r,pval] = corr(thr(ind),errdist(ind));
    line(thrlim,polyval(p,thrlim),'color','r');
    caxis(clim);
    title(['R^2 = ',num2str(round(r,3)),' pval = ' num2str(round(pval,3))]);
subplot2(3,3,3,2);
    hist2([tpow(ind,chanRad),errdist(ind)],linspace([tpowRngRad,30]),linspace(-200,200,20),norm);
    [p,s]  = polyfit(tpow(ind,chanRad),errdist(ind),1);
    [r,pval] = corr(tpow(ind,chanRad),errdist(ind));
    line(tpowRngRad,polyval(p,tpowRngRad),'color','r');
    caxis(clim);
    title(['R^2 = ',num2str(round(r,3)),' pval = ' num2str(round(pval,3))]);
subplot2(3,3,3,3);
    hist2([tpow(ind,chanLm),errdist(ind)],linspace([tpowRngLm,30]),linspace(-200,200,20),norm);
    [p,s]  = polyfit(tpow(ind,chanLm),errdist(ind),1);
    [r,pval] = corr(tpow(ind,chanLm),errdist(ind));
    line(tpowRngLm,polyval(p,tpowRngLm),'color','r');
    caxis(clim);
    title(['R^2 = ',num2str(round(r,3)),' pval = ' num2str(round(pval,3))]);    

% COMPUTE correlation bettween laminar theta power ratios vs errdist
rThrEcom = nan([32,32]);
pvalThrEcom = nan([32,32]);
errdist = nan([size(xyz,1),1]);    
errdist(dc.ind) = dc.(efet)(:,1);;
for c1 = 1:32,
    for c2 = c1:32,   
        if c1==c2
            thp = lys.copy();
            thp.data = mean(log10(lys(:,lfs>6&lfs<12,c1)),2);
            thp.resample(vxy);
            [rThrEcom(c2,c1),pvalThrEcom(c2,c1)] = corr(thp(ind),errdist(ind));
        else
            thr = lys.copy();
            thr.data = mean(log10(lys(:,lfs>6&lfs<12,c1)),2)./mean(log10(lys(:,lfs>6&lfs<12,c2)),2);
            thr.resample(vxy);
            [rThrEcom(c2,c1),pvalThrEcom(c2,c1)] = corr(thr(ind),errdist(ind));
        end
    end
end

figure,
subplot(121);
imagesc(rThrEcom');
colormap('jet');
subplot(122);
imagesc(1./pvalThrEcom');
caxis([20,100])
title('Theta Pwr Ratio vs EgoFwd (pause)');
title('Theta Pwr Ratio vs EgoLat (pause)');

figure,plot(nunity(tpow(ind,chanLm)),nunity(tpow(ind,chanRad)),'.');

figure,plot3(tpow(ind,17),tpow(ind,9),tpow(ind,9)./tpow(ind,17),'.');

errdist = nan([size(xyz,1),1]);    
errdist(dc.ind) = sqrt(sum(dc.ecom(:,1:2).^2,2));
errdist(dc.ind) = dc.ecom(:,2);
errdist(dc.ind) = dc.ecom(:,1);
tinds = [discretize(tpow(ind,chanLm),linspace(5.2,6.5,11)),discretize(tpow(ind,chanRad),linspace(4.8,6.2,11))];
evals = errdist(ind);
evals = log10(vxy(ind,2));
ninds = nniz(tinds)&nniz(evals);
out =accumarray(tinds(ninds,:),evals(ninds),[10,10],@mean);
out(out==0) = nan;

figure();
    set(pcolor(linspace(5.2,6.5,10),linspace(4.8,6.2,10),out'),'EdgeColor','none');
    axis('xy');
    colormap('jet');

    
figure();
hist2([thr(ind),log10(vxy(ind,1))],linspace(0.86,0.96,40),linspace(-2,1.6,40),'');

figure();
hist2([errdist(ind),log10(vxy(ind,1))],linspace(-300,300,40),linspace(-2,1.6,40),'xprob');




mazeCntrDist = sqrt(sum(xyz(:,'hcom',[1,2]).^2,3));

ind = false([size(xyz,1),1]);
ind(dc.ind(dc.stcm(:,8)==8&dc.stcm(:,1)~=1&dc.ucnt>2&dc.ucnt>2)) = true;
ind(dc.ind(dc.stcm(:,2)==2&dc.stcm(:,1)==1)) = true;
ind(dc.ind((dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)==1 ...
           & dc.ucnt>2 & mazeCntrDist(dc.ind)<350 )) = true;
ind(dc.ind((dc.stcm(:,4)==4|dc.stcm(:,6)==6)&dc.stcm(:,1)==1 & dc.ucnt>2 &  mazeCntrDist(dc.ind)<350)) = true;
ind(dc.ind((dc.stcm(:,6)==6)&dc.stcm(:,1)==1 & dc.ucnt>2 & mazeCntrDist(dc.ind)<350)) = true;
ind(dc.ind((dc.stcm(:,4)==4)&dc.stcm(:,1)==1)) = true;
ind(dc.ind((dc.stcm(:,6)==6)&dc.stcm(:,1)==1)) = true;
%ind(dc.ind(dc.stcm(:,3)==3&dc.stcm(:,1)==1)) = true;

thp = tpow(:,chanRad);
thp = tpow(:,chanLm);
thplim = [5,7];
%thplim = [4,6];
thpBin = linspace([thplim,40]);
thpInd = discretize(thp,thpBin);

figure
out = hist2([thp(ind),erf(ind)],thpBin,erfBin,norm);
out(out(:)==0) = nan;
set(pcolor(thpBin(1:end-1),erfBin(1:end-1),out'),'EdgeColor','none');    
[p,s]  = polyfit(thp(ind),erf(ind),1);    [r,pval] = corr(thp(ind),erf(ind));
line(thplim,polyval(p,thplim),'color','r');
caxis(clim);
colorbar();        
title(['R^2 = ',num2str(round(r,3)),' pval = ' num2str(round(pval,3))]);

chanPyr = 7;
chanRad = 12;
chanLm = 21;

norm = 'yprob';
norm = '';
clim = [0,0.15];
clim = 'auto';
thnlim = [-0.12,0.12];
%thnlim = [-0.05,0.1];
thnBin = linspace([thnlim,40]);
thnInd = discretize(thn.data,thnBin);


erd = nan([size(xyz,1),1]);    
erd(dc.ind) = sqrt(sum(dc.(efet).^2,2));
erdlim = [0,600];
erdBin = linspace([erdlim,40]);
erdInd = discretize(erd,erdBin);

erf = nan([size(xyz,1),1]);    
erf(dc.ind) = dc.(efet)(:,1);
erflim = [-200,300];
erfBin = linspace([erflim,40]);
erfInd = discretize(erf,erfBin);

erl = nan([size(xyz,1),1]);    
erl(dc.ind) = dc.(efet)(:,2);
erllim = [-200,200];
erlBin = linspace([erllim,40]);
erlInd = discretize(erl,erlBin);

mcdlim = [0,350];
%vxylim = [0,15];
vxylim = [0,30];
%thnlim = [0.875,0.98];
%tpowRngRad = [4.5,6.5];
efet = 'ecom';
nx = 6
ny = 2
figure(); 
subplot2(ny,nx,1,1);
    [out,xb,yb,pp] = hist2([thn(ind),erf(ind)],thnBin,erfBin,norm);
    outc = out;
    gind = nan([size(pp,1),1]);
    gind(nniz(pp)) = outc(sub2ind(size(outc),pp(nniz(pp),1),pp(nniz(pp),2)));
    gind(gind<=5) = nan;
    aind = nan([size(xyz,1),1]);
    aind(ind) = gind;
    aind = nniz(aind);
    out(outc(:)<=5) = nan;
    pax = pcolor(thnBin(1:end-1),erfBin(1:end-1),out');
    set(pax,'EdgeColor','none');    
    [p,s]  = polyfit(thn(ind),erf(ind),1);    [r,pval] = corr(thn(ind&aind),erf(ind&aind));
    line(thnlim,polyval(p,thnlim),'color','r');
    caxis(clim);
    colorbar();        
    title(['R^2 = ',num2str(round(r,3)),' pval = ' num2str(round(pval,3))]);
subplot2(ny,nx,1,2);
% cond mean mazeCntrDist ( thn, errfwd )
    nind = nniz([thnInd,erfInd]);
    out = accumarray([thnInd(ind&nind&aind),erfInd(ind&nind&aind)],...
                     mazeCntrDist(ind&nind&aind),...
                     [numel(thnBin)-1,numel(erfBin)-1],...
                     @(x) mean(x,'omitnan'));
    out(out(:)==0) = nan;
    set(pcolor(thnBin(1:end-1),erfBin(1:end-1),out'),'EdgeColor','none');
    caxis(mcdlim);
    colorbar();
    title(['cond mean mazeCntrDist']);
subplot2(ny,nx,1,3);
    nind = nniz([thnInd,erfInd]);
    out = accumarray([thnInd(ind&nind&aind),erfInd(ind&nind&aind)],...
                     vxy(ind&nind&aind,2),...
                     [numel(thnBin)-1,numel(erfBin)-1],...
                     @(x) mean(x,'omitnan'));
    out(out(:)==0) = nan;
    set(pcolor(thnBin(1:end-1),erfBin(1:end-1),out'),'EdgeColor','none');
    caxis(vxylim);
    colorbar();    
    title(['cond mean xy speed']);
subplot2(ny,nx,1,4);
    nind = nniz([thnInd,erfInd]);
    out = accumarray([thnInd(ind&nind&aind),erfInd(ind&nind&aind)],...
                     tpow(ind&nind&aind,chanPyr),...
                     [numel(thnBin)-1,numel(erfBin)-1],...
                     @(x) mean(x,'omitnan'));
    out(out(:)==0) = nan;
    set(pcolor(thnBin(1:end-1),erfBin(1:end-1),out'),'EdgeColor','none');
    %caxis(thrlim);
    colorbar();    
    title(['tpow pyr']);
subplot2(ny,nx,1,5);
    nind = nniz([thnInd,erfInd]);
    out = accumarray([thnInd(ind&nind&aind),erfInd(ind&nind&aind)],...
                     tpow(ind&nind&aind,chanRad),...
                     [numel(thnBin)-1,numel(erfBin)-1],...
                     @(x) mean(x,'omitnan'));
    out(out(:)==0) = nan;
    set(pcolor(thnBin(1:end-1),erfBin(1:end-1),out'),'EdgeColor','none');
    %caxis(thrlim);
    colorbar();    
    title(['tpow rad']);
subplot2(ny,nx,1,6);
    nind = nniz([thnInd,erfInd]);
    out = accumarray([thnInd(ind&nind&aind),erfInd(ind&nind&aind)],...
                     tpow(ind&nind&aind,chanLm),...
                     [numel(thnBin)-1,numel(erfBin)-1],...
                     @(x) mean(x,'omitnan'));
    out(out(:)==0) = nan;
    set(pcolor(thnBin(1:end-1),erfBin(1:end-1),out'),'EdgeColor','none');
    %caxis(thrlim);
    colorbar();    
    title(['tpow lm']);
subplot2(ny,nx,2,1);
    [out,xb,yb,pp] = hist2([thn(ind),erl(ind)],thnBin,erlBin,norm);
    outc = out;
    gind = nan([size(pp,1),1]);
    gind(nniz(pp)) = outc(sub2ind(size(outc),pp(nniz(pp),1),pp(nniz(pp),2)));
    gind(gind<=5) = nan;
    aind = nan([size(xyz,1),1]);
    aind(ind) = gind;
    aind = nniz(aind);
    out(outc(:)<=5) = nan;
    set(pcolor(thnBin(1:end-1),erlBin(1:end-1),out'),'EdgeColor','none');    
    [p,s]  = polyfit(thn(ind),erl(ind),1);    [r,pval] = corr(thn(ind),erl(ind));
    line(thnlim,polyval(p,thnlim),'color','r');
    caxis(clim);
    colorbar();        
    title(['R^2 = ',num2str(round(r,3)),' pval = ' num2str(round(pval,3))]);
subplot2(ny,nx,2,2);
    nind = nniz([thnInd,erlInd]);
    out = accumarray([thnInd(ind&nind&aind),erlInd(ind&nind&aind)],...
                     mazeCntrDist(ind&nind&aind),...
                     [numel(thnBin)-1,numel(erlBin)-1],...
                     @(x) mean(x,'omitnan'));
    out(out(:)==0) = nan;
    set(pcolor(thnBin(1:end-1),erlBin(1:end-1),out'),'EdgeColor','none');
    caxis(mcdlim);    
    colorbar();        
    title(['cond mean mazeCntrDist']);
subplot2(ny,nx,2,3);
    nind = nniz([thnInd,erlInd]);
    out = accumarray([thnInd(ind&nind&aind),erlInd(ind&nind&aind)],...
                     vxy(ind&nind&aind,2),...
                     [numel(thnBin)-1,numel(erlBin)-1],...
                     @(x) mean(x,'omitnan'));
    out(out(:)==0) = nan;
    set(pcolor(thnBin(1:end-1),erlBin(1:end-1),out'),'EdgeColor','none');
    caxis(vxylim);
    colorbar();    
    title(['cond mean xy speed']);
subplot2(ny,nx,2,4);
    nind = nniz([thnInd,erlInd]);
    out = accumarray([thnInd(ind&nind&aind),erlInd(ind&nind&aind)],...
                     tpow(ind&nind&aind,chanPyr),...
                     [numel(thnBin)-1,numel(erlBin)-1],...
                     @(x) mean(x,'omitnan'));
    out(out(:)==0) = nan;
    set(pcolor(thnBin(1:end-1),erlBin(1:end-1),out'),'EdgeColor','none');
    %caxis(thrlim);
    colorbar();    
    title(['tpow pyr']);
subplot2(ny,nx,2,5);
    nind = nniz([thnInd,erlInd]);
    out = accumarray([thnInd(ind&nind&aind),erlInd(ind&nind&aind)],...
                     tpow(ind&nind&aind,chanRad),...
                     [numel(thnBin)-1,numel(erlBin)-1],...
                     @(x) mean(x,'omitnan'));
    out(out(:)==0) = nan;
    set(pcolor(thnBin(1:end-1),erlBin(1:end-1),out'),'EdgeColor','none');
    %caxis(thrlim);
    colorbar();    
    title(['tpow rad']);
subplot2(ny,nx,2,6);
    nind = nniz([thnInd,erlInd]);
    out = accumarray([thnInd(ind&nind&aind),erlInd(ind&nind&aind)],...
                     tpow(ind&nind&aind,chanLm),...
                     [numel(thnBin)-1,numel(erlBin)-1],...
                     @(x) mean(x,'omitnan'));
    out(out(:)==0) = nan;
    set(pcolor(thnBin(1:end-1),erlBin(1:end-1),out'),'EdgeColor','none');
    %caxis(thrlim);
    colorbar();    
    title(['tpow lm']);
colormap('jet');


chanOri = 4;
chanPyr = 7;
chanRad = 12;
chanLm  = 21;

tpowOri = tpow(:,chanOri);
tpowOrilim = [4,6.2];
tpowOriBin = linspace([tpowOrilim,15]);
tpowOriInd = discretize(tpowOri,tpowOriBin);
tpowPyr = tpow(:,chanPyr);
tpowPyrlim = [4,6.2];
tpowPyrBin = linspace([tpowPyrlim,15]);
tpowPyrInd = discretize(tpowPyr,tpowPyrBin);
tpowLm = tpow(:,chanLm);
tpowLmlim = [5,7];
tpowLmBin = linspace([tpowLmlim,15]);
tpowLmInd = discretize(tpowLm,tpowLmBin);

figure,
subplot(221);
    plot(tpow(ind,chanPyr),tpow(ind,chanLm),'.')
    outcLm = hist2([tpowPyr(ind),tpowLm(ind)],tpowPyrBin,tpowLmBin);
    xlim(tpowPyrlim);
    ylim(tpowLmlim);
subplot(223);
    plot(tpow(ind,chanPyr),tpow(ind,chanOri),'.')
    outcOri = hist2([tpowPyr(ind),tpowOri(ind)],tpowPyrBin,tpowOriBin);    
    xlim(tpowPyrlim);
    ylim(tpowOrilim);
subplot(222);
    nind = nniz([tpowLmInd,tpowPyrInd]);
    out = accumarray([tpowPyrInd(ind&nind),tpowLmInd(ind&nind)],...
                     erd(ind&nind,1),...
                     [numel(tpowPyrBin)-1,numel(tpowLmBin)-1],...
                     @(x) mean(x,'omitnan'));
    out(out(:)==0|outcLm(:)<40) = nan;
    set(pcolor(tpowPyrBin(1:end-1),tpowLmBin(1:end-1),out'),'EdgeColor','none');
    %caxis([-50,100]);
    caxis([0,400]);
    colormap('jet');
    colorbar();
subplot(224);
    nind = nniz([tpowOriInd,tpowPyrInd]);
    out = accumarray([tpowPyrInd(ind&nind),tpowOriInd(ind&nind)],...
                     erd(ind&nind,1),...
                     [numel(tpowPyrBin)-1,numel(tpowOriBin)-1],...
                     @(x) mean(x,'omitnan'));
    out(out(:)==0|outcOri(:)<40) = nan;
    set(pcolor(tpowPyrBin(1:end-1),tpowOriBin(1:end-1),out'),'EdgeColor','none');
    %caxis([-50,100]);
    caxis([0,400]);
    colormap('jet');
    colorbar();    


% plot spectrograms for CA1{ori,pyr,rad,lm} and DG
figure();
sax = tight_subplot(2,1,[0.01,0.01],[0.1,0.1],[0.1,0.1]);
axes(sax(1));
hold('on');
plot(lts,sq(mean(log10(lys(:,lfs>5&lfs<12,[3])),2)));
plot(lts,sq(mean(log10(lys(:,lfs>5&lfs<12,[6])),2)));
plot(lts,sq(mean(log10(lys(:,lfs>5&lfs<12,[8])),2)));
plot(lts,sq(mean(log10(lys(:,lfs>5&lfs<12,[13])),2)));
plot([1:size(vxy,1)]./vxy.sampleRate,errdist(:,1)./1000+4);
%plot(lts,mean(sq(mean(log10(lys(:,lfs>20,[2,6,13])),2)),2)');
plot([1:size(vxy,1)]./vxy.sampleRate,vxy(:,2)./20+4);
plot([1:size(vxy,1)]./vxy.sampleRate,vxy(:,1)./20+4);
legend({'Pyr','Rad','Lm','DG','err','hv','bv'})
axes(sax(2));
plotSTC(Trial.stc,1,'text',states,stateColors);
ylim(sax(end),[0,9]);
linkaxes(sax,'x');


figure();
sax = tight_subplot(6,1,[0.01,0.01],[0.1,0.1],[0.1,0.1]);
for c = 1:5
imagesc(sax(c),lts,lfs,log10(lys(:,:,c))');
axis(sax(c),'xy');
colormap('jet');
caxis([4,6.7]);
end
gca(sax(end));
plotSTC(Trial.stc,1,'text',states,stateColors);
hold(sax(end),'on');
plot([1:size(thr,1)]./thr.sampleRate,10.^thr(:,1)-7);
plot([1:size(thr,1)]./thr.sampleRate,10.^thr(:,1)-7);
ylim(sax(end),[0,9]);
linkaxes(findobj(gcf(),'Type','Axes'),'x');


figure
for p = 2:2:24
ind = false([size(xyz,1),1]);
ind(dc.ind(dc.stcm(:,8)==8&dc.stcm(:,1)~=1&(dc.iphz==p-1|dc.iphz==p))) = true;
ind(dc.ind(dc.ucnt<3)) = false;
ind(~nniz(xInd)|~nniz(yInd)) = false;
ind = find(ind);
out = accumarray([xInd(ind),yInd(ind)],errdist(ind),[numel(xBin),numel(yBin)],@mean);
outs = accumarray([xInd(ind),yInd(ind)],errdist(ind),[numel(xBin),numel(yBin)],@std);
subplot2(12,3,p/2,1);
    imagesc(xBin,yBin,out');
    axis('xy');
    colormap('jet');
    colorbar();
    caxis([0,600]);
subplot2(12,3,p/2,2);
    imagesc(xBin,yBin,outs');
    axis('xy');
    colormap('jet');
    colorbar();
    caxis([0,350]);    
subplot2(12,3,p/2,3);    
    hist2([xFet(ind,1),yFet(ind,1)],xBin,yBin,'');
end



figure;hist2([intFet(ind,1),lwp(ind,1)],linspace(-3,3,50),linspace(-3,3,50));
figure;hist2([thpLm(ind,1),intFet(ind,1)],linspace(0,5,50),linspace(-3,3,50),norm);


figure,
subplot(121);
hist2([thr(nniz(vxy),1),log10(vxy(nniz(vxy),2))],linspace(0.85,1.1,50),linspace(-2.5,1.7,50),norm);
subplot(122);
hist2([thr(nniz(vxy),1),log10(vxy(nniz(vxy),1))],linspace(0.85,1.1,50),linspace(-2.5,1.7,50),norm);


figure,
subplot(131);
hist2([thpPyr(ind,1),lwp(ind,1)],linspace(0,5,50),linspace(-3,3,50),norm);
xlabel('Theta Delta Ratio (pyr)');
ylabel('mean power < 20Hz (pyr)');
subplot(132);
hist2([thpRad(ind,1),lwp(ind,1)],linspace(0,5,50),linspace(-3,3,50),norm);
xlabel('Theta Delta Ratio (Rad)');
ylabel('mean power < 20Hz (pyr)');
subplot(133);
hist2([thpLm(ind,1),lwp(ind,1)],linspace(0,5,50),linspace(-3,3,50),norm);
xlabel('Theta Delta Ratio (Lm)');
ylabel('mean power < 20Hz (pyr)');



% Immobile place representation segmentation
% Good Vars
%  mean power in pyr bellow 20Hz - Why is it so flat?
%  theta power in the lm         - Is this single input processing
%  rate difference between asc and dsc interneurons - see above


% COMPUTE lfp spectra
% $$$ glfp = lfp.copy();
% $$$ glfp.data = glfp.data(:,2);
% $$$ specArgsTheta = struct('nFFT',2^8,...
% $$$                   'Fs',  lfp.sampleRate,...
% $$$                   'WinLength',2^7,...
% $$$                   'nOverlap',2^7*0.875,...
% $$$                   'NW',3,...                  'Detrend',[],...
% $$$                   'nTapers',[],...
% $$$                   'FreqRange',[30,200]);
% $$$    
% $$$ [gys,gfs,gts] = fet_spec(Trial,lfp,[],[],[],specArgsTheta);

% $$$ 
% $$$ figure,
% $$$ for c = 1:5
% $$$ subplot(5,1,c)
% $$$ imagesc(gts,gfs,log10(gys(:,:,c))')
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ caxis([2.25,4]);
% $$$ end
% $$$ ForAllSubplots('xlim([4250,4250+5]);')

cmapLimsTheta = [0.8,2.6;0.8,2.6;1,3.2;1,3.5;1,3.5];
cmapLimsTheta = [3.5,6;...
                 3.5,6;...
                 4.5,6.7;...
                 4.5,6.7;...
                 4.5,6.7];
states = {'theta','sit','groom','lpause','lloc','hpause','hloc','rear'};
stateColors = 'kymbbggr';

%% START FIG

[hfig,fig,fax,sax] = set_figure_layout(figure(666008),'1080p','landscape',[],1100,50,10,10);
globalXOffset = 0;
globalYOffset = 0;


    
iax = gobjects([0,1]);
vline = gobjects([0,1]);
timeWindow = 12;
step = 0.008;
markers = {'spine_lower','pelvis_root','spine_middle','spine_upper','head_back','head_left','head_front','head_right'};

for ii = 1:numel(xyz.model.Connections)
    xyz.model.Connections{ii}.color = [0,0,1];
    fxyz.model.Connections{ii}.color = [0,0,1];
end

for ii = 1:numel(xyz.model.Markers)
    xyz.model.Markers{ii}.color = xyz.model.Markers{ii}.color./255;
    fxyz.model.Markers{ii}.color = fxyz.model.Markers{ii}.color./255;
end

%%%<<< PLOT unit raster timeseries
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(8, 1, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','pixels',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width,                        ...
                              fig.subplot.height*10],                    ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
spkRasterHax = sax(end);
unitSet = [fliplr(unitsInt),unitsPyr];
unitClr = [unitsIntColor;unitsPyrColor];
sWidth = 0.4./lfp.sampleRate;
xlim(sax(end),[0,lts(end)]);
ylim(sax(end),[-2*pi,numel(unitSet)+1]);
sphz = copy(phz);
sphz.data(sphz.data<0) = sphz.data(sphz.data<0)+2*pi;
plot([1:size(phz,1)]./phz.sampleRate,sphz(:)-2*pi);
Lines([],1:numel(unitSet),'w');
spkRasterHax.Color = [1,1,1];
vline{end+1} = line(sax(end),[timeWindow/2].*[1,1],ylim());
[yind, yOffSet, xind, xOffSet] = deal(8, 1, 2, 0);

[yind, yOffSet, xind, xOffSet] = deal(8, 1, 2, 0);
% CREATE subplot axes
thetaHax = axes('Units','pixels',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              (fig.page.ypos(yind)+yOffSet+globalYOffset)+fig.subplot.height*10*((2*pi+1)/(numel(unitsInt)+numel(unitsPyr)+2*pi)),...
                              20,                        ...
                              fig.subplot.height*10.*(numel(unitsInt)/(numel(unitsInt)+numel(unitsPyr)+2*pi))],...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
plot(cos(linspace(-pi,pi,100)),linspace(-pi,pi,100));
thetaHax.Visible = 'off';
%thetaHax.Visible = 'on';
thetaHax.YAxisLocation = 'right';
axis(thetaHax,'tight');

% Connect Interneuron raster to theta phase

delete(findobj(fax,'Type','line'))
for intInd = 1:numel(unitsInt)
    line(fax,                                                                        ... axes
         [sum(spkRasterHax.Position([1,3])),                                         ... x1 right edge of spike raster
          sum(spkRasterHax.Position([1,3])) +                                        ... x2 right edge of spike raster
          fig.subplot.horizontalPadding +                                            ... x2 horizontalPadding
          abs(((cos(mpv(intInd)-pi)+1)./2).*thetaHax.Position(3))],                  ... x2 cos(phz) to pixels
         [(fig.page.ypos(yind)+yOffSet+globalYOffset) +                              ... y1 right edge of spike raster
          fig.subplot.height*10*((numel(unitsInt)+1-intInd+2*pi)/(numel(unitsInt)+numel(unitsPyr)+2*pi)),... y1 bot offset from spike raster
          thetaHax.Position(2) +                                                     ... y2 bottom edge of thetaHax
          thetaHax.Position(4)*(mpv(intInd)./diff(ylim(thetaHax)))],         ... y2 offset from bottom spike raster
         'Color','k');
end

%%%<<< PLOT states
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(10, 1, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','Pixels',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width,                        ...
                              fig.subplot.height*2],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
plotSTC(Trial.stc,1,'text',states,stateColors);
ylim(sax(end),[1,9]);
sax(end).XTickLabel = {};
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/7+1,'r');
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/7+1,'g');
box(sax(end),'on');
sax(end).Color = [0.9,0.9,0.9];
vline{end+1} = line(sax(end),[timeWindow/2].*[1,1],ylim());

%%%<<< PLOT ori layer timeseries of low freq spectra
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(11, 1, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','pixels',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width,                        ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
iax{end+1} = imagesc();
axis(sax(end),'xy');
caxis(sax(end),cmapLimsTheta(1,:));
colormap(sax(end),'jet');
axis(sax(end),'tight');
sax(end).XTickLabel = {};
box(sax(end),'on');
ylabel({'ori','Hz'});
ylim(lfs([1,end]));
vline{end+1} = line(sax(end),[timeWindow/2].*[1,1],ylim());

%%%<<< PLOT pyramidal layer timeseries of low freq spectra
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(12, 1, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','pixels',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width,                        ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
iax{end+1} = imagesc();
axis(sax(end),'xy');
caxis(sax(end),cmapLimsTheta(2,:));
colormap(sax(end),'jet');
axis(sax(end),'tight');
sax(end).XTickLabel = {};
box(sax(end),'on');
ylabel({'pyr','Hz'});
ylim(lfs([1,end]));
vline{end+1} = line(sax(end),[timeWindow/2].*[1,1],ylim());

%%%<<< PLOT RAD layer timeseries of low freq spectra
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(13, 1, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','pixels',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width,                        ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
iax{end+1} = imagesc();
axis(sax(end),'xy');
caxis(sax(end),cmapLimsTheta(3,:));
colormap(sax(end),'jet');
axis(sax(end),'tight');
sax(end).XTickLabel = {};
box(sax(end),'on');
ylabel({'rad','Hz'});
ylim(lfs([1,end]));
vline{end+1} = line(sax(end),[timeWindow/2].*[1,1],ylim());

%%%<<< PLOT LM layer timeseries of low freq spectra
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(14, 1, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','pixels',                                     ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width,                        ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
iax{end+1} = imagesc();
axis(sax(end),'xy');
caxis(sax(end),cmapLimsTheta(4,:));
colormap(sax(end),'jet');
axis(sax(end),'tight');
sax(end).XTickLabel = {};
box(sax(end),'on');
ylabel({'lm','Hz'});
ylim(lfs([1,end]));
vline{end+1} = line(sax(end),[timeWindow/2].*[1,1],ylim());

%%%<<< PLOT LM layer timeseries of low freq spectra
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(15, 1, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','pixels',                                     ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width,                        ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
iax{end+1} = imagesc();
axis(sax(end),'xy');
caxis(sax(end),cmapLimsTheta(4,:));
colormap(sax(end),'jet');
axis(sax(end),'tight');
xlabel(sax(end),'Seconds');
box(sax(end),'on');
ylabel({'dg','Hz'});
ylim(lfs([1,end]));
vline{end+1} = line(sax(end),[timeWindow/2].*[1,1],ylim());
%% LINK spectra and raster
linkaxes(sax,'x');
%%%>>>




%%%<<< PLOT 3D stick-skelleton behaving rat
[yind, yOffSet, xind, xOffSet] = deal(9, 20+fig.subplot.height.*5, 2, 80);
%[yind, yOffSet, xind, xOffSet] = deal(9, 20+fig.subplot.height.*5, 2, 40);
% CREATE subplot axes
sax(end+1) = axes('Units','pixels',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width./5,                        ...
                              fig.subplot.height.*5],                    ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
skeletonHax = sax(end);
grid(skeletonHax,'on');
%%%>>>

%%%<<< PLOT behaving rat maze view
%[yind, yOffSet, xind, xOffSet] = deal(11, 20-fig.subplot.height.*5, 2, 150+ fig.subplot.width./5);
[yind, yOffSet, xind, xOffSet] = deal(9, 20+fig.subplot.height.*5, 2, 190+ fig.subplot.width./5);
%[yind, yOffSet, xind, xOffSet] = deal(9, 20+fig.subplot.height.*5, 2, 150+ fig.subplot.width./5);
% CREATE subplot axes
sax(end+1) = axes('Units','pixels',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width./5,                        ...
                              fig.subplot.height.*5],                    ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
mazeHax = sax(end);
grid(mazeHax,'on');
% inner circle
patch(mazeHax,                   ... a
      175.*cos(linspace(-pi,pi,100)), ... x
      175.*sin(linspace(-pi,pi,100)), ... y
      -1.*ones([1,100]),              ... z
      'FaceColor','k',                ... c
      'EdgeColor','none',             ...
      'FaceAlpha',0.3);
abins = linspace(-pi,pi,9);
% outer sectors
for p = 1:8
patch(mazeHax,                   ... a
      [175.*cos(linspace(abins(p),abins(p+1),100)),fliplr(450.*cos(linspace(abins(p),abins(p+1),100)))], ... x
      [175.*sin(linspace(abins(p),abins(p+1),100)),fliplr(450.*sin(linspace(abins(p),abins(p+1),100)))], ... y
      -1.*ones([1,200]),              ... z
      'FaceColor',ucolors(p,:),                ... c
      'EdgeColor','none',             ...
      'FaceAlpha',0.3);
end
halfSpkWindow = 0.300; % Seconds
[posterior,unitCount] = compute_bayesian_posterior_ufr_boxcar(Trial,1000,sampleRate,ufr,ratemaps,halfSpkWindow);
posterior(~mask(:)) = nan;
posterior(posterior<0.001) = nan;
paxOH = pcolor(mazeHax,                                            ...
               pfs{1}.adata.bins{1}-pfs{1}.parameters.binDims(1)/2,...
               pfs{1}.adata.bins{2}-pfs{1}.parameters.binDims(2)/2,...
               reshape(posterior,pfs{1}.adata.binSizes')');
paxOH.EdgeColor = 'none';
caxis(mazeHax,[1.0e-5,0.2]);
axis(mazeHax,'xy');
daspect(mazeHax,[1,1,1]);
%%%>>>


%%%<<< PLOT behaving rate with bayes decoded position (300ms)
[yind, yOffSet, xind, xOffSet] = deal(15, 60+fig.subplot.height.*5, 2, 80);
% CREATE subplot axes
sax(end+1) = axes('Units','pixels',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width./5,                        ...
                              fig.subplot.height.*5],                    ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
mazeBayesHax = sax(end);
% inner circle
patch(mazeBayesHax,                   ... a
      175.*cos(linspace(-pi,pi,100)), ... x
      175.*sin(linspace(-pi,pi,100)), ... y
      -1.*ones([1,100]),              ... z
      'FaceColor','k',                ... c
      'EdgeColor','none',             ...
      'FaceAlpha',0.3);
abins = linspace(-pi,pi,9);
% outer sectors
for p = 1:8
patch(mazeBayesHax,                   ... a
      [175.*cos(linspace(abins(p),abins(p+1),100)),fliplr(450.*cos(linspace(abins(p),abins(p+1),100)))], ... x
      [175.*sin(linspace(abins(p),abins(p+1),100)),fliplr(450.*sin(linspace(abins(p),abins(p+1),100)))], ... y
      -1.*ones([1,200]),              ... z
      'FaceColor',ucolors(p,:),                ... c
      'EdgeColor','none',             ...
      'FaceAlpha',0.3);
end
halfSpkWindow = 0.300; % Seconds
[posterior,unitCount] = compute_bayesian_posterior_ufr_boxcar(Trial,1000,sampleRate,ufr,ratemaps,halfSpkWindow);
posterior(~mask(:)) = nan;
posterior(posterior<0.001) = nan;
pax = pcolor(pfs{1}.adata.bins{1}-pfs{1}.parameters.binDims(1)/2,...
             pfs{1}.adata.bins{2}-pfs{1}.parameters.binDims(2)/2,...
             reshape(posterior,pfs{1}.adata.binSizes')');
pax.EdgeColor = 'none';
caxis(mazeBayesHax,[1.0e-5,0.2]);
axis(mazeBayesHax,'xy');
xlim(sax(end),[-500,500]);
ylim(sax(end),[-500,500]);
sax(end).XTick = [-400:200:400];
sax(end).YTick = [-400:200:400];
sax(end).XTickLabel = [];
sax(end).YTickLabel = [];
hold(sax(end),'on');
daspect(sax(end),[1,1,1]);
title(mazeBayesHax,'Bayesian Decoding: Window 300ms');
%%%>>>


%%%<<< PLOT behaving rat with bayes decoded position (40ms)
[yind, yOffSet, xind, xOffSet] = deal(15, 60+fig.subplot.height.*5, 2, 190+ fig.subplot.width./5);
% CREATE subplot axes
sax(end+1) = axes('Units','pixels',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width./5,                        ...
                              fig.subplot.height.*5],                    ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
mazeBayesFineHax = sax(end);
grid(sax(end),'on');
% inner circle
patch(mazeBayesFineHax,                   ... a
      175.*cos(linspace(-pi,pi,100)), ... x
      175.*sin(linspace(-pi,pi,100)), ... y
      -1.*ones([1,100]),              ... z
      'FaceColor','k',                ... c
      'EdgeColor','none',             ...
      'FaceAlpha',0.3);
abins = linspace(-pi,pi,9);
% outer sectors
for p = 1:8
patch(mazeBayesFineHax,                   ... a
      [175.*cos(linspace(abins(p),abins(p+1),100)),fliplr(450.*cos(linspace(abins(p),abins(p+1),100)))], ... x
      [175.*sin(linspace(abins(p),abins(p+1),100)),fliplr(450.*sin(linspace(abins(p),abins(p+1),100)))], ... y
      -1.*ones([1,200]),              ... z
      'FaceColor',ucolors(p,:),                ... c
      'EdgeColor','none',             ...
      'FaceAlpha',0.3);
end
halfSpkWindow = 0.04; % Seconds
[posterior,unitCount] = compute_bayesian_posterior_ufr_boxcar(Trial,1000,sampleRate,ufr,ratemaps,halfSpkWindow);
posterior(~mask(:)) = nan;
posterior(posterior<0.001) = nan;
paxFine = pcolor(mazeBayesFineHax,...
                 pfs{1}.adata.bins{1}-pfs{1}.parameters.binDims(1)/2,...
             pfs{1}.adata.bins{2}-pfs{1}.parameters.binDims(2)/2,...
             reshape(posterior,pfs{1}.adata.binSizes')');
paxFine.EdgeColor = 'none';
caxis(mazeBayesFineHax,[1.0e-5,0.2]);
axis(mazeBayesFineHax,'xy');
xlim(sax(end),[-500,500]);
ylim(sax(end),[-500,500]);
sax(end).XTick = [-400:200:400];
sax(end).YTick = [-400:200:400];
sax(end).XTickLabel = [];
sax(end).YTickLabel = [];
hold(sax(end),'on');
daspect(sax(end),[1,1,1]);
title(mazeBayesFineHax,'Bayesian Decoding: Window 40ms');
%%%>>>


%%%<<< PLOT behavior decoding 300ms
[yind, yOffSet, xind, xOffSet] = deal(15, 280-fig.subplot.height.*5, 2, 80);
sax(end+1) = axes('Units','pixels',                                     ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width./5,                     ...
                              fig.subplot.height.*5],                   ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
bhvHax = sax(end);
grid(bhvHax,'on');
halfSpkWindow = 0.300; % Seconds
[posterior,unitCount] = compute_bayesian_posterior_ufr_boxcar(Trial,1500,sampleRate,ufr,bhvRatemaps,halfSpkWindow);
posterior(~validDims(:)) = nan;
posterior(posterior<0.001) = nan;
bhvMask = double(validDims);
pchGrid = cell([2,1]);
[pchGrid{:}] = ndgrid(bfs{1}.adata.bins{:});
paxBhvBkgr = contour(bhvHax,pchGrid{:},reshape(bhvMask,bfs{1}.adata.binSizes'),[0.5,0.5],'-k');
paxBhv = pcolor(bhvHax,                                             ... Axes handle
                bfs{1}.adata.bins{1}-bfs{1}.parameters.binDims(1)/2,... x bins
                bfs{1}.adata.bins{2}-bfs{1}.parameters.binDims(2)/2,... y bins
                reshape(posterior,bfs{1}.adata.binSizes')');          % bhv posterior
paxBhv.EdgeColor = 'none';
xlim([-1.8,0.5]);
ylim([-0.5,1.8]);
pitchMarker = scatter3(bhvHax,pch(1500,1),pch(1500,2),1,20,'r','Filled');
xlabel(bhvHax,'Head-Body Pitch');
ylabel(bhvHax,'Body Pitch');
daspect(bhvHax,[1,1,1]);
%%%>>>


%%%<<< PLOT behavior decoding 40ms
[yind, yOffSet, xind, xOffSet] = deal(15, 280-fig.subplot.height.*5, 2, 190+ fig.subplot.width./5);
sax(end+1) = axes('Units','pixels',                                     ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width./5,                     ...
                              fig.subplot.height.*5],                   ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
bhvFineHax = sax(end);
grid(bhvFineHax,'on');
halfSpkWindow = 0.040; % Seconds
[posterior,unitCount] = compute_bayesian_posterior_ufr_boxcar(Trial,1500,sampleRate,ufr,bhvRatemaps,halfSpkWindow);
posterior(~validDims(:)) = nan;
posterior(posterior<0.001) = nan;
bhvMask = double(validDims);
pchGrid = cell([2,1]);
[pchGrid{:}] = ndgrid(bfs{1}.adata.bins{:});
paxBhvFineBkgr = contour(bhvFineHax,pchGrid{:},reshape(bhvMask,bfs{1}.adata.binSizes'),[0.5,0.5],'-k');
paxFineBhv = pcolor(bhvFineHax,                                             ... Axes handle
                bfs{1}.adata.bins{1}-bfs{1}.parameters.binDims(1)/2,... x bins
                bfs{1}.adata.bins{2}-bfs{1}.parameters.binDims(2)/2,... y bins
                reshape(posterior,bfs{1}.adata.binSizes')');          % bhv posterior
paxFineBhv.EdgeColor = 'none';
xlim([-1.8,0.5]);
ylim([-0.5,1.8]);
xlabel(bhvFineHax,'Head-Body Pitch');
ylabel(bhvFineHax,'Body Pitch');
daspect(bhvFineHax,[1,1,1]);
pitchFineMarker = scatter3(bhvFineHax,pch(1500,1),pch(1500,2),1,20,'r','Filled');
%%%>>>


% RECORD video while incrementing time
% $$$ vidObj = VideoWriter('/storage/share/Projects/BehaviorPlaceCode/explore_spk_lfp_state','Archival');
%vidPicDir = create_directory('/storage/share/Projects/BehaviorPlaceCode/explore_spk_lfp_state_some');
% $$$ vidObj.Height = fig.page.height;
% $$$ vidObj.Width = fig.page.width;

% $$$ open(vidObj);
%start = -timeWindow/2+0.2;
start = 3805;
x = start;
% $$$ x = 4304-timeWindow-5;
% $$$ x = 4230;
fstart = x;
fend = x+timeWindow;
xlim(spkRasterHax,[fstart,fend]);
lh = gobjects([0,1]);
for u = unitSet
    uind = find(u==unitSet);
    res = spk(u);
    res = res(res./spk.sampleRate>(fstart)&res./spk.sampleRate<fend);    
    %res = res(res./spk.sampleRate<timeWindow);
    if ~isempty(res)
        for r = 1:numel(res)
            lh(end+1) = line(spkRasterHax,repmat(res(r),[1,2])./spk.sampleRate,...
                             [uind,uind+1],'Color',unitClr(uind,:),'LineWidth',1);
        end
    end
end

figFrame = getframe(hfig);
figFrame = repmat(figFrame,[50,1]);

segment = 1;
frame = 1;

% $$$ for x = 0.1:step:10%round(lts(end))  
for x = start:step:start+300%round(lts(end))      
    tic
    fstart = x;
    fend = x+timeWindow;
    xlim(spkRasterHax,[fstart,fend]);
    oldLH = arrayfun(@(l) l.XData(1)<fstart,lh);
    delete(lh(oldLH));
    lh(oldLH) = [];
    for u = unitSet
        uind = find(u==unitSet);
        res = spk(u);
        res = res(res./spk.sampleRate>(fend-step) & res./spk.sampleRate<fend);
        if ~isempty(res)
            for r = 1:numel(res)
            lh(end+1) = line(spkRasterHax,repmat(res(r),[1,2])./spk.sampleRate,...
                 [uind,uind+1],'Color',unitClr(uind,:),'LineWidth',1);
            end
        end
    end

    tindex = round((x+timeWindow/2).*xyz.sampleRate);
    
% UPDATE side perspective
    [mkrs,stks,mkrc] = plotSkeletonLine(Trial,fxyz,tindex,'line',ang,markers,[],[],skeletonHax);
    view(skeletonHax,circ_rad2ang(ang(tindex,'spine_lower','spine_upper',1)+pi/4),20);
    [xp,yp] = deal(mkrs{4}.XData(1),mkrs{4}.YData(1));
    xlim(skeletonHax,xp(1)+[-180,180]);
    ylim(skeletonHax,yp(1)+[-180,180]);
    zlim(skeletonHax,[0,300]);      
    
% UPDATE top perspective fixed
    index = round((x+timeWindow/2)*ufr.sampleRate);
    [posterior,unitCount] = compute_bayesian_posterior_ufr_boxcar(Trial,index,sampleRate,ufr,ratemaps,0.15);
    posterior(~mask(:)) = nan;
    posterior(posterior<prctile(posterior,95)) = nan;
    caxis(mazeHax,[min(posterior),max(posterior)]);
    if ~isempty(posterior)
        paxOH.CData = reshape(posterior,pfs{1}.adata.binSizes')';
        paxOH.ZData = ones(size(paxOH.ZData));        
    else
        paxOH.CData = nan(pfs{1}.adata.binSizes');
    end    
    plotSkeletonLine(Trial,xyz,tindex,'line',ang,markers,[],[],mazeHax);
    xlim(mazeHax,xp(1)+[-150,150]);
    ylim(mazeHax,yp(1)+[-150,150]);
    view(0,90);
    
    
% UPDATE spectra
    trange = xlim(iax{1}.Parent);
    trange = round(trange(1).*lys.sampleRate):round(trange(2).*lys.sampleRate);
    trange(trange<1) = 1;
    for cc = 1:5,
        iax{cc}.XData = lts(trange);
        iax{cc}.YData = lfs;
        iax{cc}.CData = log10(lys(trange,:,cc))';
    end

% UPDATE position decoding overhead
    index = round((x+timeWindow/2)*ufr.sampleRate);
    [posterior,unitCount] = compute_bayesian_posterior_ufr_boxcar(Trial,index,sampleRate,ufr,ratemaps,0.15);
    posterior(~mask(:)) = nan;
    posterior(posterior<prctile(posterior,95)) = nan;
    caxis(mazeBayesHax,[min(posterior),max(posterior)]);
    if ~isempty(posterior)
        pax.CData = reshape(posterior,pfs{1}.adata.binSizes')';
        pax.ZData = ones(size(pax.ZData));        
        pax.Visible = 'on';        
    else
        pax.CData = nan(pfs{1}.adata.binSizes');
        pax.Visible = 'off';
    end    
    plotSkeletonLine(Trial,fxyz,index,'line',ang,markers,[],[],mazeBayesHax);

% UPDATE position decoding overhead
    index = round((x+timeWindow/2)*ufr.sampleRate);
    [posterior,unitCount] = compute_bayesian_posterior_ufr_boxcar(Trial,index,sampleRate,ufr,ratemaps,0.02);
    posterior(~mask(:)) = nan;
    posterior(posterior<prctile(posterior,95)) = nan;
    caxis(mazeBayesFineHax,[min(posterior),max(posterior)]);
    if ~isempty(posterior)
        paxFine.CData = reshape(posterior,pfs{1}.adata.binSizes')';
        paxFine.ZData = ones(size(paxFine.ZData));        
        paxFine.Visible = 'on';
    else
        paxFine.CData = nan(pfs{1}.adata.binSizes');
        paxFine.Visible = 'off';        
    end    
    plotSkeletonLine(Trial,fxyz,index,'line',ang,markers,[],[],mazeBayesFineHax);

% UPDATE Behavior Decoding
    [posterior,unitCount] = compute_bayesian_posterior_ufr_boxcar(Trial,index,sampleRate,ufr,bhvRatemaps,0.15);
    posterior(~validDims(:)) = nan;
    caxis(bhvHax,[min(posterior),max(posterior)]);
    if ~isempty(posterior)
        paxBhv.CData = reshape(posterior,bfs{1}.adata.binSizes')';
        paxBhv.ZData = ones(size(paxBhv.ZData));        
        paxBhv.Visible = 'on';
    else
        paxBhv.CData = nan(pfs{1}.adata.binSizes');
        paxBhv.Visible = 'off';
    end    
    pitchMarker.XData = pch(index,1);
    pitchMarker.YData = pch(index,2);    


% UPDATE Fine Behavior Decoding
    [posterior,unitCount] = compute_bayesian_posterior_ufr_boxcar(Trial,index,sampleRate,ufr,bhvRatemaps,0.02);
    posterior(~validDims(:)) = nan;
    caxis(bhvFineHax,[min(posterior),max(posterior)]);
    if ~isempty(posterior)
        paxFineBhv.CData = reshape(posterior,bfs{1}.adata.binSizes')';
        paxFineBhv.ZData = ones(size(paxFineBhv.ZData));        
        paxFineBhv.Visible = 'on';
    else
        paxFineBhv.CData = nan(pfs{1}.adata.binSizes');
        paxFineBhv.Visible = 'off';        
    end    
    pitchFineMarker.XData = pch(index,1);
    pitchFineMarker.YData = pch(index,2);    
    
% SET vertical line time position
    set([vline{:}],'XData',[1,1].*(x+timeWindow/2));
    
    drawnow();

    figFrame(frame) = getframe(hfig);
    %imwrite(frame2im(getframe(hfig)),fullfile(vidPicDir,['explore_spk_lfp_',num2str(frame,'%08d'),'.tiff']),'tiff');    

    delete(findobj(skeletonHax,'Type','line'));
    delete(findobj(mazeHax,'Type','line'));    
    delete(findobj(mazeBayesHax,'Type','line'));
    delete(findobj(mazeBayesFineHax,'Type','line'));    
    
    disp(num2str([frame,toc],'[INFO] frame %d ... %d'))
    frame = frame+1;

    if frame==51
        tic
        cd(vidPicDir);
        for f = 1:50,
            imwrite(frame2im(figFrame(f)),fullfile(vidPicDir,['explore_spk_lfp_',num2str(f,'%08d'),'.tiff']),'tiff');
        end

        disp(num2str([toc],'[INFO] writing images  ... %d'));
        tic
        system(['ffmpeg -framerate 30 -i ''explore_spk_lfp_%08d.tiff''  rest_jg05-20120312_seg',num2str(segment,'%08d'),'.mp4 && rm ./*.tiff']);
        disp(num2str([toc],'[INFO] writing video segment  ... %d'));
        segment = segment+1;
        frame = 1;
    end
    
end


system('ffmpeg -safe 0 -f concat -i <(find . -type f -name ''*'' -printf "file ''$PWD/%p''\n" | sort) -c copy active_example.mp4')
!mv active_example.mp4 ../videos/
system(['rm ./*_seg*.mp4']);

%ffmpeg -framerate 30 -i 'explore_spk_lfp_%08d.tiff'  rest_jg05-20120312.mp4


% $$$ close(vidObj);
% $$$ 
% $$$ figure()
% $$$ hold('on');
% $$$ plot(mean(log10(lys(stc{'s'},lfs>5&lfs<10,3)),2)./mean(log10(lys(stc{'s'},lfs<4|(lfs>12&lfs<15),3)),2),...
% $$$      mean(log10(lys(stc{'s'},lfs>5&lfs<10,4)),2)./mean(log10(lys(stc{'s'},lfs<4|(lfs>12&lfs<15),4)),2),...,...
% $$$      '.');
% $$$ hold('on');
% $$$ plot(mean(log10(lys(stc{'w'},lfs>5&lfs<10,3)),2)./mean(log10(lys(stc{'w'},lfs<4|(lfs>12&lfs<15),3)),2),...
% $$$      mean(log10(lys(stc{'w'},lfs>5&lfs<10,4)),2)./mean(log10(lys(stc{'w'},lfs<4|(lfs>12&lfs<15),4)),2),...,...
% $$$      '.g');
% $$$ caxisLims = cat(1,xlim(),ylim());
% $$$ 
% $$$ 
% $$$ figure()
% $$$ hold('on');
% $$$ plot(mean(log10(lys(stc{'s'},lfs>5&lfs<10,3)),2)./mean(log10(lys(stc{'s'},lfs<4,2)),2),...
% $$$      mean(log10(lys(stc{'s'},lfs>5&lfs<10,4)),2)./mean(log10(lys(stc{'s'},lfs<4,4)),2),...
% $$$      '.');
% $$$ hold('on');
% $$$ plot(mean(log10(lys(stc{'w'},lfs>5&lfs<10,3)),2)./mean(log10(lys(stc{'w'},lfs<4,2)),2),...
% $$$      mean(log10(lys(stc{'w'},lfs>5&lfs<10,4)),2)./mean(log10(lys(stc{'w'},lfs<4,4)),2),...
% $$$      '.g');
% $$$ caxisLims = cat(1,xlim(),ylim());
% $$$ 
% $$$ 
% $$$ figure()
% $$$ hold('on');
% $$$ plot(mean(log10(lys(stc{'s'},lfs>5&lfs<10,4)),2),...
% $$$      mean(log10(lys(stc{'s'},lfs>5&lfs<10,5)),2),...
% $$$      '.');
% $$$ hold('on');
% $$$ plot(mean(log10(lys(stc{'w'},lfs>5&lfs<10,4)),2),...
% $$$      mean(log10(lys(stc{'w'},lfs>5&lfs<10,5)),2),...
% $$$      '.g');
% $$$ caxisLims = cat(1,xlim(),ylim());
% $$$ line([1,4],[1,4],'Color','k');
% $$$ xlim(caxisLims(1,:));
% $$$ ylim(caxisLims(2,:));
% $$$ 
% $$$ figure()
% $$$ hist2([mean(log10(lys(stc{'w+s'},lfs>5&lfs<10,3)),2),...
% $$$      mean(log10(lys(stc{'w+s'},lfs>5&lfs<10,4)),2)],...
% $$$      linspace(caxisLims(1,1),caxisLims(1,2),100),...
% $$$      linspace(caxisLims(2,1),caxisLims(2,2),100));
% $$$ 
% $$$ 
