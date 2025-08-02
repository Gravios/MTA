% $$$ req20220506
% $$$     Tags:  thetarc interneurons
% $$$     Status: Active
% $$$     Type: Analysis
% $$$     Author: Justin Graboski
% $$$     Final_Forms: NA
% $$$     Project: General
% $$$     Description: theta return current feature relation to interneurons


sessionList = get_session_list_v2('MjgER2016');

Trial = MTATrial.validate(sessionList(20));

unitsInts = cf(@(T)  T.spk.get_unit_set(T,'interneurons'),  Trials); 

spk = Trial.load('spk',1250,'',units{20},'');

xyz = preproc_xyz(Trial,'trb');

% THETA Periods
tper = resample(cast([stc{'t-m-s',1250}],'TimeSeries'),phz);
% REM Periods
rper = resample(cast([stc{'t&s',1250}],'TimeSeries'),phz);

% EXAMPLE Spike Phase Histogram
u = 80;
res =  spk(u);
sphz = phz(res);
figure,
subplot(211);
histcirc(sphz(tper(res)==1))
subplot(212);
histcirc(sphz(rper(res)==1))

%int = Trial.load('spk',1250,'',unitsInts{20},'');

mphz = copy(phz);
mphz.data = [];
dphz = copy(phz);
dphz.data = [];
ires = discretize(spk.res,1:5000:size(phz,1));
uires = unique(ires)';
uires(isnan(uires)) = [];
for ind = uires
    mphz.data(ind) = circ_mean(phz(spk.res(ires == ind)));
end
mphz.data(mphz.data<0) = mphz.data(mphz.data<0)+2*pi;

% LOAD LFP channel
lfp = Trial.load('lfp',Trial.meta.channelGroup.theta);
% LOAD RC channels 
rfp = Trial.load('lfp',Trial.meta.channelGroup.thetarc);
rfp.data = diff(rfp.data,1,2);


tpw = copy(lfp);
tpw.filter('ButFilter',4,[6,11],'bandpass');
tpws = sqrt(sum(tpw.segs(1:5000:size(phz,1),5000).^2))';

dpw = copy(lfp);
dpw.filter('ButFilter',4,[3,5],'bandpass');
dpws = sqrt(sum(dpw.segs(1:5000:size(phz,1),5000).^2))';



specFet = copy(lfp);
specFet.data = [specFet.data,rfp.data];
swOrd = 12;
specArgs = struct('nFFT'     ,2^(swOrd),...
                  'Fs'       ,specFet.sampleRate,...
                  'WinLength',2^(swOrd-1),...
                  'nOverlap' ,2^(swOrd-1)*.875,...
                  'FreqRange',[.5,20]);
[ys,fs,ts] = fet_spec(Trial,specFet,'mtcsdglong',false,[],specArgs,[],true);

vxyz = vel(resample(copy(xyz),ys),{'bcom','hcom'},[1:3]);

figure,
imagesc(ts,fs,log10(sqrt(ys(:,:,1,1).^2))');
axis('xy');
colormap('jet');

figure,
subplot(6,1,[1,2]);
imagesc(ts,fs,log10(sqrt(ys(:,:,1,1).^2))');
axis('xy');
colormap('jet');
caxis([3,6]);
subplot(6,1,[3,4]);
%imagesc(ts,fs,log10(sqrt(ys(:,:,2,2).^2))');
imagesc(ts,fs,angle(ys(:,:,1,2))');
axis('xy');
colormap(gca(),'hsv');
%caxis([3,6]);
caxis([-pi,pi]);
subplot(6,1,[5]);
hold('on');
plot([0:5000:size(phz,1)]./1250,[0;mphz.data']);
plot(ts,vxyz(:,1)./10,'m')
subplot(6,1,[6]);
plotSTC(Trial.stc,1,[],{'theta','walk','rear','turn','pause','groom','sit'},'kbrgcmy');
linkx();

tdrat = copy(ys);
tdrat.data = mean(log10(sqrt(ys(:,fs>6&fs<10,1,1).^2)),2)./mean(log10(sqrt(ys(:,fs<4|(fs>12&fs<14),1,1).^2)),2);

ttrat = copy(ys);
ttrat.data = mean(log10(sqrt(ys(:,fs>6&fs<10,1,1).^2)),2)./mean(log10(sqrt(ys(:,fs>6&fs<10,2,2).^2)),2);

figure();
inds = {':',runPer,remPer};
for s = 1:3
    subplot(3,1,s);
    hist2([angle(ys(inds{s},14,1,2)),tdrat(inds{s},1)],linspace(-pi,pi,64),linspace(0.8,1.8,64))
end
linkx();


figure,
subplot(211);
histcirc(angle(ys([Trial.stc{'t-m-s'}],tfs,1,2))',64);
xlim([0,720]);
subplot(212);
histcirc(angle(ys([Trial.stc{'t&s'}],tfs,1,2))',64);
xlim([0,720]);

figure,
subplot(211);
rose(angle(ys([Trial.stc{'t-m-s'}],tfs,1,2))',64);
subplot(212);
rose(angle(ys([Trial.stc{'t&s'}],tfs,1,2))',64);

runPer = [Trial.stc{'t-m-s',ys.sampleRate}];
remPer = [Trial.stc{'t&s',ys.sampleRate}];

meanPhzDiffRun = [];
meanTTRatRun = [];
meanSpkPhzRun = [];
meanIntPhzRun = [];
tfs = find(fs>8,1,'first');
for per = 1:size(runPer.data,1)
    if log10(mean(vxyz(runPer(per,:),1))) > -0.5 && diff(runPer(per,:))>10    
        meanPhzDiffRun(per) = circ_mean(angle(ys(runPer(per,:),tfs,1,2)));
        meanTTRatRun(per) = mean(ttrat(runPer(per,:),1));
        meanSpkPhzRun(per) = mean(phz(spk.res(WithinRanges(spk.res,runPer(per,:)./ys.sampleRate* ...
                                                      spk.sampleRate))));
        meanIntPhzRun(per) = mean(phz(int.res(WithinRanges(int.res,runPer(per,:)./ys.sampleRate* ...
                                                      int.sampleRate))));
    else 
        meanPhzDiffRun(per) = nan;
        meanTTRatRun(per) = nan;
        meanSpkPhzRun(per) = nan;
        meanIntPhzRun(per) = nan;
    end
end
meanPhzDiffRem = [];
meanSpkPhzRem = [];
meanTTRatRem = [];
meanIntPhzRem = [];
for per = 1:size(remPer.data,1)
    if log10(mean(vxyz(remPer(per,:),1))) < -0.5 && diff(remPer(per,:))>10
        meanPhzDiffRem(per) = circ_mean(angle(ys(remPer(per,:),tfs,1,2)));
        meanTTRatRem(per) = mean(ttrat(remPer(per,:),1));
        meanSpkPhzRem(per) = mean(phz(spk.res(WithinRanges(spk.res,remPer(per,:)./ys.sampleRate*spk.sampleRate))));
            
        meanIntPhzRem(per) = mean(phz(int.res(WithinRanges(int.res,remPer(per,:)./ys.sampleRate*int.sampleRate))));
    else
        meanTTRatRem(per) = nan;        
        meanPhzDiffRem(per) = nan;
        meanSpkPhzRem(per) = nan;
        meanIntPhzRem(per) = nan;        
    end
end


figure();
hold('on');
scatter(meanPhzDiffRun,meanTTRatRun,20,'k','Filled')
scatter(meanPhzDiffRem,meanTTRatRem,20,'r','Filled')

figure();
hold('on');
scatter(meanPhzDiffRun,meanSpkPhzRun,20,'k','Filled')
scatter(meanPhzDiffRem,meanSpkPhzRem,20,'r','Filled')

figure();
hold('on');
scatter(meanTTRatRun,meanIntPhzRun,20,'k','Filled')
scatter(meanTTRatRem,meanIntPhzRem,20,'r','Filled')


figure();
hold('on');
scatter(meanPhzDiffRun,meanIntPhzRun,20,'k','Filled')
scatter(meanPhzDiffRem,meanIntPhzRem,20,'r','Filled')
daspect([1,1,1]);

figure();
hold('on');
scatter(meanSpkPhzRun,meanIntPhzRun,20,'k','Filled')
scatter(meanSpkPhzRem,meanIntPhzRem,20,'r','Filled')
daspect([1,1,1]);

figure,
subplot(211);
histcirc(runMeanPhzDiff,50)
xlim([0,540]);
subplot(212);
histcirc(remMeanPhzDiff,50)
xlim([0,540]);

% Need the full lfp file and theta periods
% load binary


alfp = LoadBinary(Trial.lfp.fpath,...
                  [Trial.meta.channelGroup.thetarc,Trial.meta.channelGroup.theta],...
                  LoadPar(fullfile(Trial.spath, [Trial.name '.xml'])).nChannels)';
alfp = [alfp(:,3),diff(alfp(:,[1,2]),1,2)];

% Theta periods of size greater than 5 seconds mean interneuron phase and


[ays,afs,ats] = mtcsdglong(alfp,2^12,1250,2^11,2^11.*0.875,[],[],[],[0.5,20]);


attrat = mean(log10(sqrt(ays(:,afs>6&afs<10,1,1).^2)),2)./mean(log10(sqrt(ays(:,afs>6&afs<10,2,2).^2)),2);

hfig = figure();
subplot(411);
imagesc(ats,afs,log10(sqrt(ays(:,:,1,1).^2))');
axis('xy');
caxis([4,6]);
colormap(hfig.CurrentAxes,'jet');
subplot(412);
imagesc(ats,afs,log10(sqrt(ays(:,:,2,2).^2))');
axis('xy');
caxis([4,6]);
colormap(hfig.CurrentAxes,'jet');
subplot(413);
imagesc(ats,afs,angle(ays(:,:,1,2))');
axis('xy');
caxis([-pi,pi]);
colormap(hfig.CurrentAxes,'hsv');
subplot(414);
plot(ats,attrat);
Lines([],1,'k');
ylim([0.8,1.2])
linkx();


sampleRate = 250;

rsIndex = round((1:round(size(alfp,1)./1250.*sampleRate))./sampleRate.*1250);
rsIndex(rsIndex>size(alfp,1)) = [];
rslfp = ButFilter(alfp,4,[sampleRate/2.1]./(1250.*0.5),'low');
rslfp = rslfp(rsIndex,:);

[Res, Clu, Map] = LoadCluRes(fullfile(Trial.spath, Trial.name));

ResPyr = round(Res(ismember(Clu,units{20}))./Trial.sampleRate.*sampleRate);
CluPyr = Clu(ismember(Clu,units{20}));
CluPyr(ResPyr<1|ResPyr>size(rslfp,1)) = [];
ResPyr(ResPyr<1|ResPyr>size(rslfp,1)) = [];

ResInt = round(Res(ismember(Clu,unitsInts{20}))./Trial.sampleRate.*sampleRate);
CluInt = Clu(ismember(Clu,unitsInts{20}));


Phz = Shilbert(ButFilter(rslfp(:,1),4,[5,11]./(sampleRate.*0.5),'bandpass'));
Phz = angle(Phz);

PhzRc = Shilbert(ButFilter(rslfp(:,2),4,[5,11]./(sampleRate.*0.5),'bandpass'));
PhzRc = angle(PhzRc);



UcntInt = accumarray(ResInt,ones(size(ResInt)),size(Phz),@sum);
figure,plot(UcntInt)

UcntPyr = accumarray(ResPyr,ones(size(ResPyr)),size(Phz),@sum);

figure,plot(UcntPyr)

clear('i');
% SETUP sliding window
halfWindowWidth = 64;
slidingWindow =-halfWindowWidth:halfWindowWidth;

meanPhzCpxInt = complex(zeros(size(Phz)));
for index = [1+halfWindowWidth:size(Phz,1)-halfWindowWidth]
    meanPhz = UcntInt(index+slidingWindow) .* exp(i.*Phz(index+slidingWindow));
    meanPhzCpxInt(index) = sum(meanPhz) ./ sum(UcntInt(index+slidingWindow));
end


meanPhzCpxPyr = complex(zeros(size(Phz)));
for index = [1+halfWindowWidth:size(Phz,1)-halfWindowWidth]
    meanPhz = UcntPyr(index+slidingWindow) .* exp(i.*Phz(index+slidingWindow));
    meanPhzCpxPyr(index) = sum(meanPhz) ./ sum(UcntPyr(index+slidingWindow));
end

figure();
hold('on');
plot(angle(meanPhzCpxInt))

tper = load(fullfile(Trial.spath,[Trial.name '.sts.theta']));
tper = round(tper./1250.*sampleRate);

tper = [Trial.stc{'t-m-s'}];
%tper = [Trial.stc{'l+h&t'}];
%tper = [Trial.stc{'l&t'}];
tper = round((tper.origin + tper.data./1250).*sampleRate);

%remPer = [Trial.stc{'h&t'}];
remPer = [Trial.stc{'t&s'}];
remPer = round((remPer.origin + remPer.data./1250).*sampleRate);

ratUsm = interp1(ats+diff(ats(1:2))/2,attrat,[1:size(rslfp,1)]./sampleRate)';


size(meanPhzCpxInt)


figure();
hold('on');
subplot(321);
hist2([SelectPeriods(Phz,tper,'c',1),SelectPeriods(angle(meanPhzCpxInt),tper,'c',1)],...
      linspace(-pi,pi,64),linspace(-pi,pi,64));
subplot(322);
hist2([SelectPeriods(Phz,remPer,'c',1),SelectPeriods(angle(meanPhzCpxInt),remPer,'c',1)],...
      linspace(-pi,pi,64),linspace(-pi,pi,64));
subplot(323);
hist2([SelectPeriods(PhzRc,tper,'c',1),SelectPeriods(angle(meanPhzCpxInt),tper,'c',1)],...
      linspace(-pi,pi,64),linspace(-pi,pi,64));
subplot(324);
hist2([SelectPeriods(PhzRc,remPer,'c',1),SelectPeriods(angle(meanPhzCpxInt),remPer,'c',1)],...
      linspace(-pi,pi,64),linspace(-pi,pi,64));
subplot(325);
hist2([SelectPeriods(angle(meanPhzCpxPyr),tper,'c',1),...
       SelectPeriods(angle(meanPhzCpxInt),tper,'c',1)],...
      linspace(-pi,pi,64),linspace(-pi,pi,64));
subplot(326);
hist2([SelectPeriods(angle(meanPhzCpxPyr),remPer,'c',1),...
       SelectPeriods(angle(meanPhzCpxInt),remPer,'c',1)],...
      linspace(-pi,pi,64),linspace(-pi,pi,64));
colormap('jet');



figure
normax = 'yprob';
normax = '';
subplot(221);
hist2([SelectPeriods(ratUsm,tper,'c',1),...
       SelectPeriods(angle(meanPhzCpxInt),tper,'c',1)],...
      linspace(0.85,1.05,21),linspace(-pi,pi,21),...
      normax);
Lines([],1.5,'r');
Lines(0.95,[],'r');
subplot(222);
hist2([SelectPeriods(ratUsm,remPer,'c',1),...
       SelectPeriods(angle(meanPhzCpxInt),remPer,'c',1)],...
      linspace(0.85,1.05,21),linspace(-pi,pi,21),...
      normax);
Lines([],1.5,'r');
Lines(0.95,[],'r');
subplot(223);
hist2([SelectPeriods(ratUsm,tper,'c',1),...
       SelectPeriods(circ_dist(angle(meanPhzCpxPyr),0),tper,'c',1)],...
      linspace(0.85,1.05,21),linspace(-pi,pi,21),...
      normax);
subplot(224);
hist2([SelectPeriods(ratUsm,remPer,'c',1),...
       SelectPeriods(circ_dist(angle(meanPhzCpxPyr),0),remPer,'c',1)],...
      linspace(0.85,1.05,21),linspace(-pi,pi,21),...
      normax);


figure
subplot(221);
hist2([SelectPeriods(ratUsm,tper,'c',1),...
       SelectPeriods(circ_dist(PhzRc,Phz+pi),tper,'c',1)],... linspace(0.85,1.05,32),linspace(-pi,pi,32));
subplot(222);
hist2([SelectPeriods(ratUsm,remPer,'c',1),...
       SelectPeriods(circ_dist(PhzRc,Phz+pi),remPer,'c',1)],...
      linspace(0.85,1.05,32),linspace(-pi,pi,32));

figure
subplot(221);
hist2([SelectPeriods(angle(meanPhzCpxInt),tper,'c',1),...      
       SelectPeriods(angle(meanPhzCpxPyr),tper,'c',1)],...      
       linspace(-pi,pi,32),linspace(-pi,pi,32));
subplot(222);
hist2([SelectPeriods(angle(meanPhzCpxInt),remPer,'c',1),...      
       SelectPeriods(angle(meanPhzCpxPyr),remPer,'c',1)],...      
       linspace(-pi,pi,32),linspace(-pi,pi,32));

figure();
subplot(211);
histcirc(SelectPeriods(angle(meanPhzCpxPyr),tper,'c',1),32)
subplot(212);
histcirc(SelectPeriods(angle(meanPhzCpxPyr),remPer,'c',1),32)

figure();
subplot(211);
histogram(SelectPeriods(ratUsm,tper,'c',1),linspace(0.85,1.1,32));
subplot(212);
histogram(SelectPeriods(ratUsm,remPer,'c',1),linspace(0.85,1.1,32));


% Mean value of 

%[p,stats] = circ_wwtest(remMeanPhzDiff(nniz(remMeanPhzDiff')),runMeanPhzDiff(nniz(runMeanPhzDiff')))


figure
subplot(211);
ind = [Trial.stc{'t-m-s'}];
hist2([circ_dist(angle(ys(ind,14,1,2)),-pi/2),ttrat(ind)],linspace(-pi,pi,33),linspace(0.8,1.1,40));
subplot(212);
ind = [Trial.stc{'t&s'}];
hist2([circ_dist(angle(ys(ind,14,1,2)),-pi/2),ttrat(ind)],linspace(-pi,pi,33),linspace(0.8,1.1,40));








figure,
hold('on');
plot(0:5000:size(phz,1),[0;mphz.data']);
plot(0:5000:size(phz,1),[0;dphz.data']);
plot(0:5000:size(phz,1),(log10([tpws])'-4).*10);
plot(0:5000:size(phz,1),(log10([dpws])'-4).*10);



